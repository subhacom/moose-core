# reader.py --- clean-room SBML reader for MOOSE
#
# Reads standard SBML (targeting the curated BioModels corpus) and builds a
# MOOSE chemical model that can be simulated. Reactions whose kinetic law is
# recognized as mass-action / Michaelis-Menten map to native Reac / MMenz;
# everything else (arbitrary rate laws, rate rules, assignment rules) is
# compiled to MOOSE ``Function`` objects that feed the same ODE solver. See
# the package docstring for the overall design.

import os

import libsbml
from moose import _moose

from ..base import ModelLoadError
from . import units
from .common import participants
from .identify import identify
from .mathconv import to_exprtk, UnsupportedMath
from .normalize import normalize
from .report import LoadReport


class SBMLValidationError(Exception):
    pass


class _Context:
    def __init__(self, root, report):
        self.root = root
        self.report = report
        self.compartment = {}       # sid -> CubeMesh
        self.comp_size = {}         # sid -> SBML size value (SBML units)
        self.comp_vol_m3 = {}       # sid -> compartment volume in m^3
        self.species = {}           # sid -> Pool/BufPool
        self.species_info = {}      # sid -> {'comp', 'substance_units', 'boundary', 'constant'}
        self.const_param = {}       # sid -> float (constant global parameters)
        self.variable = {}          # sid -> Pool/BufPool proxy holding a raw scalar value
        self.assign_targets = set() # ids on the LHS of an assignment rule
        self.rate_targets = set()   # ids on the LHS of a rate rule
        self.subst_factor = units.AVOGADRO


# ----------------------------------------------------------------------
# per-Function input bookkeeping
# ----------------------------------------------------------------------
class _Inputs:
    """Collects the ordered list of MOOSE elements feeding one Function's
    x[] inputs, and hands out the exprtk token for each referenced symbol."""

    def __init__(self, ctx):
        self.ctx = ctx
        self.order = []          # list of moose elements (connected via nOut)
        self._index = {}         # id(elem) -> x index

    def _slot(self, elem):
        key = elem.path
        if key not in self._index:
            self._index[key] = len(self.order)
            self.order.append(elem)
        return self._index[key]

    def species_token(self, sid):
        info = self.ctx.species_info[sid]
        k = self._slot(self.ctx.species[sid])
        amount = 'x%d/%r' % (k, self.ctx.subst_factor)  # substance units
        if info['substance_units']:
            return '(%s)' % amount
        size = self.ctx.comp_size[info['comp']]
        return '((%s)/%r)' % (amount, size)  # concentration in SBML units

    def variable_token(self, sid):
        k = self._slot(self.ctx.variable[sid])
        return 'x%d' % k  # raw scalar (parameter-like), no unit scaling


def _resolver(inputs):
    ctx = inputs.ctx

    def resolve(name):
        if name in ctx.species:
            return inputs.species_token(name)
        if name in ctx.const_param:
            return repr(ctx.const_param[name])
        if name in ctx.compartment:
            return repr(ctx.comp_size[name])
        if name in ctx.variable:
            return inputs.variable_token(name)
        raise UnsupportedMath('unresolved symbol %r' % name)

    return resolve


def _make_function(path, expr, inputs):
    fn = _moose.Function(path)
    fn.expr = expr
    for i, elem in enumerate(inputs.order):
        _moose.connect(elem, 'nOut', fn.x[i], 'input')
    return fn


# ----------------------------------------------------------------------
# passes
# ----------------------------------------------------------------------
def _create_compartments(model, ctx):
    for i in range(model.getNumCompartments()):
        comp = model.getCompartment(i)
        sid = comp.getId()
        mesh = _moose.CubeMesh('%s/%s' % (ctx.root.path, sid))
        mesh.volume = units.volume(comp)
        ctx.compartment[sid] = mesh
        ctx.comp_size[sid] = comp.getSize() if comp.isSetSize() else 1.0
        ctx.comp_vol_m3[sid] = mesh.volume


def _create_species(model, ctx):
    for i in range(model.getNumSpecies()):
        sp = model.getSpecies(i)
        sid = sp.getId()
        comp_sid = sp.getCompartment()
        parent = ctx.compartment[comp_sid]
        boundary = sp.getBoundaryCondition()
        constant = sp.getConstant()
        # BufPool if the value is held/assigned (constant, boundary, or set by
        # an assignment rule); Pool if it is integrated (by reactions or a rate
        # rule). A rate rule wins even over boundaryCondition.
        held = (constant or boundary or sid in ctx.assign_targets)
        cls = _moose.BufPool if held and sid not in ctx.rate_targets else _moose.Pool
        pool = cls('%s/%s' % (parent.path, sid))
        pool.nInit = units.species_ninit(sp, model.getCompartment(comp_sid), ctx.subst_factor)
        ctx.species[sid] = pool
        ctx.species_info[sid] = {
            'comp': comp_sid,
            'substance_units': sp.getHasOnlySubstanceUnits(),
            'boundary': boundary,
            'constant': constant,
        }


def _create_parameters(model, ctx):
    """Constant parameters fold to literals; a parameter driven by a rule
    becomes a 'variable' proxy pool (Pool if rate-ruled so it integrates,
    BufPool if assignment-ruled so it is overwritten each step). Proxies live
    inside the first compartment so they fall within a solver's scope."""
    host = next(iter(ctx.compartment.values())).path
    for i in range(model.getNumParameters()):
        prm = model.getParameter(i)
        sid = prm.getId()
        value = prm.getValue() if prm.isSetValue() else 0.0
        if sid in ctx.rate_targets:
            proxy = _moose.Pool('%s/%s' % (host, sid))
            proxy.nInit = value
            ctx.variable[sid] = proxy
        elif sid in ctx.assign_targets:
            proxy = _moose.BufPool('%s/%s' % (host, sid))
            proxy.nInit = value
            ctx.variable[sid] = proxy
        else:
            ctx.const_param[sid] = value


def _classify_rules(model, ctx):
    for i in range(model.getNumRules()):
        rule = model.getRule(i)
        if rule.isAssignment():
            ctx.assign_targets.add(rule.getVariable())
        elif rule.isRate():
            ctx.rate_targets.add(rule.getVariable())


def _reaction_home(reac, ctx):
    """(mesh, compartment-id) the reaction is created in: its first substrate's,
    else first product's compartment."""
    subs, prds, _ = participants(reac)
    refs = subs or prds
    if not refs:
        return None, None
    info = ctx.species_info.get(refs[0][0])
    if info is None:
        return None, None
    return ctx.compartment[info['comp']], info['comp']


def _inv(sid, ctx, substance_scale):
    """Factor relating a species' rate-law value to its MOOSE concentration:
    value = inv * conc. (Inverse of the value->conc scale.)"""
    info = ctx.species_info[sid]
    comp = info['comp']
    if info['substance_units']:            # value is an amount
        return ctx.comp_vol_m3[comp] / substance_scale
    size_scale = ctx.comp_vol_m3[comp] / ctx.comp_size[comp]  # m^3 per size unit
    return size_scale / substance_scale


def _si_massaction(result, reac, home_cid, ctx):
    """Convert value-space Kf/Kb to SI (mM-based) constants."""
    ss = ctx.subst_factor / units.AVOGADRO           # substance_scale (mol/unit)
    pref = ss / ctx.comp_vol_m3[home_cid]            # substance_scale / V_home[m^3]
    subs, prds, _ = participants(reac)
    kf = result['Kf_val'] * pref
    kb = result['Kb_val'] * pref
    for sid, st in subs:
        kf *= _inv(sid, ctx, ss) ** st
    for sid, st in prds:
        kb *= _inv(sid, ctx, ss) ** st
    for c in result.get('catalysts', set()):
        kf *= _inv(c, ctx, ss)
        kb *= _inv(c, ctx, ss)
    return kf, kb


def _si_mmenz(result, subs0, home_cid, ctx):
    """Convert value-space kcat/Km to SI (kcat 1/s, Km in mM)."""
    ss = ctx.subst_factor / units.AVOGADRO
    kcat = result['kcat_val'] * _inv(result['enzyme'], ctx, ss) * ss / ctx.comp_vol_m3[home_cid]
    km = result['Km_val'] / _inv(subs0, ctx, ss)
    return kcat, km


def _wire(mreac, reac, ctx):
    """Connect substrate/product pools to a Reac/MMenz, honoring stoichiometry."""
    subs, prds, _ = participants(reac)
    for field, refs in (('sub', subs), ('prd', prds)):
        for sid, st in refs:
            pool = ctx.species.get(sid)
            if pool is not None:
                for _ in range(int(st)):
                    _moose.connect(mreac, field, pool, 'reac', 'OneToOne')


def _reac_id(reac, index):
    return reac.getIdAttribute() if reac.isSetIdAttribute() else 'reaction_%d' % index


def _create_reactions(model, ctx):
    """Recognize mass-action / Michaelis-Menten kinetics numerically and emit
    native Reac / MMenz; everything else goes to the Function fallback."""
    variable_names = set(ctx.variable)
    fallback = []
    for r in range(model.getNumReactions()):
        reac = model.getReaction(r)
        if reac.getFast():
            ctx.report.unsupported_add('fast reaction %r (no QSS handling)' % reac.getId())
        home, home_cid = _reaction_home(reac, ctx)
        result = None
        if home is not None:
            result = identify(reac, ctx.const_param, ctx.comp_size, variable_names)
        if result is None:
            fallback.append(reac)
        elif result['kind'] == 'massaction':
            mr = _moose.Reac('%s/%s' % (home.path, _reac_id(reac, r)))
            mr.Kf, mr.Kb = _si_massaction(result, reac, home_cid, ctx)
            _wire(mr, reac, ctx)
            # A catalyst (modifier that enters the rate multiplicatively) is
            # wired as substrate AND product so it drives the rate without
            # being consumed.
            for cat in result.get('catalysts', ()):
                pool = ctx.species.get(cat)
                if pool is not None:
                    _moose.connect(mr, 'sub', pool, 'reac', 'OneToOne')
                    _moose.connect(mr, 'prd', pool, 'reac', 'OneToOne')
            ctx.report.reactions_native += 1
        elif result['kind'] == 'mmenz' and result['enzyme'] in ctx.species:
            enzpool = ctx.species[result['enzyme']]
            mm = _moose.MMenz('%s/%s' % (enzpool.path, _reac_id(reac, r)))
            mm.kcat, mm.Km = _si_mmenz(result, reac.getReactant(0).getSpecies(),
                                       home_cid, ctx)
            _moose.connect(enzpool, 'nOut', mm, 'enzDest')
            _wire(mm, reac, ctx)
            ctx.report.reactions_native += 1
        else:
            fallback.append(reac)
    _build_fallback_functions(fallback, ctx)


def _build_fallback_functions(reactions, ctx):
    """General-ODE assembly for reactions with non-native rate laws: per
    reacting species build one Function giving dn/dt = subst_factor *
    Sum_j netstoich_ij * KL_j."""
    contrib = {}  # sid -> list of (netstoich, reaction)
    for reac in reactions:
        subs, prds, _ = participants(reac)
        net = {}
        for sid, st in subs:
            net[sid] = net.get(sid, 0.0) - st
        for sid, st in prds:
            net[sid] = net.get(sid, 0.0) + st
        for sid, ns in net.items():
            if ns == 0.0:
                continue
            info = ctx.species_info.get(sid)
            if info is None or info['boundary'] or info['constant']:
                continue  # not altered by reactions
            contrib.setdefault(sid, []).append((ns, reac))

    for sid, terms in contrib.items():
        inputs = _Inputs(ctx)
        resolve = _resolver(inputs)
        try:
            pieces = []
            for ns, reac in terms:
                kl = reac.getKineticLaw()
                if kl is None or kl.getMath() is None:
                    raise UnsupportedMath('reaction %r has no kinetic law' % reac.getId())
                pieces.append('%r*(%s)' % (ns, to_exprtk(kl.getMath(), resolve)))
        except UnsupportedMath as e:
            ctx.report.unsupported_add(
                'species %r derivative not built: %s' % (sid, e))
            continue
        expr = '%r*(%s)' % (ctx.subst_factor, ' + '.join(pieces))
        pool = ctx.species[sid]
        fn = _make_function('%s/dot' % pool.path, expr, inputs)
        _moose.connect(fn, 'valueOut', pool, 'increment')
        ctx.report.reactions_function += 1


def _create_rules(model, ctx):
    for r in range(model.getNumRules()):
        rule = model.getRule(r)
        if rule.isAlgebraic():
            ctx.report.unsupported_add('algebraic rule (no DAE solver in MOOSE)')
            continue
        var = rule.getVariable()
        math = rule.getMath()
        if math is None:
            continue
        inputs = _Inputs(ctx)
        try:
            body = to_exprtk(math, _resolver(inputs))
        except UnsupportedMath as e:
            ctx.report.unsupported_add('rule for %r: %s' % (var, e))
            continue

        target, scale, _is_species = _rule_target(var, ctx)
        if target is None:
            ctx.report.unsupported_add('rule target %r not found' % var)
            continue

        if rule.isAssignment():
            expr = body if scale == 1.0 else '%r*(%s)' % (scale, body)
            fn = _make_function('%s/assign' % target.path, expr, inputs)
            _moose.connect(fn, 'valueOut', target, 'setN')
            ctx.report.assignment_rules += 1
        elif rule.isRate():
            expr = body if scale == 1.0 else '%r*(%s)' % (scale, body)
            fn = _make_function('%s/rate' % target.path, expr, inputs)
            _moose.connect(fn, 'valueOut', target, 'increment')
            ctx.report.rate_rules += 1


def _rule_target(var, ctx):
    """Resolve a rule variable to (moose element, n-scale factor, is_species).
    The scale converts an SBML value/derivative into MOOSE number units."""
    if var in ctx.species:
        info = ctx.species_info[var]
        scale = ctx.subst_factor
        if not info['substance_units']:
            scale = ctx.subst_factor * ctx.comp_size[info['comp']]
        return ctx.species[var], scale, True
    if var in ctx.variable:
        return ctx.variable[var], 1.0, False
    return None, 1.0, False


def _detect_unsupported(model, ctx):
    if model.getNumEvents() > 0:
        ctx.report.unsupported_add(
            '%d event(s): MOOSE has no discrete-event support' % model.getNumEvents())
    if model.getNumCompartments() > 1:
        ctx.report.unsupported_add(
            '%d compartments: cross-compartment coupling not yet supported '
            '(single-compartment models only)' % model.getNumCompartments())


def _setup_solver(ctx, solver):
    """One Ksolve+Stoich per compartment (see moose-examples
    crossComptSimpleReac). Each pool is then solved with its own compartment's
    volume, which a single shared Stoich got wrong for multi-compartment
    models."""
    if solver == 'ee':
        return
    # Cross-compartment Function coupling is not yet handled and crashes the
    # multi-Stoich solver; leave such models unsolved (flagged in the report).
    if len(ctx.compartment) > 1:
        return
    for comp in ctx.compartment.values():
        if solver == 'gssa':
            ksolve = _moose.Gsolve('%s/ksolve' % comp.path)
        else:
            ksolve = _moose.Ksolve('%s/ksolve' % comp.path)
            # Many BioModels are stiff; LSODA (like RoadRunner's CVODE) handles
            # them where the default explicit RKF45 diverges.
            try:
                ksolve.method = 'lsoda'
            except Exception:
                pass
        stoich = _moose.Stoich('%s/stoich' % comp.path)
        stoich.compartment = comp
        stoich.ksolve = ksolve
        stoich.reacSystemPath = comp.path + '/##'


# ----------------------------------------------------------------------
# public handler
# ----------------------------------------------------------------------
class SBMLHandler:
    extensions = ('.xml', '.sbml')

    def __init__(self):
        self.report = None

    def read(self, filepath, loadpath, solver='gsl', validate=True):
        if not os.path.isfile(filepath):
            raise FileNotFoundError('No such file: %s' % filepath)
        doc = libsbml.readSBML(filepath)
        if doc is None:
            raise ModelLoadError('Empty SBML doc', filepath, loadpath)
        if validate and doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR) > 0:
            errs = [doc.getError(i).getMessage()
                    for i in range(doc.getNumErrors())
                    if doc.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_ERROR]
            raise SBMLValidationError('\n'.join(errs))

        report = LoadReport(filepath=filepath, loadpath=loadpath)
        report.normalized, skipped = normalize(doc)
        report.warnings += skipped

        model = doc.getModel()
        if model is None:
            raise ModelLoadError('Invalid SBML: no model element', filepath, loadpath)
        if model.getNumCompartments() == 0:
            raise ModelLoadError('Model has no compartments', filepath, loadpath)

        root = _moose.Neutral(loadpath)
        ctx = _Context(root, report)
        # SI: n = (SBML amount) * NA * substance_scale, volumes in m^3 (units.py).
        ctx.subst_factor = units.AVOGADRO * units.substance_scale(model)

        _classify_rules(model, ctx)
        _create_compartments(model, ctx)
        _create_species(model, ctx)
        _create_parameters(model, ctx)
        _create_reactions(model, ctx)
        _create_rules(model, ctx)
        _detect_unsupported(model, ctx)
        _setup_solver(ctx, solver)

        self.report = report
        return _moose.element(loadpath)

    def write(self, modelpath, filepath, **options):
        raise NotImplementedError('SBML writing not yet implemented')

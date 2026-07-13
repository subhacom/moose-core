"""Numeric identification of native kinetics hidden in arbitrary rate laws.

Many BioModels write mass-action or Michaelis-Menten kinetics in non-canonical
algebraic forms (expanded polynomials, ``pow`` for stoichiometry, an explicit
compartment factor, reversible ``k1*A*B - k2*C`` collapsed, ...). Rather than
pattern-match the syntax, we evaluate the rate law at many random operating
points and least-squares-fit the mass-action / MM functional form. If the fit
is exact the reaction is that kinetics regardless of how it was written, and we
recover the constants. Anything that does not fit falls back to a Function.

Recognizing native forms matters because native Reac/Enz/MMenz (a) are faster,
(b) support the stochastic (GSSA) solver, and (c) are the ONLY way to couple
reactions across compartments (via fixXreacs); a Function cannot read a pool in
another compartment's solver.
"""

import math

import libsbml
import numpy as np

from . import symbolic

_UNARY = {
    libsbml.AST_FUNCTION_EXP: math.exp,
    libsbml.AST_FUNCTION_LN: lambda x: math.log(x) if x > 0 else -1e30,
    libsbml.AST_FUNCTION_ABS: abs,
    libsbml.AST_FUNCTION_SIN: math.sin,
    libsbml.AST_FUNCTION_COS: math.cos,
    libsbml.AST_FUNCTION_TAN: math.tan,
    libsbml.AST_FUNCTION_TANH: math.tanh,
    libsbml.AST_FUNCTION_FLOOR: math.floor,
    libsbml.AST_FUNCTION_CEILING: math.ceil,
}


def _eval(node, values):
    t = node.getType()
    if t == libsbml.AST_INTEGER:
        return float(node.getInteger())
    if t in (libsbml.AST_REAL, libsbml.AST_REAL_E, libsbml.AST_RATIONAL):
        return node.getValue()
    if t == libsbml.AST_NAME:
        return values.get(node.getName(), 1.0)
    if t == libsbml.AST_NAME_AVOGADRO:
        return 6.022140857e23
    if t == libsbml.AST_NAME_TIME:
        return values.get('t', 0.0)
    if t == libsbml.AST_CONSTANT_PI:
        return math.pi
    if t == libsbml.AST_CONSTANT_E:
        return math.e
    n = node.getNumChildren()
    c = [_eval(node.getChild(i), values) for i in range(n)]
    if t == libsbml.AST_PLUS:
        return sum(c) if c else 0.0
    if t == libsbml.AST_TIMES:
        r = 1.0
        for x in c:
            r *= x
        return r
    if t == libsbml.AST_MINUS:
        return -c[0] if n == 1 else c[0] - c[1]
    if t == libsbml.AST_DIVIDE:
        return c[0] / c[1]
    if t in (libsbml.AST_POWER, libsbml.AST_FUNCTION_POWER):
        return c[0] ** c[1]
    if t == libsbml.AST_FUNCTION_ROOT:
        return c[-1] ** (1.0 / c[0]) if n == 2 else math.sqrt(c[0])
    if t in _UNARY:
        return _UNARY[t](c[0])
    raise ValueError('uneval node type %d' % t)


def referenced_names(node, out):
    if node.getType() == libsbml.AST_NAME:
        out.add(node.getName())
    for i in range(node.getNumChildren()):
        referenced_names(node.getChild(i), out)


def identify(reac, model, const_params, comp_sizes, variable_names, home_volume):
    """Recognize native kinetics in a reaction's rate law.

    Symbolic decomposition (``symbolic.analyze``) is the primary recognizer --
    it recovers the reaction structure and exact constants. The result is then
    numerically self-verified: we reconstruct the rate the native objects would
    produce and sample-check it against the original rate law, so a bad
    extraction falls back to the (validated) Function path rather than
    silently mis-simulating.

    Returns a mapping dict or ``None``:

    * ``{'kind':'massaction', 'Kf', 'Kb', 'catalysts': set()}``
    * ``{'kind':'mmenz', 'kcat', 'Km', 'enzyme'}``
    """
    kl = reac.getKineticLaw()
    if kl is None or kl.getMath() is None:
        return None

    names = set()
    referenced_names(kl.getMath(), names)
    # A rate law referencing a rule-driven (time-varying) symbol is not a
    # constant-coefficient native form; treat it via the Function fallback.
    if names & variable_names:
        return None

    subst = dict(const_params)
    subst.update(comp_sizes)

    # Per-species factor relating a rate-law value to a MOOSE concentration:
    # amount (hasOnlySubstanceUnits) species carry value = conc*volume, so
    # their factor is the compartment volume; concentration species use 1.
    scale = _species_scale(reac, model, comp_sizes)

    mapping = symbolic.analyze(reac, subst, home_volume, scale)
    if mapping is None:
        return None
    if not _verify(reac, mapping, subst, home_volume, scale):
        return None
    return mapping


def _species_scale(reac, model, comp_sizes):
    scale = {}

    def add(sid):
        sp = model.getSpecies(sid)
        if sp is not None and sp.getHasOnlySubstanceUnits():
            scale[sid] = comp_sizes.get(sp.getCompartment(), 1.0)
        else:
            scale[sid] = 1.0

    for i in range(reac.getNumReactants()):
        add(reac.getReactant(i).getSpecies())
    for i in range(reac.getNumProducts()):
        add(reac.getProduct(i).getSpecies())
    for i in range(reac.getNumModifiers()):
        add(reac.getModifier(i).getSpecies())
    return scale


def _verify(reac, mapping, subst, V, scale, n=24, tol=1e-6):
    """Sample the original rate law and the rate the chosen native objects
    would reproduce (in MOOSE concentration = value/scale); accept only if they
    agree everywhere."""
    ast = reac.getKineticLaw().getMath()
    subs = [(reac.getReactant(i).getSpecies(),
             reac.getReactant(i).getStoichiometry() or 1.0)
            for i in range(reac.getNumReactants())]
    prds = [(reac.getProduct(i).getSpecies(),
             reac.getProduct(i).getStoichiometry() or 1.0)
            for i in range(reac.getNumProducts())]
    mods = [reac.getModifier(i).getSpecies() for i in range(reac.getNumModifiers())]
    species = sorted(set([s for s, _ in subs] + [s for s, _ in prds] + mods))

    base = dict(subst)
    kl = reac.getKineticLaw()
    for i in range(kl.getNumParameters()):
        base[kl.getParameter(i).getId()] = kl.getParameter(i).getValue()

    rng = np.random.default_rng(20240713)
    worst = 0.0
    denom = 1e-30
    for _ in range(n):
        v = dict(base)
        for s in species:
            v[s] = rng.uniform(0.05, 5.0)
        try:
            r_orig = _eval(ast, v)
        except (ValueError, OverflowError, ZeroDivisionError):
            return False
        if not np.isfinite(r_orig):
            return False
        # MOOSE concentration of each species = rate-law value / scale.
        def conc(sid):
            return v[sid] / scale.get(sid, 1.0)

        if mapping['kind'] == 'massaction':
            cats = mapping.get('catalysts', set())
            cat = 1.0
            for c in cats:
                cat *= conc(c)
            fwd = cat
            for s, st in subs:
                fwd *= conc(s) ** st
            bwd = cat
            for s, st in prds:
                bwd *= conc(s) ** st
            r_recon = V * (mapping['Kf'] * fwd - mapping['Kb'] * bwd)
        elif mapping['kind'] == 'mmenz':
            e = conc(mapping['enzyme'])
            s = conc(subs[0][0])
            r_recon = V * mapping['kcat'] * e * s / (mapping['Km'] + s)
        else:
            return False
        worst = max(worst, abs(r_orig - r_recon))
        denom = max(denom, abs(r_orig))
    return worst / denom < tol

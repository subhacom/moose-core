"""Symbolic decomposition of SBML rate laws into MOOSE reaction structure.

We build a sympy expression directly from the libsbml AST (no slow string
``sympify``), substitute the numeric parameters/compartment sizes, and inspect
the algebra:

* a rate law that is a polynomial in the species with <=2 monomials is
  mass-action -- the monomial coefficients are the exact forward/backward rate
  constants, and any species appearing in a monomial beyond its stoichiometric
  role is a *catalyst* (mapped in MOOSE as a buffered substrate+product);
* a rate law that is rational in the single substrate with a degree-1
  numerator and denominator is Michaelis-Menten -- we read off kcat and Km.

Symbolic analysis (unlike black-box numeric fitting) recovers the reaction
*structure*: which species are substrates/products/catalysts, and the exact
constants. Callers should still numerically verify the reconstruction (see
identify.verify) before trusting it.
"""

import libsbml
import sympy as sp


def _num(x):
    """Return x as a plain float, or None if it still contains free symbols
    (an unresolved parameter) or isn't real -- in which case the reaction is
    not a constant-coefficient native form and must fall back."""
    x = sp.sympify(x)
    if x.free_symbols:
        return None
    try:
        val = complex(x)
    except (TypeError, ValueError):
        return None
    if abs(val.imag) > 1e-12:
        return None
    return val.real

_FN = {
    libsbml.AST_FUNCTION_EXP: sp.exp,
    libsbml.AST_FUNCTION_LN: sp.log,
    libsbml.AST_FUNCTION_ABS: sp.Abs,
    libsbml.AST_FUNCTION_SIN: sp.sin,
    libsbml.AST_FUNCTION_COS: sp.cos,
    libsbml.AST_FUNCTION_TAN: sp.tan,
    libsbml.AST_FUNCTION_TANH: sp.tanh,
}


def ast_to_sympy(node, syms):
    """Build a sympy expression from a libsbml ASTNode. ``syms`` caches
    Symbol objects by name (created on demand)."""
    t = node.getType()
    if t == libsbml.AST_INTEGER:
        return sp.Integer(node.getInteger())
    if t in (libsbml.AST_REAL, libsbml.AST_REAL_E, libsbml.AST_RATIONAL):
        return sp.Float(node.getValue())
    if t == libsbml.AST_NAME:
        n = node.getName()
        return syms.setdefault(n, sp.Symbol(n))
    if t == libsbml.AST_NAME_AVOGADRO:
        return sp.Float(6.022140857e23)
    if t == libsbml.AST_NAME_TIME:
        return syms.setdefault('t', sp.Symbol('t'))
    if t == libsbml.AST_CONSTANT_PI:
        return sp.pi
    if t == libsbml.AST_CONSTANT_E:
        return sp.E
    n = node.getNumChildren()
    c = [ast_to_sympy(node.getChild(i), syms) for i in range(n)]
    if t == libsbml.AST_PLUS:
        return sp.Add(*c) if c else sp.Integer(0)
    if t == libsbml.AST_TIMES:
        return sp.Mul(*c) if c else sp.Integer(1)
    if t == libsbml.AST_MINUS:
        return -c[0] if n == 1 else c[0] - c[1]
    if t == libsbml.AST_DIVIDE:
        return c[0] / c[1]
    if t in (libsbml.AST_POWER, libsbml.AST_FUNCTION_POWER):
        return c[0] ** c[1]
    if t == libsbml.AST_FUNCTION_ROOT:
        return c[-1] ** (sp.Integer(1) / c[0]) if n == 2 else sp.sqrt(c[0])
    if t in _FN:
        return _FN[t](c[0])
    raise ValueError('unhandled AST node type %d' % t)


def analyze(reac, subst, home_volume, species_scale=None):
    """Return a native-mapping dict or None.

    ``species_scale`` maps a species id to the factor relating its rate-law
    value to a MOOSE concentration: 1 for a concentration species, its
    compartment volume for an amount (hasOnlySubstanceUnits) species (whose
    value is an amount = conc*volume). It is folded into the extracted Kf/Kb so
    amount species map to a native Reac correctly.

    * ``{'kind':'massaction', 'Kf', 'Kb', 'catalysts': set()}`` -- Kf/Kb in
      MOOSE concentration units.
    * ``{'kind':'mmenz', 'kcat', 'Km', 'enzyme'}``.
    """
    if species_scale is None:
        species_scale = {}
    kl = reac.getKineticLaw()
    if kl is None or kl.getMath() is None:
        return None

    syms = {}
    for i in range(kl.getNumParameters()):
        syms[kl.getParameter(i).getId()] = sp.Symbol(kl.getParameter(i).getId())
    try:
        expr = ast_to_sympy(kl.getMath(), syms)
    except (ValueError, RecursionError):
        return None

    smap = {syms[k]: sp.Float(v) for k, v in subst.items() if k in syms}
    for i in range(kl.getNumParameters()):
        p = kl.getParameter(i)
        smap[syms[p.getId()]] = sp.Float(p.getValue())
    expr = expr.xreplace(smap)

    subs = [(reac.getReactant(i).getSpecies(),
             int(reac.getReactant(i).getStoichiometry() or 1))
            for i in range(reac.getNumReactants())]
    prds = [(reac.getProduct(i).getSpecies(),
             int(reac.getProduct(i).getStoichiometry() or 1))
            for i in range(reac.getNumProducts())]
    mods = [reac.getModifier(i).getSpecies() for i in range(reac.getNumModifiers())]

    spset = set([s for s, _ in subs] + [s for s, _ in prds] + mods)
    gens = [syms[s] for s in spset if s in syms and syms[s] in expr.free_symbols]
    if not gens:
        return None  # zero-order / constant -- left to the caller / fallback
    sub_st = {syms[s]: st for s, st in subs if s in syms}
    prd_st = {syms[s]: st for s, st in prds if s in syms}

    ma = _mass_action(expr, gens, sub_st, prd_st, subs, prds, home_volume, species_scale)
    if ma is not None:
        return ma
    return _michaelis_menten(expr, syms, subs, prds, mods, home_volume)


def _mass_action(expr, gens, sub_st, prd_st, subs, prds, V, scale):
    try:
        poly = sp.Poly(expr, *gens)
    except (sp.PolynomialError, sp.GeneratorsError):
        return None
    terms = poly.terms()
    if len(terms) > 2:
        return None
    gidx = {g: i for i, g in enumerate(poly.gens)}
    Kf = Kb = None
    catalysts = set()
    for exps, coeff in terms:
        cf = _num(coeff)
        if cf is None:
            return None  # symbolic coefficient -> not constant mass-action
        ge = {g: exps[gidx[g]] for g in poly.gens}
        # Fold the per-species value->concentration factor into the constant:
        # value^exp = (conc*scale)^exp, so the conc-based rate constant picks up
        # prod(scale^exp) over the monomial.
        factor = 1.0
        for g in poly.gens:
            factor *= scale.get(str(g), 1.0) ** ge.get(g, 0)
        cf_eff = cf * factor
        # forward term: contains every substrate at >= its stoichiometry, coeff>0
        fwd = bool(subs) and all(ge.get(g, 0) >= st for g, st in sub_st.items()) and cf > 0
        bwd = bool(prds) and all(ge.get(g, 0) >= st for g, st in prd_st.items()) and cf < 0
        if fwd:
            Kf = cf_eff
            for g in poly.gens:
                if ge.get(g, 0) > sub_st.get(g, 0):
                    catalysts.add(str(g))
        elif bwd:
            Kb = -cf_eff
            for g in poly.gens:
                if ge.get(g, 0) > prd_st.get(g, 0):
                    catalysts.add(str(g))
        else:
            return None
    if Kf is None and Kb is None:
        return None
    return {'kind': 'massaction',
            'Kf': (Kf or 0.0) / V,
            'Kb': (Kb or 0.0) / V,
            'catalysts': catalysts}


def _michaelis_menten(expr, syms, subs, prds, mods, V):
    # MOOSE MMenz requires exactly the enzyme+single-substrate form and at
    # least one product.
    if len(subs) != 1 or not mods or not prds:
        return None
    S = syms.get(subs[0][0])
    if S is None:
        return None
    num, den = sp.fraction(sp.together(expr))
    try:
        pn, pd = sp.Poly(num, S), sp.Poly(den, S)
    except (sp.PolynomialError, sp.GeneratorsError):
        return None
    if pd.degree() != 1 or pn.degree() != 1:
        return None
    a1, a0 = pd.coeff_monomial(S), pd.coeff_monomial(1)
    b1, b0 = pn.coeff_monomial(S), pn.coeff_monomial(1)
    if b0 != 0 or a1 == 0:
        return None
    Km = _num(a0 / a1)
    if Km is None or Km <= 0:
        return None
    gain = b1 / a1  # = V * kcat * E
    for enz in mods:
        e = syms.get(enz)
        if e is not None and e in gain.free_symbols:
            kcat = _num(gain / e)
            if kcat is not None and kcat > 0:
                return {'kind': 'mmenz', 'enzyme': enz, 'kcat': kcat / V, 'Km': Km}
    return None

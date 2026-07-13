"""Unit conventions between SBML and MOOSE.

MOOSE stores a molecule number ``n`` per pool and derives concentration as
``conc = n / (NA * V)``, where ``V`` is the pool's compartment volume. We map
SBML onto this with one consistent internal scale so that **MOOSE's reported
concentration equals the SBML concentration value directly**:

    V      := the SBML compartment size (in the model's own volume units)
    n      := (SBML amount) * NA
    conc   =  n / (NA * V) = amount / size    == SBML concentration

NA here is just a fixed proportionality constant tying ``n`` to amount; we do
not reinterpret the model into SI. This keeps the deterministic ODE exact and
makes results directly comparable to reference simulators. (Stochastic/GSSA
runs, which need true molecule counts, would additionally fold in the
substance-unit -> mole factor; that refinement is deferred.)
"""

AVOGADRO = 6.022140857e23


def volume(comp):
    """MOOSE compartment volume = SBML compartment size (own units)."""
    return comp.getSize() if comp.isSetSize() else 1.0


def species_ninit(species, comp):
    """Initial ``n`` for a species: (initial amount) * NA, resolving
    initialConcentration to an amount via the compartment size."""
    if species.isSetInitialAmount():
        amount = species.getInitialAmount()
    elif species.isSetInitialConcentration():
        size = comp.getSize() if comp.isSetSize() else 1.0
        amount = species.getInitialConcentration() * size
    else:
        amount = 0.0
    return amount * AVOGADRO

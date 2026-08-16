"""Clean-room SBML reader for MOOSE.

Design
------
Standard SBML is the baseline; MOOSE-specific annotations are optional
enrichment, never required. The reader is a general SBML->ODE compiler:

1. **Normalize** (``normalize.py``) via libsbml converters -- inline function
   definitions, promote local parameters, fold initial assignments, sort
   rules -- so the mapper sees one canonical shape.
2. **Map** (``reader.py``): compartments -> CubeMesh, species -> Pool/BufPool,
   reactions and rules -> MOOSE objects. Each kinetic law is decomposed
   **symbolically** (``symbolic.py`` -- build a sympy expression from the AST
   and inspect its algebra): a polynomial with <=2 monomials is mass-action
   (exact Kf/Kb; a modifier appearing multiplicatively is a *catalyst*, mapped
   as a buffered substrate+product), and a degree-1 rational in one substrate
   is Michaelis-Menten (kcat/Km) -> native Reac / MMenz. The extraction is
   numerically self-verified (``identify.py``); anything not recognized (or
   that fails verification) is compiled to a MOOSE ``Function`` (exprtk via
   ``mathconv.py``) driving the same Ksolve ODE system. Rate rules and
   assignment rules likewise become Functions.
3. **Report** (``report.py``): every construct that could not be represented
   faithfully (events, algebraic rules, unresolved symbols) is recorded, so
   the reader never silently mis-simulates.

Multiple compartments
---------------------
A model with several compartments gets one Ksolve+Dsolve+Stoich per
compartment. Reactions that span compartments (a reactant and product in
different compartments -- i.e. transport) are placed in the smaller-volume
compartment and then split by ``fixXreacs`` into per-compartment halves coupled
through proxy pools; the two solvers exchange those proxies each tick, and a
Dsolve mesh junction couples any species carrying a diffusion constant. Because
that exchange is operator-split at the solver tick, multi-compartment models
need a fine simulation ``dt`` for accuracy (and use the default explicit
integrator -- LSODA cannot cross the junction).

A cross-compartment coupling that is *not* native mass-action is reported
unsupported and no solver is built: the Function fallback reads pools by index
within one solver's scope and cannot reach across compartments.

Hard limits (reported, not faked): SBML discrete events and algebraic rules
have no MOOSE equivalent; non-mass-action cross-compartment coupling; and
genuine diffusion needs a diffusion constant, which core SBML does not encode.
"""

from .reader import SBMLHandler, SBMLValidationError
from .report import LoadReport

__all__ = ['SBMLHandler', 'SBMLValidationError', 'LoadReport']

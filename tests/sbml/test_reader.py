"""Tests for the clean-room SBML reader (moose.io.sbml).

These exercise the general SBML->ODE path end-to-end: load a model, build the
MOOSE object graph, simulate, and check the trajectory. Where libRoadRunner is
installed we compare against it as an independent reference; otherwise we fall
back to analytic / conservation checks that need no external simulator.
"""
import os

import numpy as np
import pytest

import moose
from moose.io.sbml import SBMLHandler

HERE = os.path.dirname(os.path.abspath(__file__))


def _sim(root, species_path, runtime, npts=200):
    el = moose.element(species_path)
    tab = moose.Table2(root + '/_tab_' + el.name)
    moose.connect(tab, 'requestOut', el, 'getConc')
    moose.setClock(tab.tick, runtime / npts)
    moose.reinit()
    moose.start(runtime)
    v = np.array(tab.vector)
    return np.linspace(0, runtime, len(v)), v


def test_massaction_decay_and_conservation():
    """00001-sbml-l3v1: S1 -> S2 with rate k1*S1 (k1=1). S1 must decay as
    exp(-t) and S1+S2 must be conserved."""
    h = SBMLHandler()
    h.read(os.path.join(HERE, '00001-sbml-l3v1.xml'), '/decay')
    assert h.report.fully_supported

    t, s1 = _sim('/decay', '/decay/compartment/S1', 5.0)
    _, s2 = _sim('/decay', '/decay/compartment/S2', 5.0)
    # S1 decays as exp(-k1*t) with k1 = 1: log(S1) is linear with slope -1.
    slope = np.polyfit(t, np.log(s1), 1)[0]
    assert abs(slope - (-1.0)) < 1e-2
    # mass conservation: S1 + S2 constant
    total = s1 + s2
    assert np.allclose(total, total[0], rtol=1e-3)


def test_kholodenko_loads_and_runs():
    """A real MAPK-cascade model must load, build pools, and simulate to a
    finite, non-negative trajectory."""
    h = SBMLHandler()
    h.read(os.path.join(HERE, 'Kholodenko.sbml'), '/mapk')
    pools = moose.wildcardFind('/mapk/##[ISA=PoolBase]')
    assert len(pools) > 0
    _, v = _sim('/mapk', pools[0].path, 100.0)
    assert np.all(np.isfinite(v))
    assert np.all(v >= -1e-9)


def test_report_flags_events_unsupported(tmp_path):
    """Events have no MOOSE equivalent and must be reported, not silently
    dropped."""
    sbml = '''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">
  <model id="ev">
    <listOfCompartments><compartment id="c" size="1"/></listOfCompartments>
    <listOfSpecies>
      <species id="S" compartment="c" initialConcentration="1"/>
    </listOfSpecies>
    <listOfEvents>
      <event id="e"><trigger><math xmlns="http://www.w3.org/1998/Math/MathML">
        <apply><gt/><csymbol encoding="text"
          definitionURL="http://www.sbml.org/sbml/symbols/time">t</csymbol>
          <cn>5</cn></apply></math></trigger>
        <listOfEventAssignments>
          <eventAssignment variable="S"><math
            xmlns="http://www.w3.org/1998/Math/MathML"><cn>0</cn></math>
          </eventAssignment>
        </listOfEventAssignments>
      </event>
    </listOfEvents>
  </model>
</sbml>'''
    p = tmp_path / 'ev.xml'
    p.write_text(sbml)
    h = SBMLHandler()
    h.read(str(p), '/ev')
    assert any('event' in u.lower() for u in h.report.unsupported)


def test_massaction_maps_to_native():
    """The mass-action decay model must be recognized symbolically and built as
    a native Reac, not the Function fallback."""
    h = SBMLHandler()
    h.read(os.path.join(HERE, '00001-sbml-l3v1.xml'), '/nat')
    assert h.report.reactions_native >= 1
    assert h.report.reactions_function == 0
    assert moose.wildcardFind('/nat/##[ISA=Reac]')


def test_symbolic_catalyst_extraction():
    """A reaction S -> P with rate k*E*S (E a modifier) must be recognized as
    mass-action with E identified as a catalyst and Kf extracted exactly."""
    import libsbml
    from moose.io.sbml import symbolic

    doc = libsbml.SBMLDocument(3, 1)
    m = doc.createModel()
    c = m.createCompartment()
    c.setId('c'); c.setConstant(True); c.setSize(1.0); c.setSpatialDimensions(3)
    for sid in ('S', 'P', 'E'):
        sp = m.createSpecies()
        sp.setId(sid); sp.setCompartment('c'); sp.setInitialConcentration(1.0)
        sp.setConstant(False); sp.setBoundaryCondition(False)
        sp.setHasOnlySubstanceUnits(False)
    r = m.createReaction()
    r.setId('r'); r.setReversible(False)
    ref = r.createReactant(); ref.setSpecies('S'); ref.setStoichiometry(1); ref.setConstant(True)
    prd = r.createProduct(); prd.setSpecies('P'); prd.setStoichiometry(1); prd.setConstant(True)
    r.createModifier().setSpecies('E')
    kl = r.createKineticLaw()
    kl.setMath(libsbml.parseL3Formula('k * E * S'))
    p = kl.createLocalParameter(); p.setId('k'); p.setValue(0.5)

    res = symbolic.analyze(r, {'c': 1.0})
    assert res is not None and res['kind'] == 'massaction'
    assert res['catalysts'] == {'E'}
    assert abs(res['Kf_val'] - 0.5) < 1e-9
    assert abs(res['Kb_val']) < 1e-12


@pytest.mark.parametrize('fname', ['00001-sbml-l3v1.xml', 'Kholodenko.sbml'])
def test_matches_roadrunner(fname):
    """When libRoadRunner is available, MOOSE trajectories must agree with it
    for these mass-action / Michaelis-Menten models."""
    rr = pytest.importorskip('roadrunner')
    rr.Logger.setLevel(rr.Logger.LOG_FATAL)
    path = os.path.join(HERE, fname)

    import libsbml
    model = libsbml.readSBML(path).getModel()
    sids = [model.getSpecies(i).getId() for i in range(model.getNumSpecies())]
    runtime = 100.0

    r = rr.RoadRunner(path)
    r.timeCourseSelections = ['time'] + ['[%s]' % s for s in sids]
    ref = r.simulate(0, runtime, 201)

    # MOOSE reports concentration in SI mM; RoadRunner in the model's own units.
    # Convert: mM = value * substance_scale / size_scale.
    from moose.io.sbml import units as sbml_units
    ss = sbml_units.substance_scale(model)

    root = '/rr_' + fname.split('.')[0].replace('-', '_')
    SBMLHandler().read(path, root)
    grid = np.linspace(0, runtime, 40)
    for i, sid in enumerate(sids):
        el = [e for e in moose.wildcardFind('%s/##/%s[0]' % (root, sid))
              if 'Pool' in e.className]
        if not el:
            continue
        comp = model.getCompartment(model.getSpecies(sid).getCompartment())
        cfac = ss / sbml_units.size_scale(comp)
        t, v = _sim(root, el[0].path, runtime)
        rvals = np.interp(grid, ref[:, 0], ref[:, i + 1]) * cfac
        mvals = np.interp(grid, t, v)
        denom = max(np.max(np.abs(rvals)), 1e-12)
        assert np.max(np.abs(rvals - mvals)) / denom < 2e-2, sid

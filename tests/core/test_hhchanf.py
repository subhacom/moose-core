# Filename: test_hhchanf.py
# Description:
# Author: Subhasis Ray
# Created: Wed Jan 29 13:15:18 2025 (+0530)
#

"""Tests HHChannelF class
Usage: pytest test_hhchanf.py
"""
import moose
import math
import pytest
from ephys import create_voltage_clamp, setup_step_command


@pytest.fixture
def container():
    """
    Setup:
        Create a Neutral object as model container and ce to it before
        yielding.

    Teardown:
        ce back to parent and delete the model
        container.

    """
    ret = moose.Neutral('/test')
    moose.ce(ret)
    yield ret
    moose.ce('..')
    moose.delete(ret)


@pytest.fixture
def channel(container):
    """
    Setup:
        Create a HHChannelF object.

    Teardown:
        Nothing. The container gets deleted
    """
    return moose.HHChannelF('ch')


def test_hhgatef_creation(container):
    channel = moose.HHChannelF('ch0')
    channel.Xpower = 1
    channel.Ypower = 2
    assert math.isclose(channel.Xpower, 1), 'Xpower not set'
    assert moose.exists('ch0/gateX'), 'gateX object does not exist'
    assert math.isclose(channel.Ypower, 2), 'Ypower not set'
    assert moose.exists('ch0/gateY'), 'gateY object does not exist'

    
def test_hhgatef_alpha_beta(channel):
    """Test set/get alpha and beta expressions"""
    channel.Xpower = 1
    xgate = moose.element(f'{channel.path}/gateX')
    alpha = '0.01 * (10 - v) / exp((10 - v)/10 - 1)'
    beta = '0.125 * exp( - v/80)'
    xgate.alpha = alpha
    xgate.beta = beta
    assert xgate.alpha == alpha, 'alpha not set'
    assert xgate.beta == beta, 'beta not set'
    assert xgate.tau == '', 'tau not reset'
    assert xgate.inf == '', 'inf not reset'
    

def test_hhgatef_tau_inf(channel):
    """Test set/get tau and inf expressions"""
    channel.Xpower = 1
    xgate = moose.element(f'{channel.path}/gateX')
    alpha = '0.01 * (10 - v) / exp((10 - v)/10 - 1)'
    beta = '0.125 * exp( - v/80)'
    xgate.tau = alpha
    xgate.inf = beta
    assert xgate.tau == alpha, 'tau not set'
    assert xgate.inf == beta, 'inf not set'
    assert xgate.alpha == '', 'alpha not reset'
    assert xgate.beta == '', 'beta not reset'

def test_complex_expr(channel):    
    channel.Xpower = 3
    channel.Ypower = 1
    channel.Ek = 55e-3
    mgate = moose.element(f'{channel.path}/gateX')
    hgate = moose.element(f'{channel.path}/gateY')
    # These are based on NeuroML - for some reason the A parameter is
    # double that in the paper
    mtau ='~(alpha:=1500 * exp(81 *((v - 10e-3) - (-39e-3))), beta:=1500 * exp(-66 * (v - 10e-3) - (-39e-3)), alpha + beta == 0? 0: min(1/(alpha+beta), 5e-5))'
    mgate.tau = mtau
    minf = '~(alpha:=1500 * exp(81 *((v - 10e-3) - (-39e-3))), beta:=1500 * exp(-66 * (v - 10e-3) - (-39e-3)), alpha/(alpha+beta))'
    mgate.inf = minf
    htau = '~(alpha:=120*exp(((v - 10e-3) - (-0.04))/(-0.01123596)), beta:=120 * exp(((v - 10e-3) - (-0.04)) / 0.01123596), alpha + beta == 0? 0: min(1/(alpha+beta), 2.25e-4))'
    hgate.tau = htau
    hinf = '~(alpha:=120*exp(((v - 10e-3) - (-0.04))/(-0.01123596)), beta:=120 * exp(((v - 10e-3) - (-0.04)) / 0.01123596), alpha/(alpha+beta))'
    hgate.inf = hinf
    assert mgate.tau == mtau, 'tau not set'
    assert mgate.inf == minf, 'inf not set'
    assert hgate.tau == htau, 'tau not set'
    assert hgate.inf == hinf, 'inf not set'
    
    
def test_hh_k_vclamp(container, steptime=5.0):
    """Simulate a voltage clamp with Hodhkin and Huxley's K+ channel.

    Parameters
    ----------
    steptime: float
        Time of the voltage step
    """
    comp = moose.Compartment('comp0')
    chan = moose.HHChannelF(f'{comp.path}/K')
    moose.connect(chan, 'channel', comp, 'channel')
    chan.Gbar = 36.0
    chan.Ek = -12.0
    chan.Xpower = 4
    n_gate = moose.element(f'{chan.path}/gateX')
    n_gate.alpha = '0.01 * (10 - v) / (exp((10-v)/10) - 1)'
    n_gate.beta = '0.125 * exp(-v/80)'
    comp.Em = 0  # Hodgkin and Huxley used resting voltage as 0
    comp.Vm = 0
    comp.initVm = 0
    comp.Cm = 1
    comp.Rm = 1 / 0.3  # G_leak is 0.3 mS/cm^2
    dt = 0.01
    for tick in range(8):
        moose.setClock(tick, dt)
    # vclamp needs to be setup after clock dt because it uses dt
    vclamp, command, commandtab = create_voltage_clamp(comp)
    # Hodgkin and Huxley 1952, "Currents carried by sodium and
    # potassium ions ...", find maximum conductance try 0-100 mV at 10 mV steps
    simtime = 100.0 + steptime
    # Precalculated steady state K conductance
    vm_gk = {
        0: 0.3666444556069115,
        # 10: 1.8401,
        20: 5.287062775435651,
        30: 10.176967265218654,
        40: 15.220230442170081,
        50: 19.596740114488156,
        60: 23.1009379487905,
        70: 25.821065780738117,
        80: 27.91721528657564,
        90: 29.537324370047777,
        100: 30.79812392487768,
    }

    for vstep, gk in vm_gk.items():
        setup_step_command(command, 0.0, delay=steptime, level=vstep)
        moose.reinit()
        moose.start(simtime)
        assert math.isclose(
            gk, chan.Gk, abs_tol=1e-6
        ), f'Vm={vstep} Gk={chan.Gk}, expected={gk}'


def test_hh_na_vclamp(container, steptime=5.0):
    """Test the evaluation of hhchannel conductance using Hodgkin and
    Huxleys Na channel model"""
    comp = moose.Compartment('comp0')
    chan = moose.HHChannelF(f'{comp.path}/Na')
    moose.connect(chan, 'channel', comp, 'channel')
    chan.Gbar = 120.0
    chan.Ek = 115.0
    chan.Xpower = 3
    chan.Ypower = 1
    m_gate = moose.element(f'{chan.path}/gateX')
    h_gate = moose.element(f'{chan.path}/gateY')
    m_gate.alpha = '0.1 * (25 - v) / (exp((25 - v) / 10) - 1)'
    m_gate.beta = '4 * exp(- v / 18)'
    h_gate.alpha = '0.07 * exp(- v / 20)'
    h_gate.beta = '1 / (exp((30 - v) / 10) + 1)'
    comp.Em = 0  # Hodgkin and Huxley used resting voltage as 0
    comp.Vm = 0
    comp.initVm = 0
    comp.Cm = 1
    comp.Rm = 1 / 0.3  # G_leak is 0.3 mS/cm^2
    dt = 0.01
    for tick in range(8):
        moose.setClock(tick, dt)
    # This must be after setting clock dt because parameters depend on dt
    vclamp, command, commandtab = create_voltage_clamp(comp)
    # Hodgkin and Huxley 1952, "Currents carried by sodium and
    # potassium ions ...", find maximum conductance
    simtime = 100.0 + steptime
    # Precalculated steady state Na conductance
    vm_gna = {
        0: 0.010609192838829854,
        10: 0.12443211458767735,
        20: 0.527787746095522,
        30: 0.8966173682113173,
        40: 0.8361209504788413,
        50: 0.5983994566094422,
        60: 0.3893933499995894,
        70: 0.24432304663465443,
        80: 0.15080726000137287,
        90: 0.09232255449417581,
        100: 0.0562750224683144,
    }

    for vstep, gna in vm_gna.items():
        setup_step_command(command, 0.0, delay=steptime, level=vstep)
        moose.reinit()
        moose.start(simtime)
        assert math.isclose(
            gna, chan.Gk, abs_tol=1e-6
        ), f'Vm={vstep} Gk={chan.Gk}, expected={gna}'
    

def test_hhchanf_eval(container):
    """Test the evaluation of hhchannel conductance using Hodgkin and
    Huxleys Na channel model"""
    comp = moose.Compartment('comp0')
    chan = moose.HHChannelF(f'{comp.path}/Na')
    moose.connect(chan, 'channel', comp, 'channel')
    chan.Gbar = 120.0
    chan.Ek = 115.0
    chan.Xpower = 3
    chan.Ypower = 1
    m_gate = moose.element(f'{chan.path}/gateX')
    h_gate = moose.element(f'{chan.path}/gateY')
    m_gate.alpha = '0.1 * (25 - v) / (exp((25 - v) / 10) - 1)'
    m_gate.beta = '4 * exp(- v / 18)'
    h_gate.alpha = '0.07 * exp(- v / 20)'
    h_gate.beta = '1 / (exp((30 - v) / 10) + 1)'
    comp.Em = 0  # Hodgkin and Huxley used resting voltage as 0
    comp.Vm = 0
    comp.initVm = 0
    comp.Cm = 1
    comp.Rm = 1 / 0.3  # G_leak is 0.3 mS/cm^2
    moose.reinit()
    t = 0
    tdelta = comp.dt * 1000
    for ii in range(10):
        moose.start(tdelta)
        t += tdelta
        print(f'time={t} Gk={chan.Gk} Vm={comp.Vm}')
    comp.inject = 1.0
    print('Starting current injection...')
    for ii in range(10):
        moose.start(tdelta)
        t += tdelta
        print(f'time={t} Gk={chan.Gk} Vm={comp.Vm}')


if __name__ == '__main__':
    test_hhgatef_creation()
    test_hhgatef_alpha_beta()
    test_hhgatef_tau_inf()
    test_hh_k_vclamp()
    test_hh_na_vclamp()
    test_hhchanf_eval()
#
# test_hhchanf.py ends here

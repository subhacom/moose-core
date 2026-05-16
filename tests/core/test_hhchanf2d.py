# Filename: test_hhchanf2d.py
# Description:
# Author: Subhasis Ray
# Created: Wed Jan 29 13:15:18 2025 (+0530)
#
# Last updated: Sat May 16 09:53:25 IST 2026
# Updated by: Subhasis Ray and Claude 4.6
#



"""Tests HHChannelF2D class.

Usage: pytest test_hhchanf2d.py
"""
import moose
import math
from ephys import create_voltage_clamp, setup_step_command



def test_hhgatef2d_creation():
    cwe = moose.getCwe()
    container = moose.Neutral('/test')
    moose.ce(container)
    chan = moose.HHChannelF2D('ch')

    chan.Xpower = 1
    chan.Ypower = 2
    assert math.isclose(chan.Xpower, 1), 'Xpower not set'
    assert moose.exists(f'{chan.path}/gateX'), 'gateX object does not exist'
    assert math.isclose(chan.Ypower, 2), 'Ypower not set'
    assert moose.exists(f'{chan.path}/gateY'), 'gateY object does not exist'

    moose.ce(cwe)
    moose.delete(container)


def test_alpha_beta():
    """Test set/get alpha and beta expressions"""
    cwe = moose.getCwe()
    container = moose.Neutral('/test')
    moose.ce(container)
    chan = moose.HHChannelF2D('ch')
    chan.Xpower = 1
    xgate = moose.element(f'{chan.path}/gateX')
    # These are from KCa chan in Eric DeSchutter's granule cell model
    alpha = '2500 / (1 + 1.5e-3 * exp(-85*v)/c)'
    beta = '1500 / (1 + c / (1.5e-4 * exp (-77*v)))'
    xgate.alphaExpr = alpha
    xgate.betaExpr = beta
    assert xgate.alphaExpr == alpha, 'alpha not set'
    assert xgate.betaExpr == beta, 'beta not set'
    assert xgate.tauExpr == '', 'tau not reset'
    assert xgate.infExpr == '', 'inf not reset'
    moose.ce(cwe)
    moose.delete(container)


def test_complex_expr():
    cwe = moose.getCwe()
    container = moose.Neutral('/test')
    moose.ce(container)
    chan = moose.HHChannelF2D('ch')

    chan.Xpower = 3
    chan.Ypower = 1
    chan.Ek = 55e-3
    mgate = moose.element(f'{chan.path}/gateX')
    hgate = moose.element(f'{chan.path}/gateY')
    # These are based on NeuroML - for some reason the A parameter is
    # double that in the paper
    mtau ='~(alpha:=1500 * exp(81 *((v - 10e-3) - (-39e-3))), beta:=1500 * exp(-66 * (v - 10e-3) - (-39e-3)), alpha + beta == 0? 0: min(1/(alpha+beta), 5e-5))'
    mgate.tau = tau
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

    moose.ce(cwe)
    moose.delete(container)


def test_tau_inf():
    """Test set/get tau and inf expressions"""
    cwe = moose.getCwe()
    container = moose.Neutral('/test')
    moose.ce(container)
    chan = moose.HHChannelF2D('ch')

    chan.Xpower = 1
    xgate = moose.element(f'{chan.path}/gateX')
    alpha = '2500 / (1 + 1.5e-3 * exp(-85*v)/c)'
    beta = '1500 / (1 + c / (1.5e-4 * exp (-77*v)))'
    tau = f'1/({alpha} + {beta})'
    minf = f'({alpha}) / ({alpha} + {beta})'
    xgate.tauExpr = tau
    xgate.infExpr = minf
    assert xgate.tauExpr == tau, 'tau not set'
    assert xgate.infExpr == minf, 'inf not set'
    assert xgate.alphaExpr == '', 'alpha not reset'
    assert xgate.betaExpr == '', 'beta not reset'

    moose.ce(cwe)
    moose.delete(container)


def test_vclamp(steptime=5.0):
    """Simulate a voltage clamp experiment with fixed Ca conc

    Parameters
    ----------
    steptime: float
        Time of the voltage step
    """
    cwe = moose.getCwe()
    container = moose.Neutral('/test')
    moose.ce(container)
    comp = moose.Compartment('comp0')
    sarea = 4 * math.pi * (10e-3) ** 3 / 3
    chan = moose.HHChannelF2D(f'{comp.path}/K')
    moose.connect(chan, 'channel', comp, 'channel')
    chan.Gbar = sarea * 17.9811e-3 * 1e4  # 17.9811 mS/cm^2
    print(f'Gbar={chan.Gbar}')
    chan.Ek = -90e-3
    chan.Xpower = 1
    chan.Xindex = 'VOLT_C1_INDEX'
    n_gate = moose.element(f'{chan.path}/gateX')
    n_gate.alphaExpr = '2500 / (1 + 1.5e-3 * exp(-85*(v-0.01))/c)'
    n_gate.betaExpr = '1500 / (1 + c / (1.5e-4 * exp (-77*(v-0.01))))'
    comp.Em = -65e-3  # Hodgkin and Huxley used resting voltage as 0
    comp.initVm = -65e-3
    comp.Cm = sarea * 1e-6 * 1e4  # 1 uS/cm^2
    comp.Rm = 1 / (
        sarea * 0.0330033e-3 * 1e4
    )  # passive conductance is 0.0330033 mS/cm^2
    capool = moose.CaConc(f'{comp.path}/Ca')
    capool.CaBasal = 1e-9
    moose.connect(capool, 'concOut', chan, 'concen')
    vclamp, command, commandtab = create_voltage_clamp(comp)
    simtime = 100.0 + steptime
    # Precalculated steady state K conductance
    v_commands = [
        -75e-3,
        -65e-3,
        -55e-3,
        -45e-3,
        -35e-3,
        -25e-3,
        -15e-3,
        0,
        15e-3,
    ]
    ca_commands = [0.1e-9, 0.3e-9, 0.6e-9, 0.9e-9, 1.2e-9]

    def alpha(v, c):
        return 2500 / (1 + 1.5e-3 * math.exp(-85 * (v - 0.01)) / c)

    def beta(v, c):
        return 1500 / (1 + c / (1.5e-4 * math.exp(-77 * (v - 0.01))))

    vm_ca_gk = [
        [
            chan.Gbar * alpha(v, c) / (alpha(v, c) + beta(v, c))
            for c in ca_commands
        ]
        for v in v_commands
    ]

    for ii, vstep in enumerate(v_commands):
        setup_step_command(command, comp.Em, delay=steptime, level=vstep)
        for jj, ca in enumerate(ca_commands):
            capool.CaBasal = ca
            moose.reinit()
            moose.start(simtime)
            gk = vm_ca_gk[ii][jj]
            print('V', vstep, 'Ca', ca, 'Gk: Computed', gk, 'Simulated', chan.Gk)
            assert math.isclose(
                gk, chan.Gk, abs_tol=1e-10
            ), f'Vm={vstep} [Ca]={ca} Gk={chan.Gk}, expected={gk}'
    moose.ce(cwe)
    moose.delete(container)


# ── Moczydlowski-Latorre BK channel ──────────────────────────────────────────
#
# Reference: Moczydlowski E & Latorre R (1983) J Gen Physiol 82:511-542.
# Original GENESIS script: neurokit/prototypes/MoczydKC.g (De Schutter).
#
# Single activation gate X, Xpower=1, VOLT_C1_INDEX dependency.
# Units: voltage in V, Ca2+ in mM, rates in s^-1, temperature 37 °C.
#
# Alpha and beta rates:
#   alpha(V, Ca) = 480*Ca / (Ca + 0.180*exp(-0.84*ZFbyRT*V))
#   beta(V, Ca)  = 280 / (1 + Ca / (0.011*exp(-1.00*ZFbyRT*V)))
# where ZFbyRT = 23210 / (273.15 + 37) ≈ 74.83 V^-1  (= 2*F/R/T).
#
# Non-separability: V and Ca appear additively in the denominators, so
# f(V, Ca) cannot be factored into g(V)*h(Ca).

_BK_EK     = -0.085
_BK_ZFBYRT = 23210.0 / (273.15 + 37.0)   # ≈ 74.83 V^-1

_BK_ALPHA = (
    f'480 * c / (c + 0.180 * exp(-0.84 * {_BK_ZFBYRT:.8f} * v))'
)
_BK_BETA = (
    f'280 / (1 + c / (0.011 * exp(-1.00 * {_BK_ZFBYRT:.8f} * v)))'
)


def make_BK_KC_prototype(path='/library/BK_KC_Moczyd'):
    """
    Build (or retrieve) an HHChannelF2D prototype for the Moczydlowski-Latorre
    BK (Big-conductance K+) channel.

    The call is idempotent: repeated calls with the same path return the
    existing element without rebuilding it.

    Parameters
    ----------
    path : str
        MOOSE element path; must be under ``/library``.

    Returns
    -------
    moose.element (HHChannelF2D)
        Prototype with Ek=-85 mV, Gbar=1.0 (scale factor), single X gate.
    """
    if not moose.exists('/library'):
        moose.Neutral('/library')
    if moose.exists(path):
        return moose.element(path)

    chan = moose.HHChannelF2D(path)
    chan.Ek   = _BK_EK
    chan.Gbar = 1.0           # scale factor; multiply by actual conductance at insert time

    chan.Xindex = 'VOLT_C1_INDEX'   # gate depends on Vm (dim 0) and Ca (dim 1)
    chan.Xpower = 1.0               # triggers gate creation

    gate = moose.element(f'{path}/gateX[0]')
    gate.alphaExpr = _BK_ALPHA
    gate.betaExpr  = _BK_BETA

    return moose.element(path)


# ── prototype validation tests ────────────────────────────────────────────────

def _reset_library():
    if moose.exists('/library'):
        moose.delete('/library')


def test_bk_prototype_is_hhchannelf2d():
    _reset_library()
    assert make_BK_KC_prototype().className == 'HHChannelF2D'


def test_bk_prototype_in_library():
    _reset_library()
    assert make_BK_KC_prototype().path.startswith('/library')


def test_bk_prototype_idempotent():
    _reset_library()
    p1 = make_BK_KC_prototype()
    p2 = make_BK_KC_prototype()
    assert p1.path == p2.path


def test_bk_xindex():
    _reset_library()
    assert make_BK_KC_prototype().Xindex == 'VOLT_C1_INDEX'


def test_bk_xpower():
    _reset_library()
    assert math.isclose(make_BK_KC_prototype().Xpower, 1.0)


def test_bk_gate_expressions():
    _reset_library()
    proto = make_BK_KC_prototype()
    gate = moose.element(proto.path + '/gateX[0]')
    assert gate.alphaExpr == _BK_ALPHA
    assert gate.betaExpr  == _BK_BETA


def test_bk_ek():
    _reset_library()
    assert math.isclose(make_BK_KC_prototype().Ek, _BK_EK)


# ── voltage clamp simulation ──────────────────────────────────────────────────

def test_bk_voltage_clamp(steptime=0.05):
    """
    Voltage clamp with fixed Ca2+: simulated steady-state Gk must match the
    analytical value Gbar * alpha/(alpha+beta) at each (V, Ca) pair.

    Ca2+ is held constant by setting CaConc.CaBasal (mM) before each reinit.
    The BK channel equilibrates in ~3 ms, so a 50 ms step is ample.
    """
    cwe = moose.getCwe()
    if moose.exists('/test_bk_vclamp'):
        moose.delete('/test_bk_vclamp')
    container = moose.Neutral('/test_bk_vclamp')
    moose.ce(container)

    comp = moose.Compartment('comp')
    comp.Cm     = 1e-9    # F
    comp.Rm     = 1e8     # Ω
    comp.Em     = -65e-3
    comp.initVm = -65e-3

    chan = moose.HHChannelF2D(f'{comp.path}/BK')
    chan.Ek     = _BK_EK
    chan.Gbar   = 1e-6    # S
    chan.Xindex = 'VOLT_C1_INDEX'
    chan.Xpower = 1.0
    gate = moose.element(f'{chan.path}/gateX[0]')
    gate.alphaExpr = _BK_ALPHA
    gate.betaExpr  = _BK_BETA

    moose.connect(chan, 'channel', comp, 'channel')

    capool = moose.CaConc(f'{comp.path}/Ca')
    moose.connect(capool, 'concOut', chan, 'concen')

    vclamp, command, commandtab = create_voltage_clamp(comp)

    # Voltage steps (V) and Ca2+ values (mM = numerically what the expression sees)
    v_steps  = [-75e-3, -35e-3, -15e-3, 0.0, 35e-3]
    ca_steps = [0.01, 0.05, 0.1, 0.3]    # 10, 50, 100, 300 µM

    def alpha(v, c):
        return 480 * c / (c + 0.180 * math.exp(-0.84 * _BK_ZFBYRT * v))

    def beta(v, c):
        return 280 / (1 + c / (0.011 * math.exp(-1.00 * _BK_ZFBYRT * v)))

    simtime = steptime + 0.1   # 100 ms equilibration after the step

    for vstep in v_steps:
        setup_step_command(command, comp.Em, delay=steptime, level=vstep)
        for ca in ca_steps:
            capool.CaBasal = ca
            moose.reinit()
            moose.start(simtime)
            a = alpha(vstep, ca)
            b = beta(vstep, ca)
            expected_gk = chan.Gbar * a / (a + b)
            assert math.isclose(chan.Gk, expected_gk, rel_tol=1e-4), (
                f'V={vstep*1e3:.0f} mV  Ca={ca} mM: '
                f'Gk={chan.Gk:.6g} S  expected={expected_gk:.6g} S')

    moose.ce(cwe)
    moose.delete(container)


if __name__ == '__main__':
    test_hhgatef2d_creation()
    test_alpha_beta()
    test_tau_inf()
    test_vclamp()
    test_bk_prototype_is_hhchannelf2d()
    test_bk_gate_expressions()
    test_bk_voltage_clamp()
    print('All HHChannelF2D tests passed.')
#
# test_hhchanf2d.py ends here

#
# test_hhchanf2d.py ends here

"""
compare_moose_lems.py
=====================
Compare MOOSE HHGate kinetics (inf, tau) against jNeuroML/LEMS reference
for OSB BlueBrain channel files.

The pyneuroml channel analyser runs a voltage-ramp simulation and records
inf(t) and tau(t) from the LEMS gate model.  Voltage at each time step is
obtained from the recorded `rampCellPop0[0]/v` quantity.  We then
interpolate MOOSE's gate tables at the same voltage points and compare.

For SK_E2 (Ca-only gate) the LEMS run uses a fixed Ca concentration.  The
inf curve is constant vs voltage but varies with Ca, so we compare the LEMS
value at the reference Ca_conc against MOOSE's 1D table at the same point.
"""
import os
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import moose
import pyneuroml.analysis.NML2ChannelAnalysis as nca
from moose.neuroml2.reader import NML2Reader
from neuroml.loaders import read_neuroml2_file

# Channel files live in the same directory as this script.
_CHAN_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
CHAN_FILES = {
    "NaTa_t": {"gates": ["m", "h"], "ca_only": False},
    "Ca_HVA":  {"gates": ["m", "h"], "ca_only": False},
    "Ih":      {"gates": ["m"],      "ca_only": False},
    "K_Tst":   {"gates": ["m", "h"], "ca_only": False},
    "Nap_Et2": {"gates": ["m", "h"], "ca_only": False},
    "SK_E2":   {"gates": ["z"],      "ca_only": True},
}
# Im is omitted because jLEMS crashes with fcond=-Infinity at high voltages.

TEMPERATURE = 32        # degC — matches the GranuleCell network temperature
# SK_E2 half-activation: 4.3e-10 mol/cm³ = 4.3e-4 mol/m³ = 0.43 µM.
# Test near this point for good sensitivity; LEMS ca_conc parameter is in mM,
# and 1 mM = 1 mol/m³ numerically.
CA_CONC_MM  = 4.3e-4   # mM ≡ 0.43 µM — SK_E2 half-activation, inf ≈ 0.5

# MOOSE gate table parameters
# Voltage tables are in V (passed to pint as "V").
# Concentration tables are in mol/m³ (passed to pint as "mole/meter**3").
# Numerically: 1 mM = 1 mol/m³, so mM ≡ mol/m³ as numbers.
# CMAX = 1e-3 mol/m³ = 1 µM covers SK_E2's physiological [Ca²⁺] range.
# With CDIVS=5000 the step is ~2e-7 mol/m³ (0.2 nM) — fine enough to
# resolve the SK_E2 transition at 0.43 µM.
VMIN  = -150e-3    # V
VMAX  =  100e-3    # V
VDIVS = 3000
CMIN  =    0.0     # mol/m³
CMAX  =    1e-3    # mol/m³ (= 1 µM)
CDIVS =  5000

MAX_REL_ERROR = 0.02   # 2 % — threshold for PASS/FAIL


def _lems_args(min_v=-100, max_v=100):
    return argparse.Namespace(
        v=False,
        min_v=min_v, max_v=max_v, step_target_voltage=20,
        temperature=TEMPERATURE, duration=100,
        clamp_delay=10, clamp_duration=80, clamp_base_voltage=-70,
        erev=50, scale_dt=1, ca_conc=CA_CONC_MM, dat_suffix='',
        iv_curve=False, norun=False, nogui=True, html=False, md=False,
    )


def load_lems_kinetics(chan_name, gates):
    """Run jNeuroML channel analyser and extract (v_mV, inf, tau_ms) per gate."""
    a = _lems_args()
    chan_file = os.path.join(_CHAN_DIR, f"{chan_name}.channel.nml")
    channels = nca.get_channels_from_channel_file(chan_file)
    lems_file = nca.make_lems_file(channels[0], a)
    results = nca.run_lems_file(lems_file, verbose=False)
    if not results:
        raise RuntimeError(f"LEMS returned no data for {chan_name}")

    t = np.array(results["t"])   # seconds

    # Extract voltage from the ramp cell recording
    v_key = "rampCellPop0[0]/v"
    if v_key in results:
        v_mv = np.array(results[v_key]) * 1e3   # V → mV
    else:
        # Fallback: reconstruct from known ramp parameters.
        # rampCell is set up by generate_lems_channel_analyser to go from
        # min_v to max_v over the full simulation length (no delay).
        v_mv = a.min_v + (a.max_v - a.min_v) * t / t[-1]

    # Restrict to the voltage range we care about
    mask = (v_mv >= a.min_v) & (v_mv <= a.max_v)
    v_mv = v_mv[mask]
    t = t[mask]

    kinetics = {}
    for gate in gates:
        inf_key = next((k for k in results if f"/{chan_name}/{gate}/inf" in k), None)
        tau_key = next((k for k in results if f"/{chan_name}/{gate}/tau" in k), None)
        if inf_key and tau_key:
            inf = np.array(results[inf_key])[mask]
            tau_s = np.array(results[tau_key])[mask]
            # LEMS records tau in seconds; convert to ms
            tau_ms = tau_s * 1e3
            kinetics[gate] = (v_mv, inf, tau_ms)
    return kinetics


def load_moose_kinetics(chan_name, ca_only=False):
    """Load channel into MOOSE /library and return gate (x_axis, inf, tau_ms) dict."""
    if moose.exists("/library"):
        moose.delete("/library")
    chan_file = os.path.join(_CHAN_DIR, f"{chan_name}.channel.nml")
    doc = read_neuroml2_file(chan_file, include_includes=True, verbose=False)
    r = NML2Reader()
    r._vmin = VMIN; r._vmax = VMAX; r._vdivs = VDIVS
    r._cmin = CMIN; r._cmax = CMAX; r._cdivs = CDIVS
    r.doc = doc
    r.importIonChannels(doc, vmin=VMIN, vmax=VMAX, vdivs=VDIVS,
                        cmin=CMIN, cmax=CMAX, cdivs=CDIVS)

    mchan = moose.element(moose.wildcardFind(f"/library/{chan_name}")[0])
    kinetics = {}
    for letter in ["X", "Y", "Z"]:
        power = getattr(mchan, f"{letter}power", 0)
        if power <= 0:
            continue
        mgate = moose.element(f"{mchan.path}/gate{letter}")
        tableA = np.array(mgate.tableA)
        tableB = np.array(mgate.tableB)
        # avoid divide-by-zero
        tableB = np.where(tableB == 0, np.finfo(float).tiny, tableB)
        inf = tableA / tableB
        tau_ms = 1.0 / tableB * 1e3   # s → ms

        if ca_only and letter == "Z":
            # Table indexed by [Ca²⁺]: x-axis is concentration in M
            x = np.linspace(CMIN, CMAX, len(tableA))
            kinetics["z"] = (x, inf, tau_ms)
        else:
            # Table indexed by voltage in V → convert to mV
            x_mv = np.linspace(mgate.min * 1e3, mgate.max * 1e3, len(tableA))
            kinetics[letter.lower()] = (x_mv, inf, tau_ms)
    return kinetics


def compare(chan_name, lems_kin, moose_kin, ca_only=False):
    """Return per-gate relative-error stats between MOOSE tables and LEMS data."""
    res = {}
    # Use insertion order (Python 3.7+), not alphabetical sort.
    # MOOSE creates gates in the same order as the NML2 gate list,
    # and load_lems_kinetics preserves the gate order from the input list.
    lems_gates = list(lems_kin.keys())
    moose_gates = list(moose_kin.keys())

    for lems_gate in lems_gates:
        lems_v, lems_inf, lems_tau = lems_kin[lems_gate]

        # Match MOOSE gate by name or by position
        if lems_gate in moose_kin:
            moose_gate = lems_gate
        elif lems_gates.index(lems_gate) < len(moose_gates):
            moose_gate = moose_gates[lems_gates.index(lems_gate)]
        else:
            print(f"  {chan_name}/{lems_gate}: no MOOSE gate found")
            continue

        moose_v, moose_inf, moose_tau = moose_kin[moose_gate]

        if ca_only:
            # SK_E2: LEMS uses fixed Ca at CA_CONC_MM mM.
            # MOOSE table runs over [CMIN, CMAX] in mol/m³.
            # Since 1 mM = 1 mol/m³, the mM value equals the mol/m³ value.
            ca_m = CA_CONC_MM  # mM == mol/m³ numerically
            moose_inf_at = np.interp(ca_m, moose_v, moose_inf)
            moose_tau_at = np.interp(ca_m, moose_v, moose_tau)
            # LEMS inf/tau should be nearly constant vs voltage for Ca-only gate.
            # Use median of LEMS values (stable region mid-ramp)
            lems_inf_median = np.median(lems_inf)
            lems_tau_median = np.median(lems_tau)
            eps = 1e-12
            inf_rel = abs(moose_inf_at - lems_inf_median) / (abs(lems_inf_median) + eps)
            tau_rel = abs(moose_tau_at - lems_tau_median) / (abs(lems_tau_median) + eps)
            res[lems_gate] = {
                "inf_max_rel": float(inf_rel), "inf_mean_rel": float(inf_rel),
                "tau_max_rel": float(tau_rel), "tau_mean_rel": float(tau_rel),
                "lems_v": lems_v, "lems_inf": lems_inf, "lems_tau": lems_tau,
                "moose_inf": np.full_like(lems_v, moose_inf_at),
                "moose_tau": np.full_like(lems_v, moose_tau_at),
            }
        else:
            # Interpolate MOOSE table at the LEMS voltage sample points
            moose_inf_at = np.interp(lems_v, moose_v, moose_inf)
            moose_tau_at = np.interp(lems_v, moose_v, moose_tau)
            eps = 1e-12
            inf_rel = np.abs(moose_inf_at - lems_inf) / (np.abs(lems_inf) + eps)
            tau_rel = np.abs(moose_tau_at - lems_tau) / (np.abs(lems_tau) + eps)
            res[lems_gate] = {
                "inf_max_rel": float(inf_rel.max()),
                "inf_mean_rel": float(inf_rel.mean()),
                "tau_max_rel": float(tau_rel.max()),
                "tau_mean_rel": float(tau_rel.mean()),
                "lems_v": lems_v, "lems_inf": lems_inf, "lems_tau": lems_tau,
                "moose_inf": moose_inf_at, "moose_tau": moose_tau_at,
            }
    return res


def make_plot(chan_name, comparison, ca_only=False, outdir="."):
    n_gates = len(comparison)
    if not n_gates:
        return
    fig, axes = plt.subplots(2, n_gates, figsize=(5 * n_gates, 8), squeeze=False)
    fig.suptitle(f"{chan_name} — MOOSE (red dashed) vs LEMS (blue)", fontsize=13)
    xlabel = "[Ca²⁺] (M)" if ca_only else "V (mV)"

    for col, (gate_id, r) in enumerate(comparison.items()):
        lv = r["lems_v"]
        axes[0, col].plot(lv, r["lems_inf"], "b-", lw=2, label="LEMS")
        axes[0, col].plot(lv, r["moose_inf"], "r--", lw=1.5, label="MOOSE")
        axes[0, col].set_title(f"gate {gate_id} — inf")
        axes[0, col].set_xlabel(xlabel)
        axes[0, col].set_ylabel("inf")
        axes[0, col].legend(fontsize=8)

        axes[1, col].plot(lv, r["lems_tau"], "b-", lw=2, label="LEMS")
        axes[1, col].plot(lv, r["moose_tau"], "r--", lw=1.5, label="MOOSE")
        axes[1, col].set_title(f"gate {gate_id} — tau (ms)")
        axes[1, col].set_xlabel(xlabel)
        axes[1, col].set_ylabel("tau (ms)")
        axes[1, col].legend(fontsize=8)

    plt.tight_layout()
    out = os.path.join(outdir, f"{chan_name}_comparison.png")
    plt.savefig(out, dpi=100)
    plt.close()
    return out


def main():
    os.makedirs("plots", exist_ok=True)
    all_pass = True
    print(f"\n{'Chan':12s} {'Gate':6s} {'inf_max%':10s} {'tau_max%':10s} {'Status':8s}")
    print("-" * 55)

    for chan_name, info in CHAN_FILES.items():
        gates = info["gates"]
        ca_only = info["ca_only"]
        try:
            lems_kin = load_lems_kinetics(chan_name, gates)
        except Exception as exc:
            print(f"{chan_name:12s}: LEMS FAILED — {exc}")
            all_pass = False
            continue

        try:
            moose_kin = load_moose_kinetics(chan_name, ca_only=ca_only)
        except Exception as exc:
            print(f"{chan_name:12s}: MOOSE FAILED — {exc}")
            all_pass = False
            continue

        cmp = compare(chan_name, lems_kin, moose_kin, ca_only=ca_only)
        if not cmp:
            print(f"{chan_name:12s}: no gates matched")
            all_pass = False
            continue

        make_plot(chan_name, cmp, ca_only=ca_only, outdir="plots")

        for gate_id, r in cmp.items():
            inf_pct = r["inf_max_rel"] * 100
            tau_pct = r["tau_max_rel"] * 100
            ok = (r["inf_max_rel"] < MAX_REL_ERROR and
                  r["tau_max_rel"] < MAX_REL_ERROR)
            status = "PASS" if ok else "FAIL"
            if not ok:
                all_pass = False
            print(f"{chan_name:12s} {gate_id:6s} {inf_pct:8.2f}%  {tau_pct:8.2f}%  {status}")

    print()
    print("Overall:", "PASS" if all_pass else "FAIL")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())

"""
test_osb_channels.py
====================
Pytest wrapper for the MOOSE vs jNeuroML/LEMS channel kinetics comparison.

Requires:
  - Java (jNeuroML is invoked via java -jar ...)
  - pyneuroml (pip install pyneuroml)

Skip the whole module when Java is not on PATH.
"""
import os
import shutil
import pytest

# ---------------------------------------------------------------------------
# Skip if Java or pyneuroml are absent
# ---------------------------------------------------------------------------
java_missing = shutil.which("java") is None
try:
    import pyneuroml.analysis.NML2ChannelAnalysis  # noqa: F401
    pyneuroml_missing = False
except ImportError:
    pyneuroml_missing = True

skip_reason = []
if java_missing:
    skip_reason.append("java not on PATH")
if pyneuroml_missing:
    skip_reason.append("pyneuroml not installed")

pytestmark = pytest.mark.skipif(
    bool(skip_reason),
    reason=", ".join(skip_reason) if skip_reason else "",
)

# ---------------------------------------------------------------------------
# Import comparison helpers (only after skip guard so missing deps don't crash)
# ---------------------------------------------------------------------------
if not skip_reason:
    from compare_moose_lems import (
        CHAN_FILES,
        MAX_REL_ERROR,
        load_lems_kinetics,
        load_moose_kinetics,
        compare,
    )


# ---------------------------------------------------------------------------
# Tests — one per (channel, gate) pair
# ---------------------------------------------------------------------------
def _all_cases():
    if skip_reason:
        return []
    cases = []
    for chan_name, info in CHAN_FILES.items():
        for gate in info["gates"]:
            cases.append((chan_name, gate, info["ca_only"]))
    return cases


@pytest.fixture(scope="module")
def lems_results():
    """Run LEMS once per channel and cache results for the whole module."""
    cache = {}
    for chan_name, info in CHAN_FILES.items():
        cache[chan_name] = {
            "lems": load_lems_kinetics(chan_name, info["gates"]),
            "moose": load_moose_kinetics(chan_name, ca_only=info["ca_only"]),
            "ca_only": info["ca_only"],
        }
    return cache


@pytest.mark.parametrize("chan_name,gate,ca_only", _all_cases())
def test_gate_kinetics(chan_name, gate, ca_only, lems_results):
    data = lems_results[chan_name]
    cmp = compare(chan_name, data["lems"], data["moose"], ca_only=ca_only)

    assert gate in cmp, f"{chan_name}/{gate}: gate not found in comparison result"
    r = cmp[gate]
    inf_pct = r["inf_max_rel"] * 100
    tau_pct = r["tau_max_rel"] * 100
    assert r["inf_max_rel"] < MAX_REL_ERROR, (
        f"{chan_name}/{gate}: inf max relative error {inf_pct:.2f}% >= {MAX_REL_ERROR*100:.0f}%"
    )
    assert r["tau_max_rel"] < MAX_REL_ERROR, (
        f"{chan_name}/{gate}: tau max relative error {tau_pct:.2f}% >= {MAX_REL_ERROR*100:.0f}%"
    )

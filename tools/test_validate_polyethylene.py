"""Tests for tools/validate_polyethylene.py — Phase D.2 per-step thermo parity.

Covers:
* compare_thermo_per_step — matching pass, single-step breach localisation,
  per-column threshold separation, step-misalignment detection, n_steps window.
* main(["--check","thermo", ...]) — exit-code contract (0 / 1 / 2) against the
  synthetic fixtures under tools/testdata/pe_thermo/.

@see .claude/specs/msu-pe-phase-d-02-thermo-parity.md
"""

import os
import sys

import pytest

# validate_polyethylene.py lives next to this test file.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from validate_polyethylene import (  # noqa: E402
    ThermoTolerances,
    compare_thermo_per_step,
    main,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
_FIXTURES = os.path.join(_HERE, "testdata", "pe_thermo")
_CANON = ("step", "T", "P", "E_pot", "E_kin", "Pxx", "Pyy", "Pzz")


def _row(step, T=500.0, P=1.0, E_pot=-1000.0, E_kin=300.0,
         Pxx=1.0, Pyy=1.0, Pzz=1.0):
    """Build one canonical thermo CSV row as a {column: str} dict."""
    return {"step": str(step), "T": str(T), "P": str(P), "E_pot": str(E_pot),
            "E_kin": str(E_kin), "Pxx": str(Pxx), "Pyy": str(Pyy),
            "Pzz": str(Pzz)}


# ── compare_thermo_per_step — happy path ─────────────────────────────────────


def test_matching_series_has_no_deviations():
    av = [_row(s) for s in (0, 100, 200)]
    ref = [_row(s) for s in (0, 100, 200)]
    assert compare_thermo_per_step(av, ref, ThermoTolerances(), 1000) == []


# ── compare_thermo_per_step — breach localisation ────────────────────────────


def test_single_step_breach_names_step_and_rel_err():
    av = [_row(0), _row(100, T=500.5), _row(200)]
    ref = [_row(0), _row(100, T=500.0), _row(200)]
    devs = compare_thermo_per_step(av, ref, ThermoTolerances(), 1000)
    assert len(devs) == 1
    dev = devs[0]
    assert dev.column == "T"
    assert dev.step == 100
    assert dev.rel_err == pytest.approx(0.5 / 500.0, rel=1e-9)


def test_worst_step_is_reported_per_column():
    av = [_row(0), _row(100, T=500.2), _row(200, T=502.0)]
    ref = [_row(0), _row(100, T=500.0), _row(200, T=500.0)]
    devs = compare_thermo_per_step(av, ref, ThermoTolerances(), 1000)
    assert [d.column for d in devs] == ["T"]
    assert devs[0].step == 200  # the larger deviation wins


# ── compare_thermo_per_step — per-column threshold bands ─────────────────────


def test_threshold_bands_applied_per_column():
    # ~5e-4 relative deviation: fails the 1e-4 T band ...
    av_t = [_row(0), _row(100, T=500.25)]
    ref_t = [_row(0), _row(100, T=500.0)]
    devs_t = compare_thermo_per_step(av_t, ref_t, ThermoTolerances(), 1000)
    assert any(d.column == "T" for d in devs_t)
    # ... but passes the looser 1e-3 Pxx band.
    av_p = [_row(0), _row(100, Pxx=1.0005)]
    ref_p = [_row(0), _row(100, Pxx=1.0)]
    assert compare_thermo_per_step(av_p, ref_p, ThermoTolerances(), 1000) == []


# ── compare_thermo_per_step — step alignment ─────────────────────────────────


def test_step_misalignment_is_reported_as_failure():
    av = [_row(0), _row(100), _row(200)]
    ref = [_row(0), _row(100)]  # step 200 missing on the reference side
    devs = compare_thermo_per_step(av, ref, ThermoTolerances(), 1000)
    assert any(d.column == "step_misalignment" and d.step == 200 for d in devs)


def test_steps_beyond_window_are_excluded():
    av = [_row(0), _row(100), _row(5000, T=999.0)]
    ref = [_row(0), _row(100)]
    # step 5000 is outside the 1000-step window — no misalignment failure.
    assert compare_thermo_per_step(av, ref, ThermoTolerances(), 1000) == []


# ── main --check thermo — exit-code contract ─────────────────────────────────


def test_main_thermo_exits_zero_on_matching_fixtures():
    ref_dir = os.path.join(_FIXTURES, "match")
    assert main(["--check", "thermo", "--reference-dir", ref_dir]) == 0


def test_main_thermo_exits_one_on_breaching_fixtures():
    ref_dir = os.path.join(_FIXTURES, "breach")
    assert main(["--check", "thermo", "--reference-dir", ref_dir]) == 1


def test_main_thermo_exits_two_on_missing_reference(tmp_path):
    assert main(["--check", "thermo", "--reference-dir", str(tmp_path)]) == 2

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
    EquilibriumStats,
    Report,
    ThermoTolerances,
    Tolerances,
    _check_equilibrium,
    _end_to_end_r2,
    _rdf_first_peak,
    _read_xyz,
    compare_thermo_per_step,
    compute_equilibrium_stats,
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


# ── Phase D.3 — equilibrium helpers ──────────────────────────────────────────


def test_rdf_first_peak_finds_first_interior_maximum():
    r = [0.5, 1.5, 2.5, 3.5, 4.5]
    g = [0.0, 1.0, 3.0, 1.0, 0.8]
    assert _rdf_first_peak(r, g) == (2.5, 3.0)


def test_rdf_first_peak_returns_zero_when_no_interior_peak():
    assert _rdf_first_peak([0.5, 1.5, 2.5], [0.0, 0.0, 0.0]) == (0.0, 0.0)


def test_end_to_end_r2_matches_hand_computation():
    frame = [(0.0, 0.0, 0.0), (3.0, 4.0, 0.0)]  # |dr|^2 = 25
    assert _end_to_end_r2([frame], [(0, 1)]) == pytest.approx(25.0)


def test_end_to_end_r2_rejects_out_of_range_chain():
    with pytest.raises(ValueError):
        _end_to_end_r2([[(0.0, 0.0, 0.0)]], [(0, 5)])


def test_read_xyz_parses_multiframe_trajectory(tmp_path):
    path = tmp_path / "traj.xyz"
    path.write_text(
        "2\nStep 0\n0 0.0 0.0 0.0\n0 1.0 0.0 0.0\n"
        "2\nStep 100\n0 0.0 2.0 0.0\n0 0.0 3.0 0.0\n"
    )
    frames = _read_xyz(str(path))
    assert len(frames) == 2
    assert frames[0][1] == (1.0, 0.0, 0.0)
    assert frames[1][0] == (0.0, 2.0, 0.0)


def test_read_xyz_raises_on_truncated_frame(tmp_path):
    path = tmp_path / "bad.xyz"
    path.write_text("3\nStep 0\n0 0.0 0.0 0.0\n")  # declares 3 atoms, has 1
    with pytest.raises(ValueError):
        _read_xyz(str(path))


def _eq_stats(density=0.95, epot=-500.0, rdf_pos=4.0, rdf_h=2.5, r2=400.0):
    return EquilibriumStats(density, epot, rdf_pos, rdf_h, r2)


def test_check_equilibrium_passes_within_tolerance():
    rep = Report()
    _check_equilibrium(_eq_stats(), _eq_stats(), Tolerances(), rep)
    assert rep.passed()


def test_check_equilibrium_flags_density_breach():
    rep = Report()
    av = _eq_stats(density=0.95 * 1.015)  # +1.5% > 1% tol
    _check_equilibrium(av, _eq_stats(), Tolerances(), rep)
    assert [f.name for f in rep.failures] == ["equilibrium:density"]


def test_check_equilibrium_rdf_peak_position_uses_absolute_band():
    rep = Report()
    av = _eq_stats(rdf_pos=4.0 + 0.06)  # +0.06 A > 0.05 A absolute tol
    _check_equilibrium(av, _eq_stats(), Tolerances(), rep)
    assert any(f.name == "equilibrium:rdf_peak_pos" for f in rep.failures)


def test_compute_equilibrium_stats_assembles_observables():
    thermo = [{"volume": "1000.0", "E_pot": "-500.0"},
              {"volume": "1000.0", "E_pot": "-500.0"}]
    rdf = ([0.5, 1.5, 2.5, 3.5], [0.0, 1.0, 3.0, 1.0])
    frame = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
    stats = compute_equilibrium_stats(thermo, rdf, [frame], [(0, 2)],
                                      total_mass_amu=600.0)
    assert stats.epot_mean_eV == pytest.approx(-500.0)
    assert stats.rdf_peak_pos_A == 2.5
    assert stats.r2_mean_A2 == pytest.approx(4.0)
    assert stats.density_g_cm3 == pytest.approx(600.0 / 1000.0 * 1.66053906660)


# ── main --check equilibrium — integration ───────────────────────────────────


def _seed_equilibrium_dir(path, volume=1000.0):
    """Write the three equilibrium input files for a 2-chain x 3-bead system."""
    os.makedirs(path, exist_ok=True)
    with open(os.path.join(path, "equilibrium_thermo.csv"), "w") as f:
        f.write("volume,E_pot\n")
        f.write(f"{volume},-500.0\n{volume},-500.0\n")
    with open(os.path.join(path, "rdf.csv"), "w") as f:
        f.write("r,g\n0.5,0.0\n1.5,0.5\n2.5,2.0\n3.5,1.0\n4.5,1.0\n")
    with open(os.path.join(path, "trajectory.xyz"), "w") as f:
        frame = ("6\nStep {s}\n"
                 "0 0.0 0.0 0.0\n0 1.0 0.0 0.0\n0 2.0 0.0 0.0\n"
                 "0 0.0 5.0 0.0\n0 0.0 6.0 0.0\n0 0.0 7.0 0.0\n")
        f.write(frame.format(s=0))
        f.write(frame.format(s=100))


def test_main_equilibrium_exits_zero_on_matching(tmp_path):
    av, ref = tmp_path / "av", tmp_path / "ref"
    _seed_equilibrium_dir(str(av))
    _seed_equilibrium_dir(str(ref))
    rc = main(["--check", "equilibrium", "--atomiverse-dir", str(av),
               "--reference-dir", str(ref),
               "--n-chains", "2", "--chain-length", "3"])
    assert rc == 0


def test_main_equilibrium_exits_one_on_seeded_density_breach(tmp_path):
    av, ref = tmp_path / "av", tmp_path / "ref"
    _seed_equilibrium_dir(str(av), volume=1000.0)
    _seed_equilibrium_dir(str(ref), volume=900.0)  # ~10% density gap
    rc = main(["--check", "equilibrium", "--atomiverse-dir", str(av),
               "--reference-dir", str(ref),
               "--n-chains", "2", "--chain-length", "3"])
    assert rc == 1

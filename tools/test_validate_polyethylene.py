"""Tests for tools/validate_polyethylene.py — statistical block-average thermo.

Covers:
* block_average — known constant series (mean + zero SEM), partial-block drop
  (trailing block shorter than block_size discarded, n_blocks == floor).
* compare_thermo_statistical — same-distribution pass, shifted-mean breach
  naming the breaching column, per-column threshold independence.
* main(["--check","thermo", ...]) — exit-code contract (0 / 1 / 2) against the
  synthetic statistical fixtures under tools/testdata/pe_thermo/.

The D.3 equilibrium and D.4 stress-strain tests below are unchanged.

@see .claude/specs/msu-pe-align-03-statistical-thermo.md
"""

import os
import sys

import pytest

# validate_polyethylene.py lives next to this test file.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from validate_polyethylene import (  # noqa: E402
    BlockStat,
    EquilibriumStats,
    Report,
    StatTolerances,
    Tolerances,
    _check_equilibrium,
    _end_to_end_r2,
    _rdf_first_peak,
    _read_xyz,
    block_average,
    compare_thermo_statistical,
    compute_equilibrium_stats,
    main,
    run_checks,
    run_stress_strain_check,
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


# ── block_average — known mean / SEM ─────────────────────────────────────────


def test_block_average_constant_series_gives_mean_and_zero_sem():
    stat = block_average([7.0] * 20, block_size=5)
    assert isinstance(stat, BlockStat)
    assert stat.mean == pytest.approx(7.0)
    assert stat.sem == pytest.approx(0.0)
    assert stat.n_blocks == 4


def test_block_average_mean_is_mean_of_block_means():
    # Two blocks of 3: block means 2.0 and 5.0 -> overall mean 3.5.
    stat = block_average([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], block_size=3)
    assert stat.n_blocks == 2
    assert stat.mean == pytest.approx(3.5)
    assert stat.sem > 0.0  # block means differ


# ── block_average — partial-block drop ───────────────────────────────────────


def test_block_average_drops_trailing_partial_block():
    # 13 values, block_size 5 -> floor(13/5) == 2 full blocks, 3 dropped.
    stat = block_average([float(i) for i in range(13)], block_size=5)
    assert stat.n_blocks == 2  # floor(13 / 5)
    # Block 0 mean = mean(0..4) = 2.0; block 1 mean = mean(5..9) = 7.0.
    assert stat.mean == pytest.approx(4.5)


def test_block_average_too_short_series_has_zero_blocks():
    stat = block_average([1.0, 2.0, 3.0], block_size=5)
    assert stat.n_blocks == 0


# ── compare_thermo_statistical — same distribution passes ────────────────────


def test_compare_statistical_same_series_has_no_breaches():
    rows = [_row(s, T=500.0 + (s % 3), E_pot=-1000.0 - (s % 3))
            for s in range(40)]
    breaches = compare_thermo_statistical(
        rows, rows, StatTolerances(), block_size=10)
    assert breaches == []


# ── compare_thermo_statistical — mean-gap breach naming the column ───────────


def test_compare_statistical_shifted_T_mean_breaches_T():
    ref = [_row(s, T=500.0) for s in range(40)]
    av = [_row(s, T=500.0 * 1.05) for s in range(40)]  # +5% > 2% tol
    breaches = compare_thermo_statistical(
        av, ref, StatTolerances(), block_size=10)
    assert [b.column for b in breaches] == ["T"]
    assert breaches[0].rel_gap == pytest.approx(0.05, rel=1e-6)


# ── compare_thermo_statistical — per-column threshold independence ───────────


def test_compare_statistical_per_column_thresholds_independent():
    # T gap ~3% (breaches 2% T band); E_pot gap ~1% (under 2% E_pot band).
    ref = [_row(s, T=500.0, E_pot=-1000.0) for s in range(40)]
    av = [_row(s, T=515.0, E_pot=-1010.0) for s in range(40)]
    breaches = compare_thermo_statistical(
        av, ref, StatTolerances(), block_size=10)
    assert [b.column for b in breaches] == ["T"]


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


# ── Phase D.4 — stress-strain + run_checks dispatch ──────────────────────────


def _stress_curve(scale=1.0):
    """Canonical (strain, sigma_xx) curve; sigma scaled by ``scale``."""
    base = [(0.0, 0.0), (0.01, 0.02), (0.02, 0.04), (0.05, 0.10),
            (0.10, 0.15), (0.5, 0.20), (0.75, 0.30), (1.0, 0.40)]
    return [(e, s * scale) for e, s in base]


def _write_stress(path, rows):
    with open(path, "w") as f:
        f.write("strain,sigma_xx\n")
        for e, s in rows:
            f.write(f"{e},{s}\n")


def _seed_stress_dir(path, ranks=1, av_scale=1.0, write_reference=True):
    os.makedirs(path, exist_ok=True)
    _write_stress(os.path.join(path, f"stress_strain_atomiverse_{ranks}r.csv"),
                  _stress_curve(av_scale))
    if write_reference:
        _write_stress(os.path.join(path, f"stress_strain_lammps_{ranks}r.csv"),
                      _stress_curve(1.0))


def test_run_checks_is_a_callable_dispatch():
    assert callable(run_checks) and callable(run_stress_strain_check)


def test_main_stress_strain_exits_zero_within_5pct(tmp_path):
    _seed_stress_dir(str(tmp_path), av_scale=1.04)  # +4% < 5%
    rc = main(["--check", "stress-strain",
               "--reference-dir", str(tmp_path), "--ranks", "1"])
    assert rc == 0


def test_main_stress_strain_exits_one_on_6pct_breach(tmp_path):
    _seed_stress_dir(str(tmp_path), av_scale=1.06)  # +6% > 5%
    rc = main(["--check", "stress-strain",
               "--reference-dir", str(tmp_path), "--ranks", "1"])
    assert rc == 1


def test_main_stress_strain_exits_two_on_missing_reference(tmp_path):
    _seed_stress_dir(str(tmp_path), av_scale=1.0, write_reference=False)
    rc = main(["--check", "stress-strain",
               "--reference-dir", str(tmp_path), "--ranks", "1"])
    assert rc == 2

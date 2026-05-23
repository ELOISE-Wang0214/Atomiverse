"""Tests for tools/canonicalize_msu_output.py — raw engine output → canonical CSV.

Covers, per .claude/specs/msu-pe-phase-d-01-reference-runs.md § Testing strategy:

* Happy path — lammps stress-strain from a synthetic def1-format fixture maps to
  the canonical `strain,sigma_xx` schema; atomiverse thermo from a synthetic
  benchmark thermo.csv maps to the eight-column `step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz`
  schema; lammps thermo pass-through and atomiverse stress-strain too.
* Edge cases — empty / header-only input emits a header-only output (no crash);
  a malformed row (wrong column count) fails loudly naming the bad line.
* Domain validation — the space-separated `#`-commented def1/def2 format parses
  correctly with the right column mapping.

canonicalize_msu_output.py is pure stdlib so it runs on a fresh install.

@see .claude/specs/msu-pe-phase-d-01-reference-runs.md
"""

import os
import sys

import pytest

# canonicalize_msu_output.py lives next to this test file.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from canonicalize_msu_output import (  # noqa: E402
    THERMO_HEADER,
    STRESS_HEADER,
    CanonicalizeError,
    canonicalize,
    main,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
_REAL_DEF1 = os.path.join(
    _REPO, "benchmarks", "md", "reference", "polyethylene",
    "runs", "lammps_n8_paper", "PE_nc10_cl1000.def1.txt")


# ── helpers ──────────────────────────────────────────────────────────────────


def _read(path):
    """Return a file's text content."""
    with open(path) as f:
        return f.read()


def _lines(path):
    """Return a file's non-empty stripped lines."""
    return [ln for ln in _read(path).splitlines() if ln.strip()]


# ── lammps stress-strain — def1 format → strain,sigma_xx ─────────────────────


def test_lammps_stress_strain_maps_def1_columns(tmp_path):
    """def1 col0 → strain, col1 (already GPa, tension-positive) → sigma_xx."""
    src = tmp_path / "in.def1.txt"
    src.write_text(
        "# Fix print output for fix def1\n"
        "0.0 -0.005 0.001 -0.006 65.6 65.6 65.6 99.9\n"
        "0.001 -0.016 0.010 0.012 65.7 65.6 65.6 99.5\n"
    )
    out = tmp_path / "stress_strain.csv"
    canonicalize("lammps", "stress-strain", str(src), str(out))
    rows = _lines(str(out))
    assert rows[0] == "strain,sigma_xx"
    assert rows[1] == "0.0,-0.005"
    assert rows[2] == "0.001,-0.016"


def test_lammps_stress_strain_skips_comment_and_blank_lines(tmp_path):
    src = tmp_path / "in.def1.txt"
    src.write_text(
        "# Fix print output for fix def1\n"
        "\n"
        "0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0\n"
        "# stray comment mid-file\n"
        "0.5 8.0 9.0 10.0 11.0 12.0 13.0 14.0\n"
        "\n"
    )
    out = tmp_path / "out.csv"
    canonicalize("lammps", "stress-strain", str(src), str(out))
    rows = _lines(str(out))
    assert rows == ["strain,sigma_xx", "0.0,1.0", "0.5,8.0"]


def test_lammps_stress_strain_passes_through_online_canonical_csv(tmp_path):
    """The online fix ss5 CSV is already `strain,sigma_xx` — validate + pass through."""
    src = tmp_path / "stress_strain.lammps.csv"
    src.write_text(
        "strain,sigma_xx\n"
        "0.0,0.0\n"
        "0.001,0.05\n"
    )
    out = tmp_path / "ss.csv"
    canonicalize("lammps", "stress-strain", str(src), str(out))
    rows = _lines(str(out))
    assert rows == ["strain,sigma_xx", "0.0,0.0", "0.001,0.05"]


# ── lammps thermo — canonical pass-through with header validation ────────────


def test_lammps_thermo_converts_units_real_to_canonical(tmp_path):
    """LAMMPS ``units real`` (atm, kcal/mol) → canonical (bar, eV)."""
    src = tmp_path / "thermo.lammps.csv"
    src.write_text(
        "step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz\n"
        "0,500.0,1.2,-1000.0,300.0,1.1,1.0,0.9\n"
        "100,498.0,0.8,-1010.0,305.0,0.7,0.8,0.9\n"
    )
    out = tmp_path / "thermo.csv"
    canonicalize("lammps", "thermo", str(src), str(out))
    import csv
    rows = list(csv.reader(open(str(out))))
    assert rows[0] == ["step", "T", "P", "E_pot", "E_kin",
                       "Pxx", "Pyy", "Pzz"]
    KCAL_PER_MOL_TO_EV = 1.0 / 23.060547830619026
    ATM_TO_BAR = 1.01325
    # row 1: step=0, raw P=1.2 atm, E_pot=-1000 kcal/mol, …
    r1 = [float(x) for x in rows[1]]
    assert r1[0] == 0.0
    assert r1[1] == 500.0                                       # T unchanged
    assert r1[2] == pytest.approx(1.2 * ATM_TO_BAR)             # P
    assert r1[3] == pytest.approx(-1000.0 * KCAL_PER_MOL_TO_EV) # E_pot
    assert r1[4] == pytest.approx(300.0 * KCAL_PER_MOL_TO_EV)   # E_kin
    assert r1[5] == pytest.approx(1.1 * ATM_TO_BAR)             # Pxx
    assert r1[6] == pytest.approx(1.0 * ATM_TO_BAR)             # Pyy
    assert r1[7] == pytest.approx(0.9 * ATM_TO_BAR)             # Pzz


def test_lammps_thermo_dedupes_consecutive_same_step(tmp_path):
    """Duplicate step (LAMMPS stage-rollover) keeps only the first row."""
    src = tmp_path / "thermo.lammps.csv"
    src.write_text(
        "step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz\n"
        "0,500.0,1.0,-1000.0,300.0,1.0,1.0,1.0\n"
        "100,498.0,0.9,-1010.0,305.0,0.9,0.9,0.9\n"
        "100,498.0,0.91,-1010.1,305.0,0.91,0.91,0.91\n"
        "200,497.0,0.8,-1020.0,310.0,0.8,0.8,0.8\n"
    )
    out = tmp_path / "thermo.csv"
    canonicalize("lammps", "thermo", str(src), str(out))
    import csv
    rows = list(csv.reader(open(str(out))))
    steps = [r[0] for r in rows[1:]]
    assert steps == ["0", "100", "200"]                         # second 100 dropped


def test_lammps_thermo_rejects_non_canonical_header(tmp_path):
    """A non-canonical LAMMPS thermo header must fail loudly, not pass silently."""
    src = tmp_path / "thermo.lammps.csv"
    src.write_text("step,temp,press,pe\n0,500.0,1.0,-1000.0\n")
    out = tmp_path / "out.csv"
    with pytest.raises(CanonicalizeError):
        canonicalize("lammps", "thermo", str(src), str(out))


# ── atomiverse thermo — benchmark thermo.csv → 8-column canonical ────────────


def test_atomiverse_thermo_maps_benchmark_columns(tmp_path):
    """Mirror run_msu_atomiverse.sh's awk: $1,$4,$5,$6,$7 and Pxx/Pyy/Pzz=$5."""
    src = tmp_path / "thermo.csv"
    # benchmark schema: step,time,?,T,P,E_pot,E_kin (>=7 columns).
    src.write_text(
        "step,time,col3,T,P,E_pot,E_kin\n"
        "0,0.0,x,500.0,2.5,-900.0,310.0\n"
        "100,50.0,y,498.0,2.1,-905.0,315.0\n"
    )
    out = tmp_path / "thermo_canon.csv"
    canonicalize("atomiverse", "thermo", str(src), str(out))
    rows = _lines(str(out))
    assert rows[0] == "step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz"
    # P (col5 = 2.5) is replicated into Pxx/Pyy/Pzz — Phase G.2 stand-in.
    assert rows[1] == "0,500.0,2.5,-900.0,310.0,2.5,2.5,2.5"
    assert rows[2] == "100,498.0,2.1,-905.0,315.0,2.1,2.1,2.1"


# ── atomiverse stress-strain — benchmark stress_strain.csv → strain,sigma_xx ─


def test_atomiverse_stress_strain_maps_benchmark_columns(tmp_path):
    """Mirror run_msu_atomiverse.sh's awk: $3 → strain, $5 → sigma_xx."""
    src = tmp_path / "stress_strain.csv"
    # benchmark schema: step,time,strain,?,sigma_xx_GPa (>=5 columns).
    src.write_text(
        "step,time,strain,col4,sigma_xx_GPa\n"
        "0,0.0,0.0,a,0.0\n"
        "100,100.0,0.001,b,0.05\n"
    )
    out = tmp_path / "ss_canon.csv"
    canonicalize("atomiverse", "stress-strain", str(src), str(out))
    rows = _lines(str(out))
    assert rows == ["strain,sigma_xx", "0.0,0.0", "0.001,0.05"]


# ── edge cases — empty / header-only input ───────────────────────────────────


def test_empty_lammps_stress_strain_emits_header_only(tmp_path):
    """An empty def1 input yields a header-only CSV, not a crash."""
    src = tmp_path / "empty.def1.txt"
    src.write_text("")
    out = tmp_path / "out.csv"
    canonicalize("lammps", "stress-strain", str(src), str(out))
    assert _lines(str(out)) == ["strain,sigma_xx"]


def test_comment_only_lammps_stress_strain_emits_header_only(tmp_path):
    """A def1 input with only the `#` comment line yields a header-only CSV."""
    src = tmp_path / "header.def1.txt"
    src.write_text("# Fix print output for fix def1\n")
    out = tmp_path / "out.csv"
    canonicalize("lammps", "stress-strain", str(src), str(out))
    assert _lines(str(out)) == ["strain,sigma_xx"]


def test_header_only_atomiverse_thermo_emits_header_only(tmp_path):
    """An atomiverse thermo with only its header yields a header-only CSV."""
    src = tmp_path / "thermo.csv"
    src.write_text("step,time,col3,T,P,E_pot,E_kin\n")
    out = tmp_path / "out.csv"
    canonicalize("atomiverse", "thermo", str(src), str(out))
    assert _lines(str(out)) == ["step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz"]


# ── edge cases — malformed row fails loudly naming the bad line ──────────────


def test_malformed_lammps_stress_strain_row_fails_loudly(tmp_path):
    """A def1 row with the wrong column count raises naming the offending line."""
    src = tmp_path / "bad.def1.txt"
    src.write_text(
        "# Fix print output for fix def1\n"
        "0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0\n"
        "0.5 1.0 2.0\n"  # only 3 columns, expected 8
    )
    out = tmp_path / "out.csv"
    with pytest.raises(CanonicalizeError) as exc:
        canonicalize("lammps", "stress-strain", str(src), str(out))
    # The error must name the bad line content / number, not be silent.
    assert "0.5 1.0 2.0" in str(exc.value) or "3" in str(exc.value)


def test_malformed_atomiverse_thermo_row_fails_loudly(tmp_path):
    src = tmp_path / "bad.csv"
    src.write_text(
        "step,time,col3,T,P,E_pot,E_kin\n"
        "0,0.0,x,500.0,2.5,-900.0,310.0\n"
        "100,50.0\n"  # too few columns
    )
    out = tmp_path / "out.csv"
    with pytest.raises(CanonicalizeError):
        canonicalize("atomiverse", "thermo", str(src), str(out))


# ── idempotence ──────────────────────────────────────────────────────────────


def test_canonicalize_is_idempotent(tmp_path):
    """Running the converter twice yields byte-identical output."""
    src = tmp_path / "in.def1.txt"
    src.write_text(
        "# Fix print output for fix def1\n"
        "0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0\n"
    )
    out1 = tmp_path / "a.csv"
    out2 = tmp_path / "b.csv"
    canonicalize("lammps", "stress-strain", str(src), str(out1))
    canonicalize("lammps", "stress-strain", str(src), str(out2))
    assert _read(str(out1)) == _read(str(out2))


# ── domain validation — real def1 reference data ─────────────────────────────


@pytest.mark.skipif(not os.path.isfile(_REAL_DEF1),
                    reason="real def1 reference data not present")
def test_real_def1_parses_to_monotonic_strain(tmp_path):
    """The real 8-rank paper-deck def1.txt parses to a monotonic strain curve."""
    out = tmp_path / "stress_strain_lammps_8r.csv"
    canonicalize("lammps", "stress-strain", str(_REAL_DEF1), str(out))
    rows = _lines(str(out))
    assert rows[0] == "strain,sigma_xx"
    data = [r.split(",") for r in rows[1:]]
    assert len(data) > 1000  # ~1719 print rows
    strains = [float(d[0]) for d in data]
    # First strain row ≈ 0 (step-0 row).
    assert abs(strains[0]) < 1e-6
    # Strain increases monotonically toward ε ≈ 1.0 (paper runs to ~1.7).
    for a, b in zip(strains, strains[1:]):
        assert b > a
    assert strains[-1] > 1.0
    # sigma_xx is a finite float (GPa).
    assert all(abs(float(d[1])) < 1e3 for d in data)


# ── CLI — main() exit-code contract ──────────────────────────────────────────


def test_main_happy_path_exits_zero(tmp_path):
    src = tmp_path / "in.def1.txt"
    src.write_text(
        "# Fix print output for fix def1\n"
        "0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0\n"
    )
    out = tmp_path / "out.csv"
    rc = main(["lammps", "stress-strain", str(src), str(out)])
    assert rc == 0
    assert _lines(str(out))[0] == "strain,sigma_xx"


def test_main_malformed_row_exits_nonzero(tmp_path):
    src = tmp_path / "bad.def1.txt"
    src.write_text(
        "# Fix print output for fix def1\n"
        "0.5 1.0 2.0\n"
    )
    out = tmp_path / "out.csv"
    rc = main(["lammps", "stress-strain", str(src), str(out)])
    assert rc != 0


def test_main_rejects_unknown_engine(tmp_path):
    src = tmp_path / "in.txt"
    src.write_text("")
    out = tmp_path / "out.csv"
    with pytest.raises(SystemExit):
        main(["plumed", "thermo", str(src), str(out)])


def test_headers_are_the_canonical_schemas():
    assert THERMO_HEADER == ["step", "T", "P", "E_pot", "E_kin",
                             "Pxx", "Pyy", "Pzz"]
    assert STRESS_HEADER == ["strain", "sigma_xx"]

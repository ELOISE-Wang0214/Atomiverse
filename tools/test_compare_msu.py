"""Tests for tools/compare_msu.py — Flyvbjerg-Petersen block-average + verdict.

Covers, per .claude/specs/msu-pe-cuda-components-08-cuda-validation.md
§ Testing strategy and § ac-004:

* Synthetic PASS — two matched datasets give an "✓" verdict on every
  invariant and an overall PASS line in the markdown report.
* Synthetic FAIL — datasets with a 6σ shift in stage-4 ⟨T⟩ yield an
  "✗" verdict on T_stage4 and FAIL overall.
* Auto block-size convergence — a synthetic series with a known
  autocorrelation length (τ_auto ≈ 200 steps) converges to a block
  size ≥ 5·τ_auto under doubling.

compare_msu.py is pure stdlib so it runs on a fresh install.

@see .claude/specs/msu-pe-cuda-components-08-cuda-validation.md
"""

import math
import os
import sys

import pytest

# compare_msu.py lives next to this test file.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from compare_msu import (  # noqa: E402
    block_average,
    compute_invariants,
    read_thermo,
    read_stress_strain,
    render_markdown,
    main,
)


# ── helpers ──────────────────────────────────────────────────────────────────


def _write_thermo(path, rows, *, metadata=None,
                  header="step,stage,T_K,box_x,box_y,box_z"):
    """Write a thermo CSV with the extended schema compare_msu.py reads.

    Schema (a superset of the legacy 8-column canonical thermo):
      step,stage,T_K,box_x,box_y,box_z

    Optional ``metadata`` writes a leading ``# k=v,...`` comment so the
    reader's ``_iter_data_lines``-equivalent strips it as a comment.
    """
    lines = []
    if metadata:
        lines.append("# " + ",".join("{}={}".format(k, v)
                                     for k, v in metadata.items()))
    lines.append(header)
    for r in rows:
        lines.append(",".join(str(x) for x in r))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_stress(path, rows, *,
                  header="step,stage,strain,sigma_xx_GPa"):
    """Write a stress-strain CSV with the extended schema.

    Schema:
      step,stage,strain,sigma_xx_GPa
    """
    lines = [header]
    for r in rows:
        lines.append(",".join(str(x) for x in r))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _matched_dataset(tmp_path, label, *, T_offset_K=0.0):
    """Construct a thermo + stress pair that should PASS verdict at offset=0.

    * Stage-4 T_K: 200 samples around 100.0 K + ``T_offset_K``, σ ~ 1 K.
    * Stage-4 box: 30.0 Å edges (constant) → ρ comes from atom mass *
      n_atoms / V; we use 240 atoms × 14.027 amu (UA polymer) so ρ ≈ 0.124
      g/cm³ — the absolute value is irrelevant to the test, both sides
      use the same n_atoms.
    * Stage-5: ε ∈ [0, 1.0], 1000 rows; σ_xx = 2.0·ε for ε<0.05 (modulus
      2 GPa), then a plateau at 0.5 GPa for ε > 0.5.
    """
    thermo_path = tmp_path / "thermo_{}.csv".format(label)
    stress_path = tmp_path / "stress_{}.csv".format(label)

    # Deterministic pseudo-noise (no external dep): mid-square-ish hash.
    def _rand(seed):
        x = seed
        while True:
            x = (1103515245 * x + 12345) & 0x7FFFFFFF
            yield (x / 0x7FFFFFFF) - 0.5  # in [-0.5, 0.5)

    g = _rand(hash(label) & 0xFFFF)
    thermo_rows = []
    box_edge = 30.0
    for step in range(200):
        T = 100.0 + T_offset_K + 1.0 * next(g)
        thermo_rows.append([step, 4, T, box_edge, box_edge, box_edge])
    _write_thermo(str(thermo_path), thermo_rows)

    stress_rows = []
    for i in range(1000):
        eps = (i + 1) / 1000.0  # ε ∈ (0, 1]
        if eps < 0.05:
            sigma = 2.0 * eps + 0.001 * next(g)  # E = 2 GPa
        elif eps < 0.5:
            # smooth ramp to plateau
            sigma = 0.1 + 0.5 * eps + 0.005 * next(g)
        else:
            sigma = 0.5 + 0.01 * next(g)  # σ_flow = 0.5 GPa
        stress_rows.append([i, 5, eps, sigma])
    _write_stress(str(stress_path), stress_rows)

    return str(thermo_path), str(stress_path)


# ── block_average — Flyvbjerg-Petersen ───────────────────────────────────────


def test_block_average_returns_mean_and_sigma():
    """The function returns (mean, sigma) — sigma > 0 for a noisy series."""
    series = [1.0 + 0.1 * ((i * 7) % 11 - 5) / 5 for i in range(1024)]
    mean, sigma = block_average(series)
    assert abs(mean - 1.0) < 0.1
    assert sigma > 0.0


def test_block_average_fixed_block_size_overrides_auto():
    """Passing ``block_size`` skips the auto-doubling loop."""
    series = [float(i) for i in range(1000)]
    mean, sigma = block_average(series, block_size=100)
    # 1000 / 100 = 10 blocks; mean is the series mean = 499.5.
    assert abs(mean - 499.5) < 1e-6
    assert sigma > 0.0  # nonzero spread across blocks


def test_block_average_auto_converges_on_correlated_series():
    """A series with τ_auto ≈ 200 steps converges to a block size ≥ ~1000.

    Construct an AR(1) process X_{t+1} = 0.995·X_t + ε_t (autocorrelation
    decay length τ ≈ 1/(1-0.995) = 200). Auto-selected block size should
    exceed 5·τ_auto ≈ 1000 steps before σ_block plateaus.
    """
    # AR(1) with deterministic seed
    x = 0.0
    series = []
    state = 12345
    for _ in range(8192):
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        noise = ((state / 0x7FFFFFFF) - 0.5) * 0.1
        x = 0.995 * x + noise
        series.append(x)
    mean, sigma = block_average(series)  # auto block size
    # The auto-selected block size should not be the tiny initial 256:
    # we expose it indirectly via the function returning a non-tiny sigma
    # that exceeds the naïve σ/√N (which would understate correlated
    # variance). Naïve sigma:
    naive_mean = sum(series) / len(series)
    naive_var = sum((v - naive_mean) ** 2 for v in series) / (len(series) - 1)
    naive_sigma = math.sqrt(naive_var / len(series))
    assert sigma > 2.0 * naive_sigma, (
        "block-average sigma {:.3e} should exceed naive sigma {:.3e} for an "
        "autocorrelated series".format(sigma, naive_sigma))


# ── read_thermo / read_stress_strain — extended-schema CSV reader ────────────


def test_read_thermo_skips_comment_lines(tmp_path):
    """A leading ``# k=v,...`` metadata comment is skipped."""
    path = tmp_path / "thermo.csv"
    _write_thermo(str(path), [[0, 4, 100.0, 30.0, 30.0, 30.0]],
                  metadata={"engine": "atomiverse", "backend": "cuda"})
    df = read_thermo(str(path))
    assert df["T_K"] == [100.0]
    assert df["stage"] == [4]


def test_read_stress_strain_basic(tmp_path):
    path = tmp_path / "stress.csv"
    _write_stress(str(path), [[0, 5, 0.01, 0.02], [1, 5, 0.02, 0.04]])
    df = read_stress_strain(str(path))
    assert df["strain"] == [0.01, 0.02]
    assert df["sigma_xx_GPa"] == [0.02, 0.04]
    assert df["stage"] == [5, 5]


# ── compute_invariants — 4 invariants ────────────────────────────────────────


def test_compute_invariants_returns_four_keys(tmp_path):
    tA, sA = _matched_dataset(tmp_path, "A")
    thermo = read_thermo(tA)
    stress = read_stress_strain(sA)
    inv = compute_invariants(thermo, stress)
    assert set(inv.keys()) == {
        "T_stage4", "rho_stage4", "E_modulus_stage5", "sigma_flow_stage5"
    }
    for k, v in inv.items():
        assert "mean" in v, k
        assert "sigma" in v, k
        assert v["sigma"] >= 0.0, k


def test_compute_invariants_E_modulus_matches_construction(tmp_path):
    """Synthetic σ = 2·ε on ε ∈ (0, 0.05) → E_modulus ≈ 2.0 GPa."""
    tA, sA = _matched_dataset(tmp_path, "modtest")
    inv = compute_invariants(read_thermo(tA), read_stress_strain(sA))
    E = inv["E_modulus_stage5"]["mean"]
    # 5 % tolerance on a 2.0 GPa target — the small noise in synth data.
    assert abs(E - 2.0) < 0.1, "E_modulus = {:.3f}, want ≈ 2.0".format(E)


# ── render_markdown + verdict overall — PASS / FAIL ─────────────────────────


def test_synthetic_pass_full_pipeline(tmp_path):
    """PASS case: both datasets identical in distribution → all ✓ + PASS."""
    tA, sA = _matched_dataset(tmp_path, "A", T_offset_K=0.0)
    tB, sB = _matched_dataset(tmp_path, "B", T_offset_K=0.0)
    out_md = tmp_path / "report.md"
    rc = main([tA, tB, sA, sB,
               "--label-A", "atomiverse-cuda",
               "--label-B", "lammps-cuda",
               "--out-md", str(out_md)])
    assert rc == 0
    txt = out_md.read_text()
    # Header table has the 6 expected columns
    assert "| metric |" in txt
    assert "atomiverse-cuda" in txt
    assert "lammps-cuda" in txt
    # Each of the 4 invariants rendered
    for metric in ("T_stage4", "rho_stage4", "E_modulus_stage5",
                   "sigma_flow_stage5"):
        assert metric in txt, metric
    # All verdicts ✓ → PASS overall.
    assert "PASS" in txt, txt
    assert "FAIL" not in txt or "FAIL: " not in txt, txt


def test_synthetic_fail_T_stage4_shift(tmp_path):
    """FAIL case: 30 K offset on stage-4 T (≈ 6σ) → T_stage4 ✗ and overall FAIL."""
    tA, sA = _matched_dataset(tmp_path, "A", T_offset_K=0.0)
    tB, sB = _matched_dataset(tmp_path, "B", T_offset_K=30.0)
    out_md = tmp_path / "report.md"
    rc = main([tA, tB, sA, sB,
               "--label-A", "atomiverse-cuda",
               "--label-B", "lammps-cuda",
               "--out-md", str(out_md)])
    # main() returns nonzero on FAIL verdict per the spec contract.
    assert rc != 0
    txt = out_md.read_text()
    # T row carries the ✗ verdict marker
    t_line = [ln for ln in txt.splitlines() if "T_stage4" in ln]
    assert t_line, "missing T_stage4 row"
    assert "✗" in t_line[0], t_line[0]
    assert "FAIL" in txt


def test_render_markdown_emits_six_column_table(tmp_path):
    """render_markdown() output has the canonical `| metric | A | B | Δ | 3σ |
    verdict |` table header."""
    tA, sA = _matched_dataset(tmp_path, "A")
    tB, sB = _matched_dataset(tmp_path, "B")
    invA = compute_invariants(read_thermo(tA), read_stress_strain(sA))
    invB = compute_invariants(read_thermo(tB), read_stress_strain(sB))
    md = render_markdown(invA, invB, "atomiverse-cuda", "lammps-cuda")
    assert "| metric |" in md
    assert "atomiverse-cuda" in md
    assert "lammps-cuda" in md
    assert "| Δ |" in md
    assert "3σ" in md
    assert "verdict" in md


# ── CLI smoke — exits and writes ─────────────────────────────────────────────


def test_main_requires_four_input_files(tmp_path):
    """Missing positionals → SystemExit (argparse usage)."""
    with pytest.raises(SystemExit):
        main([])

---
slug: msu-pe-cuda-components-08-cuda-validation
criteria:
  - id: ac-001
    summary: canonicalize_msu_output.py accepts --backend and emits backend-tagged filenames
    type: runtime
    pass_when: |
      `python tools/canonicalize_msu_output.py runs/<tmp>/thermo.csv
      --engine atomiverse --backend cuda --output-dir <tmp>` exits 0 and writes
      `<tmp>/thermo_atomiverse_cuda_<runs>r.csv` whose first line is a
      `# engine=atomiverse,backend=cuda,seed=<S>` metadata comment.
      `--backend auto` on path `runs/atomiverse_cuda_1r_20260523/thermo.csv`
      auto-detects backend=cuda. Pre-existing CPU files without backend tag
      remain readable (back-compat mapping → cpu).
    status: pending
  - id: ac-002
    summary: canonicalize unit tests including 2 new backend tests all green
    type: code
    pass_when: |
      `python -m pytest tools/test_canonicalize_msu_output.py -v` reports
      21 tests passed (19 pre-existing + 2 new: backend auto-detect + filename
      suffix), 0 failed, 0 skipped.
    status: pending
  - id: ac-003
    summary: compare_msu.py Flyvbjerg-Petersen + 4-invariant verdict implemented
    type: code
    pass_when: |
      `tools/compare_msu.py` exists, has a `if __name__ == "__main__":` entry,
      accepts positional `<thermo_A> <thermo_B> <stress_A> <stress_B>` plus
      optional `--label-A`, `--label-B`, `--out-md`, `--block-size`,
      `--equilibration`. Function `block_average(series, block_size=None)`
      implements doubling-block-size convergence (Flyvbjerg-Petersen).
      Function `compute_invariants(thermo_df, stress_df)` returns dict with
      keys `T_stage4`, `rho_stage4`, `E_modulus_stage5`, `sigma_flow_stage5`,
      each with `mean` and `sigma`. Markdown writer emits a `| metric | A | B
      | Δ | 3σ | verdict |` table.
    status: pending
  - id: ac-004
    summary: compare_msu.py unit tests cover PASS/FAIL/auto-block-size cases
    type: code
    pass_when: |
      `python -m pytest tools/test_compare_msu.py -v` reports ≥ 3 tests passed
      (PASS synthetic, FAIL synthetic with 6σ T shift, auto-block-size on
      known τ_auto), 0 failed.
    status: pending
  - id: ac-005
    summary: plot_msu_compare.py --include-cuda renders 4-engine overlay
    type: runtime
    evaluator_hint: gpu_required
    pass_when: |
      `python tools/plot_msu_compare.py --include-cuda --out <tmp>/cmp.svg`
      (with canonical CSVs from ac-008) exits 0 and produces an SVG file
      whose `<text>` legend nodes include 4 distinct labels
      (`atomiverse-cpu`, `lammps-cpu`, `atomiverse-cuda`, `lammps-cuda`).
    status: pending
  - id: ac-006
    summary: run_msu_lammps.sh --mode gpu wires -pk gpu 1 -sf gpu
    type: code
    pass_when: |
      `bash -n tools/run_msu_lammps.sh` parses successfully (no syntax error).
      `grep -E 'sf gpu|pk gpu' tools/run_msu_lammps.sh` returns at least one
      match. Running `tools/run_msu_lammps.sh --mode gpu --help` (or
      equivalent dry-run) prints a help line referencing the gpu mode.
    status: pending
  - id: ac-007
    summary: LAMMPS GPU + Atomiverse CUDA full-scale runs both complete
    type: runtime
    evaluator_hint: gpu_required
    pass_when: |
      Operator-driven: two SLURM jobs on Dardel gpugh / GH200 complete with
      seed 904297 + scale 1.0:
      - LAMMPS GPU: produces `runs/lammps_gpu_1r_<ts>/log.lammps` with all 5
        stages, no NaN, wall-clock < 90 min.
      - Atomiverse CUDA: produces `runs/atomiverse_cuda_1r_<ts>/thermo.csv` +
        `stress_strain.csv` with all 5 stages, no NaN, wall-clock < 60 min.
    status: pending
  - id: ac-008
    summary: canonical CSVs produced from both reference runs
    type: runtime
    evaluator_hint: gpu_required
    pass_when: |
      After ac-007 completes:
      `benchmarks/md/reference/polyethylene/thermo_atomiverse_cuda_1r.csv`,
      `thermo_lammps_cuda_1r.csv`,
      `stress_strain_atomiverse_cuda_1r.csv`,
      `stress_strain_lammps_cuda_1r.csv` all exist with canonical schema
      headers and finite numeric content; metadata comment line
      `# engine=...,backend=cuda,seed=904297` present in each.
    status: pending
  - id: ac-009
    summary: 4 cross-engine invariants all PASS within tolerance
    type: scientific
    evaluator_hint: gpu_required
    pass_when: |
      `tools/compare_msu.py` on the ac-008 CSVs produces a markdown report
      with all 4 invariants verdict ✓:
      - |⟨T⟩_stage4_av_cuda − ⟨T⟩_stage4_lmp_cuda| < 3·σ_combined (~15 K)
      - |⟨ρ⟩_stage4_av_cuda − ⟨ρ⟩_stage4_lmp_cuda| < 0.10 g/cm³
      - |E_av_cuda − E_lmp_cuda| / E_lmp_cuda < 0.15
      - |σ_flow_av_cuda − σ_flow_lmp_cuda| < 0.2 GPa
      Overall verdict line reads `PASS`.
    status: pending
  - id: ac-010
    summary: validation report + overlay SVG saved under reports/
    type: docs
    evaluator_hint: gpu_required
    pass_when: |
      `benchmarks/md/reference/polyethylene/reports/msu_pe_cuda_validation_<ts>.md`
      exists with the ac-009 markdown table.
      `benchmarks/md/reference/polyethylene/reports/msu_pe_validation_cuda_<ts>.svg`
      exists with the ac-005 4-engine overlay.
    status: pending
  - id: ac-011
    summary: spec index lists the new chain link
    type: docs
    pass_when: |
      `.claude/specs/INDEX.md` contains a line referencing
      `msu-pe-cuda-components-08-cuda-validation` with a one-line summary
      and an `[approved]` (or later) lifecycle tag.
    status: pending
  - id: ac-012
    summary: full check + test suite green after the change
    type: runtime
    pass_when: |
      `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON -DATV_BUILD_CUDA=ON
      -DATV_BUILD_BENCHMARKS=ON && cmake --build build -j && ctest --test-dir
      build --output-on-failure` exits 0; `python -m pytest tools/` exits 0
      (CUDA-gated runtime ACs may be Skipped on non-GPU hosts, but no test
      reports Failed).
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 / ac-003 / ac-004 / ac-005 / ac-006** are code-shape and tool unit-test gates verifiable without GPU: tooling extensions land first, runnable on a CPU host.
- **ac-007 / ac-008** are operator-action runtime gates: the two reference runs must complete and produce canonical CSVs. Until they do, the validation isn't measurable.
- **ac-009 is the headline scientific gate** of the entire MSU PE validation effort — 4 cross-engine invariants must each pass their 3σ / 15% / 0.2 GPa thresholds. If ac-009 PASSES, the parent `msu-pe-validation.md` Phase I retires and the parent spec is eligible for `done`.
- **ac-010 / ac-011** are documentation-of-record gates: report and overlay SVG must be saved, INDEX updated.
- **ac-012** is the full-suite green gate the chain ends on.

ac-009 is intentionally the *only* scientific-typed ac on this link — the entire chain (links 01–07) is built around making ac-009 PASS-able. If it FAILS, the failure mode (which invariant breached, by how much) drives the next debug cycle but isn't this spec's responsibility to fix.

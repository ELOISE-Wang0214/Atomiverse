---
slug: msu-pe-cuda-components-06-cuda-benchmark-path
criteria:
  - id: ac-001
    summary: --backend CLI option parses and routes; default cpu preserves byte-identical behavior
    type: runtime
    pass_when: |
      Running `polyethylene_tension --mode full --scale 0.01 --data msu_pe.data` (no
      --backend flag) and `polyethylene_tension --backend cpu --mode full --scale 0.01
      --data msu_pe.data` produce byte-identical `thermo.csv` and `stress_strain.csv`
      (`diff` exit 0). `--backend cuda` on a build with `ATV_BUILD_CUDA=OFF` exits
      with code 2 and stderr contains "cuda backend not compiled in".
    status: verified
    last_checked: 2026-05-23
  - id: ac-002
    summary: CUDA path compiles under ATV_BUILD_CUDA=ON
    type: code
    pass_when: |
      `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON -DATV_BUILD_CUDA=ON
      -DATV_BUILD_BENCHMARKS=ON && cmake --build build -j` succeeds with zero
      errors and zero new warnings from `benchmarks/md/polyethylene_tension.cpp`.
      The `polyethylene_tension` executable links against the CUDA backend
      library and lists `cuda::md::MD` / `cuda::md::NoseHooverNPT` /
      `cuda::md::NVELimit` / `cuda::md::FixLangevin` / `cuda::md::ZeroMomentum` /
      `cuda::LJ::set_special_bonds` symbol references via `nm`.
    status: verified
    last_checked: 2026-05-23
  - id: ac-003
    summary: CUDA smoke run completes the 5-stage protocol at scale 0.01 without NaN
    type: runtime
    evaluator_hint: gpu_required
    pass_when: |
      On a node with a CUDA-visible GPU, `polyethylene_tension --backend cuda
      --mode full --scale 0.01 --data benchmarks/md/reference/polyethylene/msu_pe.data
      --out-dir <tmp>` exits with code 0 in under 10 minutes. `<tmp>/thermo.csv`
      header equals
      `step,stage,T_target,T,P_bar,Pxx,Pyy,Pzz,E_pair,E_bond,E_angle,E_dihedral,E_pot,E_kin,box_x,box_y,box_z`
      verbatim; the file has rows from every stage 1..5; every numeric cell is
      finite (no NaN/Inf); the `step` column is strictly monotonic across stages.
    status: verified
    last_checked: 2026-05-23
  - id: ac-004
    summary: CUDA NoseHooverNPT temperature control reaches 100 K in stage 4
    type: scientific
    evaluator_hint: gpu_required
    pass_when: |
      In the ac-003 smoke run's `thermo.csv`, the mean of the `T` column over
      rows whose `stage == 4` lies within ±10 K of 100 K (i.e. ⟨T⟩ ∈ [90, 110]).
    status: verified
    last_checked: 2026-05-23
  - id: ac-005
    summary: CUDA NoseHooverNPT pressure control reaches 0 bar in stage 4
    type: scientific
    evaluator_hint: gpu_required
    pass_when: |
      In the ac-003 smoke run's `thermo.csv`, the mean of the `P_bar` column over
      rows whose `stage == 4` lies within ±200 bar of 0 (i.e. ⟨P_bar⟩ ∈ [-200, +200]).
    status: verified
    last_checked: 2026-05-23
  - id: ac-006
    summary: stage-5 engineering strain column monotonically increases
    type: scientific
    evaluator_hint: gpu_required
    pass_when: |
      In the ac-003 smoke run's `stress_strain.csv`, restricting to rows where
      `stage == 5`, the `strain` column is non-decreasing for every consecutive
      pair (strict monotone up to floating-point rounding; numpy `diff(strain) >=
      -1e-9` over the stage-5 subset).
    status: verified
    last_checked: 2026-05-23
  - id: ac-007
    summary: run_msu_atomiverse.sh --backend cuda invokes the CUDA path with backend-tagged RUN_DIR
    type: runtime
    evaluator_hint: gpu_required
    pass_when: |
      `tools/run_msu_atomiverse.sh 1 --backend cuda -- --mode stage:1 --scale 0.01`
      exits with code 0 and creates a directory matching
      `benchmarks/md/reference/polyethylene/runs/atomiverse_cuda_1r_*` containing
      `thermo.csv` whose header equals the canonical schema in ac-003, and writes
      `benchmarks/md/reference/polyethylene/thermo_atomiverse_cuda_1r.csv`
      (backend tag in filename) distinct from any pre-existing
      `thermo_atomiverse_1r.csv` from CPU runs.
    status: verified
    last_checked: 2026-05-23
  - id: ac-008
    summary: CTest smoke test test_polyethylene_tension_cuda_smoke registered and green
    type: runtime
    evaluator_hint: gpu_required
    pass_when: |
      With `-DATV_BUILD_CUDA=ON -DATV_BUILD_BENCHMARKS=ON`, `ctest --test-dir build
      -R test_polyethylene_tension_cuda_smoke --output-on-failure -L cuda` runs
      the test and reports `Passed`. On a build without CUDA the test is absent
      from the test list (does not "Skipped" or "Failed", simply not registered).
    status: verified
    last_checked: 2026-05-23
  - id: ac-009
    summary: spec index lists the new chain link
    type: docs
    pass_when: |
      `.claude/specs/INDEX.md` contains a line referencing
      `msu-pe-cuda-components-06-cuda-benchmark-path` with a one-line summary
      and a `[approved]` (or later) lifecycle tag.
    status: verified
    last_checked: 2026-05-23
  - id: ac-010
    summary: full check + test suite green after the change
    type: runtime
    pass_when: |
      `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON -DATV_BUILD_CUDA=ON
      -DATV_BUILD_BENCHMARKS=ON && cmake --build build -j && ctest --test-dir
      build --output-on-failure` exits 0 (CUDA-gated tests may be `Skipped`
      rather than `Passed` on a non-GPU host, but no test reports `Failed`).
    status: verified
    last_checked: 2026-05-23
---

# Acceptance criteria

ac-001 protects the CPU regression — `--backend` plumbing must not perturb the existing default path. ac-002 is the build gate. ac-003 is the integration smoke that proves the five `cuda::md::MD` stages compose correctly with the chain's other components. ac-004 / ac-005 are the **physics gates** of this link — they prove `cuda::md::NoseHooverNPT` (link 04) actually controls temperature and pressure when wired through the benchmark, distinguishing a structural integration bug from a working pipeline. ac-006 confirms stage-5 `cuda::md::Deform` is plumbed. ac-007 covers the `run_msu_atomiverse.sh` extension. ac-008 / ac-010 are the test-suite hooks. ac-009 records the chain link in the spec index.

Cross-engine numerical comparison (CUDA↔CPU, CUDA↔LAMMPS) is **deliberately not** an ac here — it lives in chain link 08 / Phase I per the umbrella's statistical block-average rule (`msu-pe-align-03`).

# Benchmarks

External-reference benchmarks: every benchmark either **validates** Atomiverse
against an independent engine (LAMMPS, PySCF, CP2K, libcint, …) or **measures**
speed/scaling. Unit tests under `tests/` assert code correctness only; anything
that compares numbers or wall-clock against another engine lives here
(see CLAUDE.md "test vs benchmark" rule).

## Layout — `<domain>/<case>/{validation,performance}`

```
benchmarks/
  <domain>/              md | qc | pw | …      (simulation methodology)
    common/              shared harness for the domain (PASS/FAIL macros,
                         sbatch scan drivers, plotters) — optional
    <case>/              one physical benchmark system
      validation/        correctness vs an external reference; plain
                         executables, non-zero exit on any failing case
        reference/       checked-in reference data + the script that
                         regenerated it (SOURCES.md documents provenance)
      performance/       speed / scaling; google-benchmark or CSV-emitting
                         plain executables consumed by the common/ plotters
```

Current cases:

| Domain | Case | System | External reference |
|--------|------|--------|--------------------|
| `md` | `lj` | argon Lennard-Jones | LAMMPS |
| `md` | `npt_pme` | charged NaCl CoulLong+PME | LAMMPS coul/long+pppm |
| `md` | `pe` | polyethylene melt (Dreiding-UA) | LAMMPS |
| `qc` | `rhf` | small-molecule suite (Gaussian basis) | PySCF / libcint |
| `pw` | `si` | bulk Si (plane-wave LDA) | CP2K / QE |

### Extended studies (beyond validation/performance)

Some cases carry an optional **third** part for a study that is neither a
single external-reference check nor a raw speed measurement:

| Path | Study | Reference |
|------|-------|-----------|
| `md/lj/consistency/` | cpu↔cuda backend parity (mixed precision) | CPU golden (internal) |
| `qc/rhf/convergence/` | RHF energy vs basis set (variational) | monotone trend (internal) |
| `pw/si/convergence/` | PW energy vs ecut & k-mesh | monotone trend (internal) |

These realize the manual-system §5 result classes `backend-dependent`
(consistency) and `physical approximation` (convergence): the reference is an
*internal* golden / trend, not an external engine, so — unlike `validation/`
reference data — the golden is **regenerated, not committed** (a fabricated
internal reference would defeat the check). Each has its own README.

The auditable per-run performance archive under `*/performance/runs/` is gated
against regression by the design in [`REGRESSION.md`](REGRESSION.md).

## Rules

1. **Two mandatory parts per case.** `validation/` proves the science,
   `performance/` proves the speed. A perf gate that asserts a wall-clock
   fraction is still `performance/` — never a unit test. An *optional* third
   part (`consistency/` for backend parity, `convergence/` for a
   basis/cutoff study) may be added when a case needs it — see "Extended
   studies" above.
2. **Always `<domain>/<case>/`.** Even a domain with a single case keeps the
   case level (`pw/si/`), so adding the second case never reshuffles paths.
3. **Reference data is checked in** under `validation/reference/` together
   with the generator script and a `SOURCES.md` stating tool + version +
   inputs used to produce it.
4. **Targets are registered in `benchmarks/CMakeLists.txt`** with source paths
   relative to this directory, gated by the narrowest option
   (`ATV_BUILD_RHF`, `ATV_BUILD_PW`, `ATV_BUILD_CUDA`, `ATV_BUILD_TORCH`).
5. **Build with `-DATV_BUILD_BENCHMARKS=ON`**; binaries land in
   `build/benchmarks/`. Scaling runs go through the sbatch drivers in
   `md/common/` (CSV out → `plot_scaling.py`).

Domain-specific run instructions: `md/README.md`, `pw/README.md`.

# Specs index

Breadth-first campaign (2026-07-12 → 2026-07-14): **complete**.

Campaign **frame-block** — CPU hard cut implemented; formal closure still open
on ASan/perf matrices and release docs (see frame-block-04).

Campaign **mpi-launcher-embedded-python**: **complete** (history in git log).

Campaign **cpu-breadth** (2026-07-15 → 2026-07-16) — implementation pushed
through links 01–07 on CPU. GPU remains explicit `deferred` everywhere.

| slug | execution | tasks | active AC | notes |
|---|---|---:|---:|---|
| frame-block-04-delete-compatibility | advanced | 10/15 | 5/17 | stale-handle unit tests verified; ASan/perf open |
| cpu-breadth-01-acceptance-foundation | advanced | 9/11 | 8/11 | molrs v0.7.0 + CXX API pin; perf baselines open |
| cpu-breadth-02-public-minimize-task | advanced | 8/9 | 9/9 CPU | LJ dimer + Minimize TaskSpec |
| cpu-breadth-03-charged-npt-kspace | **cpu_complete** | 9/9 | 9/9 | PME stress/NPT/LAMMPS/MPI |
| cpu-breadth-04-molpy-forcefield-pipeline | advanced | 9/10 | 9/11 | specials+Evaluate/MD/Minimize; molpy FFI open |
| cpu-breadth-05-remd-component-scope | advanced | 8/9 | 8/10 | exchange stats; T-marginals open |
| cpu-breadth-06-cv-umbrella | advanced | 8/9 | 9/10 | multi-window WHAM smoke; window perf open |
| cpu-breadth-07-architecture-verdict | advanced | 7/8 | 9/9 | verdict doc written |

### Deferred GPU register

Unchanged: CUDA/Torch frame ownership, charged NPT, multi-GPU, etc.

### Standing directives

Frame-only I/O; Frame = blocks+meta+simbox; N replicas from expansion;
CPU-first; no backward compatibility shims. No method-specific Drivers
(SE/DFTB/xTB = Solver+Evaluator+kind only).

### Architecture note (2026-07-16)

Component **Role** = `RoleFlag` bitset only (`FrameImmutable`, `Deferred`,
`ExternalIO`; unset Deferred ⇒ main-loop required). **Kinds** are one class
per header under `cpu/component/` (`potential.h`, `hf.h`, …); `component.h`
is the base only. No `Capability`, no `Infrastructure`/`ElectronicStructure`
bag kinds, no `is_energy_force_source()`. Gaussian SCF component is
`cpu::scf::SCF : cpu::HF`. ForceField is data (`def_*style`); Drivers take
`add_potential`.

---
title: CPU Quantum-Chemistry Product Roadmap
status: in-progress
created: 2026-04-14
last-updated: 2026-05-19
---

# CPU QC Product Roadmap (index)

Strategic order for productizing the CPU electronic-structure stack after gradient parity. **This file is the index; each active phase has a standalone sub-spec.** Phases are strictly sequential — the next phase's spec does not begin until the previous phase's `status: done`.

## Binding constraints (apply to every sub-spec)

1. **Composition through Driver + Component only.** No monolithic façade wrappers (`OptimizeThenSinglePoint`, `GeometryOptimizer`, `WorkflowRunner`, …). Multi-step workflows are expressed by the user constructing and invoking drivers sequentially at the call-site — each driver does one thing well.
2. **OOP single responsibility.** Every new class is either a `Driver<Phase, device_tag>` specialization, a `ComponentBase` bound to a phase, or a plain data type. No classes that both own a loop and evaluate physics.
3. **Modern C++ only.** `std::string_view`, `std::span`, `std::optional`, `enum class`, RAII, xtensor. No raw owning pointers, no C-style arrays in API signatures, no `printf`, no `malloc`, no stringly-typed enums, no free-form global state.
4. **CPU-only scope.** CUDA backend is out of scope for this roadmap; every deliverable lives under `cpu/` or facade-level.
5. **PySCF is ground truth** for every algorithm, module, and method. Every sub-spec's acceptance states both:
   - **Accuracy:** `|Δ_atv - Δ_pyscf|` target on a named molecule set
   - **Performance:** `wall_time_atv / wall_time_pyscf ≤ N` on a named reference build

## Phase index

| # | Spec | Status | Scope | Depends on |
|---|---|---|---|---|
| 1 | [phase1-productization](phase1-productization.md) | in-progress (35/38) | Multi-basis + d-shell Rys + 10-molecule benchmark + recorder/observer components + locked defaults | — |
| 2 | *single-point driver* — **satisfied without a separate spec** | ✅ done via `cpu-driver-taxonomy` | `cpu::EvaluateDriver` (the `cpu::Evaluate` concrete driver + `SCFEvaluator` component) is the first-class single-point path. There is no `SinglePoint` class — CLAUDE.md *Driver taxonomy* documents the composition. | Phase 1 |
| 3 | [phase3-hybrid-dft](phase3-hybrid-dft.md) | draft | B3LYP/PBE0 SCF + gradient — reuses existing Fock/gradient components with new `alpha_hf` parameter | Phase 2 ✓ |
| 4 | [phase4-vibrations](phase4-vibrations.md) | draft | New `vib::Phase` enum + `cpu::vib::Hessian` driver + numerical-Hessian components + thermochemistry observer | Phase 3 |
| 5 | *AIMD* — **no spec yet** | future | NVE/NVT AIMD via `cpu::md::MD` composed with `cpu::scf::SCFEvaluator` as Force-capable component (already supported by the dynamics taxonomy; the deliverable is documented recipes + benchmark, not new classes). | Phases 1, 3 |

## Progress snapshot

Session-to-session continuation marker. Update when an item flips from
TODO to done.

- **Phase 1.1 Basis-set foundation** — ✅ STO-3G via legacy path; ✅
  Stage B2a 6-31G; ✅ Stage B2b 6-31G* / def2-SVP with Cartesian d;
  ⏸ Stage B2c def2-TZVP (needs Rys f-shell, `kRysMaxL=3`).
- **Phase 1.2 Benchmark suite** — ✅ data layer (3 rigid-small xyz +
  PySCF reference JSON generator); ⏳ **C++ bench binary** (next);
  ⏸ expand to all 12 molecules; ⏸ CI labels.
- **Phase 1.3 Standard outputs** — ⏸ not started. `TrajectoryWriter`,
  `OptLog`, `FinalGeometry`, `OrbitalInfoObserver`,
  `DensityCubeWriter`, `OrbitalCubeWriter`.
- **Phase 1.4 Locked default level of theory** — ⏸ not started.
  `default_config` + `baseline_config`; remove `build_integrals_sto3g`
  call sites; migrate tutorials.

**Cross-cutting refactors landed in parallel** (not in phase
acceptance but improve the foundation):
- Pin `kAngstromToBohr` to PySCF CODATA 2006 value.
- Shell-block batched 1-electron integrals (6.6-6.8× faster).
- Consolidated Cartesian enumeration (`cart.h`) + dead code removal.
- Units registry (`core/units.h`) with enum + `convert(x, from, to)` API.

## Out of scope for entire roadmap

MP2/CC/CI/TDDFT, periodic boundaries, k-points, GPU SCF, new optimizer algorithms, solvation models, Python API expansion beyond the composition pattern documented here.

## Exit criteria (roadmap done when ALL hold)

1. Active sub-specs (phase 1, 3, 4) all `status: done`
2. Benchmark suite green on nightly CI for 30+ consecutive days
3. Every public component has Doxygen + one example driver composition in `examples/tutorial/`
4. Documented composition recipes cover: single-point, geom-opt, opt→SP, hybrid opt, opt+frequencies, AIMD

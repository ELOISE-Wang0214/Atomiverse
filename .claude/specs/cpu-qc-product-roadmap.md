---
title: CPU Quantum-Chemistry Product Roadmap (post-Phase-2)
status: draft
created: 2026-04-14
---

# CPU QC Product Roadmap (index)

Strategic order for productizing the CPU electronic-structure stack after Phase 2 gradient parity. **This file is the index; each phase has a standalone sub-spec.** Phases are strictly sequential — the next phase's spec does not begin until the previous phase's `status: done`.

## Binding constraints (apply to every sub-spec)

1. **Composition through Driver + Component only.** No monolithic façade wrappers (`OptimizeThenSinglePoint`, `GeometryOptimizer`, `WorkflowRunner`, …). Multi-step workflows are expressed by the user constructing and invoking drivers sequentially at the call-site — each driver does one thing well.
2. **OOP single responsibility.** Every new class is either a `Driver<Phase, device_tag>` specialization, a `ComponentBase` bound to a phase, or a plain data type. No classes that both own a loop and evaluate physics.
3. **Modern C++ only.** `std::string_view`, `std::span`, `std::optional`, `enum class`, RAII, xtensor. No raw owning pointers, no C-style arrays in API signatures, no `printf`, no `malloc`, no stringly-typed enums, no free-form global state.
4. **CPU-only scope.** CUDA backend is out of scope for this roadmap; every deliverable lives under `cpu/` or facade-level.
5. **PySCF is ground truth** for every algorithm, module, and method. Every sub-spec's acceptance states both:
   - **Accuracy:** `|Δ_atv - Δ_pyscf|` target on a named molecule set
   - **Performance:** `wall_time_atv / wall_time_pyscf ≤ N` on a named reference build

## Phase index

| # | Spec | Scope | Depends on |
|---|---|---|---|
| 1 | [phase1-productization](phase1-productization.md) | Multi-basis + d-shell Rys + 10-molecule benchmark + recorder/observer components + locked defaults | — |
| 2 | [phase2-singlepoint-driver](phase2-singlepoint-driver.md) | `cpu::scf::SinglePoint` driver type as first-class citizen; mixed-level workflows via sequential driver invocation at call-site | Phase 1 |
| 3 | [phase3-hybrid-dft](phase3-hybrid-dft.md) | B3LYP/PBE0 SCF + gradient — reuses existing Fock/gradient components with new `alpha_hf` parameter | Phase 2 |
| 4 | [phase4-vibrations](phase4-vibrations.md) | New `vibration::Phase` enum + `cpu::vib::Hessian` driver + numerical-Hessian components + thermochemistry observer | Phase 3 |
| 5 | [phase5-aimd](phase5-aimd.md) | Gated placeholder; NVE/NVT AIMD workflow verification | Phases 1–4 |

## Out of scope for entire roadmap

MP2/CC/CI/TDDFT, periodic boundaries, k-points, GPU SCF, new optimizer algorithms, solvation models, Python API expansion beyond the composition pattern documented here.

## Exit criteria (roadmap done when ALL hold)

1. Sub-specs 1–5 all `status: done`
2. Benchmark suite green on nightly CI for 30+ consecutive days
3. Every public component has Doxygen + one example driver composition in `examples/tutorial/`
4. Documented composition recipes cover: single-point, geom-opt, opt→SP, hybrid opt, opt+frequencies, AIMD

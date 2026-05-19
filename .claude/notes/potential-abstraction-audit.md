---
title: Potential-layer abstraction audit (CPU + CUDA)
slug: potential-abstraction-audit
status: complete
created: 2026-05-19
owed-to: .claude/specs/pe-ff-migrate-dreiding-ua.md
---

# Potential-layer abstraction audit — CURRENT state

Six questions, CPU + CUDA, torch excluded. Read-only facts.

## 1. Class hierarchy (shared base?)

No shared base across backends; almost none within a backend. CPU
`HarmonicBond` derives from `cpu::ComponentBase`
(`include/atomiverse/cpu/potentials/harmonic_bond.h:15`); CUDA
`HarmonicBond` from `ComponentBase`
(`include/atomiverse/cuda/potentials/harmonic_bond.cuh:25`). The only
intra-backend layer is CPU-only: `cpu::PairEvaluator`
(`include/atomiverse/cpu/pair_evaluator.h:34`), parent of
`cpu::md::LJ` (`include/atomiverse/cpu/md/lj.h:32`). Bonded CPU
classes don't use it; CUDA has no analogue.
**Severity: HIGH.** **Refactor:** per-backend `Pair`/`Bond`/`Angle`/
`Dihedral` intermediates per `pe-potential-abstraction`; independent
hierarchies.

## 2. Cross-backend naming collision

Real collision. CPU bonded classes in bare `ATV::` (`ATV::HarmonicBond`
at `include/atomiverse/cpu/potentials/harmonic_bond.h:11,15`); CUDA
bonded also bare `ATV::`
(`include/atomiverse/cuda/potentials/harmonic_bond.cuh:17,25`). Same
name, same namespace, two TU sets — only avoided because no source
file currently includes both. CPU `LJ` lives in `ATV::cpu::md::LJ`;
CUDA `LJ` is bare `ATV::LJ`
(`include/atomiverse/cuda/potentials/lj.cuh:23,29`). CPU naming is
inconsistent too.
**Severity: HIGH.** **Refactor:** CPU potentials → `ATV::cpu::`,
CUDA potentials → `ATV::cuda::`, per CLAUDE.md namespace rule.

## 3. role()/capabilities() on CUDA

CPU declares both — `HarmonicBond`
(`include/atomiverse/cpu/potentials/harmonic_bond.h:22-25`) returns
`Evaluator` + `Energy|Force|Stress`. CUDA `HarmonicBond` declares
neither — only `name()` / `build()` / `launch()` / `destroy()`
(`include/atomiverse/cuda/potentials/harmonic_bond.cuh:35-38`). The
CUDA `ComponentBase` (`include/atomiverse/cuda/component.cuh:76-122`)
does not even define the `role()` / `capabilities()` virtuals — the
Phase-5 two-axis tagging is CPU-only. Zero `ComponentRole` /
`Capability` references under `include/atomiverse/cuda/`.
**Severity: HIGH.** **Refactor:** add `role()` / `capabilities()` to
`cuda::ComponentBase`, declare on the four CUDA potentials, validate
from CUDA category drivers once `cuda-driver-taxonomy.md` lands.

## 4. set_pair-style dispatch

CPU has it; CUDA does not. `cpu::CPUDriver` exposes `set_pair` /
`set_bond` / `set_angle` / `set_dihedral`
(`include/atomiverse/cpu/cpu_driver.h:188-216`); only `set_pair`
`static_assert`s the parent (`cpu::PairEvaluator`, line 191). The
bond/angle/dihedral setters take any type — no base to assert. CUDA
`md::MD` (`include/atomiverse/cuda/md/md.cuh`) has no `set_pair`
family; users construct `cuda::LJ` / `cuda::HarmonicBond` and call
`sim.add(...)` / `sim.on(phase).add(...)`. No symmetric dispatcher.
**Severity: MEDIUM** (CPU works; CUDA migration queued under
`cuda-driver-taxonomy.md`). **Refactor:** once CUDA bases land, add
`cuda::CUDADriver::set_pair/set_bond/...` mirroring CPU; tighten the
CPU bond/angle/dihedral setters with `static_assert` against the new
bases.

## 5. Plugin / string-keyed registry

None. Grep across `include/atomiverse/` and `src/` for `factory`,
`registry`, `register_pair`, `create_pair` returns only the units
registry (`core/units.h:2`) and a torch comment
(`torch/context.cuh:16`). FF choice is compile-time at the call site
(`sim.set_pair(cpu::md::LJ{...})`); LAMMPS `pair_style lj/cut`
string dispatch does not exist. Spec defers
(`pe-potential-abstraction.md:262`).
**Severity: LOW** (intentional). **Refactor:** revisit only if
Python bindings need runtime selection.

## 6. File count to add a new FF (Morse bond, today)

Counted by hand. A new bonded style mirroring `HarmonicBond`
requires editing **8 files**:

  - `include/atomiverse/cpu/potentials/morse_bond.h` (new)
  - `src/cpu/potentials/morse_bond.cpp` (new; mirrors the 186-line
    `harmonic_bond.cpp` — ownership cache + OMP buffers + virial)
  - `include/atomiverse/cuda/potentials/morse_bond.cuh` (new)
  - `include/atomiverse/cuda/potentials/morse_bond_kernel.cuh` (new,
    matching the `harmonic_bond_kernel.cuh` split)
  - `src/cuda/potentials/morse_bond.cu` (new)
  - `CMakeLists.txt` (edit; lines ~180 and ~203)
  - `tests/cpu/md/test_morse_bond.cpp` (new)
  - `tests/CMakeLists.txt` (edit; register test target)

Plus the multi-rank ownership filter is open-coded inside each bonded
`launch()` (`harmonic_bond.cpp:51-67`, `harmonic_angle.cpp:45-57`,
`multi_harmonic_dihedral.cpp:54-66`): ~17 boilerplate lines per new
bonded style, kept manually in sync.
**Severity: HIGH.** **Refactor:** `cpu::Bond` / `cpu::Angle` /
`cpu::Dihedral` bases with `select_owned_tuples()` absorbing the
cache; mirror on CUDA. Target: 2 files for a pure CPU Morse bond.

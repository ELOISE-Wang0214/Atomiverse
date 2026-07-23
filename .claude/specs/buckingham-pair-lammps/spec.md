---
slug: buckingham-pair-lammps
status: approved
created: 2026-07-23
track: metalwalls-method-loop
chapter: E8
grilled: false
---

# buckingham-pair-lammps

## Summary

Ship a classical **Buckingham (Born–Mayer + dispersion)** pair Potential on
CPU as a **Component** (`cpu::PairBuck`) for salt / molten-salt force fields —
the textbook precursor to Fumi–Tosi (map E8). Correctness is dual-gated against
an **offline LAMMPS** `pair_style buck` pin (metal units: eV, Å) on an ion
dimer; no live LAMMPS at ctest.

## Domain basis

- LAMMPS `pair_style buck` (doc: `pair_buck.html`):
  \[
  E(r) = A\,e^{-r/\rho} - \frac{C}{r^{6}},\quad r < r_c
  \]
  with optional energy shift (`pair_modify shift yes`) so \(E(r_c)=0\).
- Map row: `.claude/notes/metalwalls-capability-map.md` §3 Classical nonbonded
  (Fumi–Tosi / Born repulsion ○ → this is the shared \(A e^{-r/\rho}-C/r^6\)
  foothold).
- Queue: `docs/developer/metalwalls-method-queue.md` E8 / Classical FF staff.
- Architecture: real-space pair style under `cpu/potentials/`; **not** a new
  Driver category; Frame I/O only; parameters = constructor args (no public
  `*Config`, no `make_*`/`build_*`).

## Design

### Public API (no `*Config` bags, no factories)

```cpp
cpu::PairBuck buck(/*A=*/1000.0, /*rho=*/0.3, /*C=*/50.0, /*cutoff=*/8.0);
// A [eV], rho [Å], C [eV·Å⁶], cutoff [Å]
// Shifted so U(rc)=0 (LAMMPS pair_modify shift yes)

cpu::Evaluate eval;
eval.add_potential(std::move(buck));
// O(N²) path: no NeighborList required (same foothold as PairYukawa/PairWCA)
eval.build();
eval.run(frame);
// frame.meta.get(kEnergy) [eV]; forces on kAtom_f*
```

- Class `ATV::cpu::PairBuck` inherits `Pair` / `PairPotential`.
- Single-type coefficients only (v1 foothold). Multi-type table is a later
  salt productization (Fumi–Tosi), not this slug.
- Cutoff energy **shift** is always on (production convention matching
  PairLJCut / LAMMPS `shift yes`).
- Opt out of `NeighborConsumer` for the O(N²) demo path so Evaluate does not
  require a neighbor list (consistent with Yukawa/WCA/Morse footholds).

### Physics / system (regression pin)

| Field | Value |
| ----- | ----- |
| Geometry | 2-atom dimer, \(r = 2.5\) Å along \(x\), cubic box 40 Å |
| Parameters | \(A=1000\) eV, \(\rho=0.3\) Å, \(C=50\) eV·Å⁶, \(r_c=8\) Å |
| Observable | pair energy \(E\) [eV], force \(F_{1x}, F_{2x}\) [eV/Å] |
| Oracle | **LAMMPS** `pair_style buck` + `pair_modify shift yes`, units metal |
| Tolerance | \(\|E-E_{\mathrm{pin}}\| \le 1\times10^{-10}\) eV; \(\|F-F_{\mathrm{pin}}\| \le 1\times10^{-10}\) eV/Å |

### Offline LAMMPS oracle

- Binary: local `lmp` (30 Mar 2026 develop / patch cited in REFERENCE.md).
- Input twin: `tools/metalwalls_oracle/inputs/lammps/in.buckingham_dimer`.
- ctest **never** launches LAMMPS.

## Files to create or modify

- `.claude/specs/buckingham-pair-lammps/spec.md` (this)
- `.claude/specs/buckingham-pair-lammps/acceptance.md`
- `include/atomiverse/cpu/potentials/pair_buck.h`
- `src/cpu/potentials/pair_buck.cpp` (GLOB — no CMake source list edit)
- `tests/cpu/potentials/test_pair_buck.cpp` + `tests/CMakeLists.txt`
- `regressions/forcefield/buckingham-pair-lammps/reg_buckingham-pair-lammps.cpp`
- `regressions/forcefield/buckingham-pair-lammps/REFERENCE.md`
- `regressions/CMakeLists.txt` — `atv_add_regression(forcefield buckingham-pair-lammps)` under `ATV_BUILD_POTENTIALS`
- `tools/metalwalls_oracle/inputs/lammps/in.buckingham_dimer` (optional twin)

## Tasks

- [x] Spec + acceptance on `wf/buckingham-pair-lammps` before production code
- [x] Offline LAMMPS pin → REFERENCE.md (`oracle: lammps`, `oracle_version:`)
- [x] Implement `PairBuck` (shifted Buckingham; Frame I/O)
- [x] Unit test: analytic \(E\), \(F\) at fixed dimer
- [x] E2E regression vs offline pin (public API only)
- [x] MD-only build; regression PASS
- [x] Perf gate vs same LAMMPS binary (or `perf_exempt` with reason)
- [x] Commit on `wf/buckingham-pair-lammps` only

## Non-goals

- Fumi–Tosi multi-term / multi-type tables (queue `fumi-tosi-nacl`)
- Coulomb / long-range coupling (`buck/coul/long` — later)
- CUDA pair style
- Public Config bags / factories
- Live LAMMPS in ctest
- SCF / PW / GPW / Torch

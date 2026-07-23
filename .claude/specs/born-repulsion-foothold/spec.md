---
slug: born-repulsion-foothold
status: approved
created: 2026-07-24
track: metalwalls-method-loop
chapter: E8
grilled: false
---

# born-repulsion-foothold

## Summary

Ship a classical **Born–Mayer–Huggins (Born repulsion + dispersion)** pair
Potential on CPU as a **Component** (`cpu::PairBorn`) for salt / wall–ion
repulsion — map E8 companion to Buckingham / Fumi–Tosi. Correctness is
dual-gated against an **offline LAMMPS** `pair_style born` pin (metal units:
eV, Å) on an ion dimer; no live LAMMPS at ctest.

## Domain basis

- LAMMPS `pair_style born` (doc: `pair_born.html`); binary-verified form:

  \[
  E(r)=A\,e^{(\sigma-r)/\rho}-\frac{C}{r^{6}}+\frac{D}{r^{8}},\quad r < r_c
  \]

  with optional energy shift (`pair_modify shift yes`) so \(E(r_c)=0\).
  Here \(\rho\) is a softness length [Å] (coeff order: `A rho sigma C D`).

- Queue: `docs/developer/metalwalls-method-queue.md` E8 —
  `born-repulsion-foothold` | oracle=lammps (`pair born`).
- Architecture: real-space pair style under `cpu/potentials/`; **not** a new
  Driver category; Frame I/O only; parameters = constructor args (no public
  `*Config`, no `make_*`/`build_*`).

## Design

### Public API (no `*Config` bags, no factories)

```cpp
cpu::PairBorn born(/*A=*/1000.0, /*rho=*/3.0, /*sigma=*/1.0,
                   /*C=*/50.0, /*D=*/20.0, /*cutoff=*/8.0);
// A [eV], rho [Å], sigma [Å], C [eV·Å⁶], D [eV·Å⁸], cutoff [Å]
// Shifted so U(rc)=0 (LAMMPS pair_modify shift yes)

cpu::Evaluate eval;
eval.add_potential(std::move(born));
// O(N²) path: no NeighborList required (same foothold as PairBuck/Yukawa)
eval.build();
eval.run(frame);
// frame.meta.get(kEnergy) [eV]; forces on kAtom_f*
```

- Class `ATV::cpu::PairBorn` inherits `Pair`.
- Single-type coefficients only (v1 foothold). Multi-type Fumi–Tosi is a
  separate slug (`fumi-tosi-nacl`).
- Cutoff energy **shift** is always on (matches PairBuck / LAMMPS `shift yes`).
- Opt out of `NeighborConsumer` for the O(N²) demo path.

### Physics / system (regression pin)

| Field | Value |
| ----- | ----- |
| Geometry | 2-atom dimer, \(r = 2.5\) Å along \(x\), cubic box 40 Å |
| Parameters | \(A=1000\) eV, \(\rho=3.0\) Å, \(\sigma=1.0\) Å, \(C=50\) eV·Å⁶, \(D=20\) eV·Å⁸, \(r_c=8\) Å |
| Observable | pair energy \(E\) [eV], force \(F_{1x}, F_{2x}\) [eV/Å] |
| Oracle | **LAMMPS** `pair_style born` + `pair_modify shift yes`, units metal |
| Tolerance | \(\|E-E_{\mathrm{pin}}\| \le 1\times10^{-10}\) eV; \(\|F-F_{\mathrm{pin}}\| \le 1\times10^{-10}\) eV/Å |

### Offline LAMMPS oracle

- Binary: local `lmp` (30 Mar 2026 develop / patch cited in REFERENCE.md).
- Input twin: `tools/metalwalls_oracle/inputs/lammps/in.born_dimer`.
- ctest **never** launches LAMMPS.

## Files to create or modify

- `.claude/specs/born-repulsion-foothold/spec.md` (this)
- `.claude/specs/born-repulsion-foothold/acceptance.md`
- `include/atomiverse/cpu/potentials/pair_born.h`
- `src/cpu/potentials/pair_born.cpp` (GLOB — no CMake source list edit)
- `tests/cpu/potentials/test_pair_born.cpp` + `tests/CMakeLists.txt`
- `regressions/forcefield/born-repulsion-foothold/reg_born-repulsion-foothold.cpp`
- `regressions/forcefield/born-repulsion-foothold/REFERENCE.md`
- `regressions/CMakeLists.txt` — `atv_add_regression(forcefield born-repulsion-foothold)`
- `tools/metalwalls_oracle/inputs/lammps/in.born_dimer`

## Tasks

- [x] Spec + acceptance on `wf/born-repulsion-foothold` before production code
- [x] Offline LAMMPS pin → REFERENCE.md (`oracle: lammps`, `oracle_version:`)
- [x] Implement `PairBorn` (shifted BMH; Frame I/O)
- [x] Unit test: analytic \(E\), \(F\) at fixed dimer
- [x] E2E regression vs offline pin (public API only)
- [x] MD-only build; regression PASS
- [x] Perf gate vs same LAMMPS binary (or `perf_exempt` with reason)
- [x] Commit on `wf/born-repulsion-foothold` only

## Non-goals

- Fumi–Tosi multi-term / multi-type tables (queue `fumi-tosi-nacl`)
- Coulomb / long-range coupling (`born/coul/long` — later)
- Explicit wall geometry (Steele wall is `steele-wall-argon`)
- CUDA pair style
- Public Config bags / factories
- Live LAMMPS in ctest
- SCF / PW / GPW / Torch

---
slug: fumi-tosi-nacl
status: approved
created: 2026-07-24
track: metalwalls-method-loop
chapter: E8
grilled: false
oracle: lammps
---

# fumi-tosi-nacl

## Summary

Ship multi-type **Born–Mayer–Huggins (Fumi–Tosi short-range)** pair Potential
on CPU as `cpu::PairBorn`, composable with `cpu::PairCoulCut` for molten-salt
Born+Coulomb (map E8). Correctness is dual-gated against an **offline LAMMPS**
`hybrid/overlay born + coul/cut` pin (metal units: eV, Å) on a multi-type NaCl
ion square; no live LAMMPS at ctest.

## Domain basis

- LAMMPS `pair_style born` (doc: `pair_born.html`):
  \[
  E(r)=A\,e^{(\sigma-r)/\rho}-\frac{C}{r^{6}}+\frac{D}{r^{8}},\quad r < r_c
  \]
  with optional energy shift (`pair_modify shift yes`) so \(E(r_c)=0\).
  Coefficient order: `A rho sigma C D` (same as LAMMPS `pair_coeff`).
- Coulomb: existing `cpu::PairCoulCut` = LAMMPS `coul/cut` (no energy shift).
  Twin uses `pair_style hybrid/overlay born … coul/cut …` with
  `pair_modify pair born shift yes`.
- Map E8 / queue: Fumi–Tosi molten salt; **dep** Buckingham pair ✓.
- Architecture: real-space pair under `cpu/potentials/`; **not** a new Driver
  category; Frame I/O only; parameters = ctor args / type table (no public
  `*Config`, no `make_*`/`build_*`).

## Design

### Public API (no `*Config` bags, no factories)

```cpp
// Single-type (ion dimer foothold)
cpu::PairBorn born(/*A=*/20.35, /*rho=*/0.317, /*sigma=*/2.755,
                   /*C=*/48.18, /*D=*/48.18, /*cutoff=*/10.0);

// Multi-type Fumi–Tosi table (Na=0, Cl=1 frame types)
cpu::PairBorn born({
    {0, 0, A_NaNa, rho, sig_NaNa, C_NaNa, D_NaNa},
    {1, 1, A_ClCl, rho, sig_ClCl, C_ClCl, D_ClCl},
    {0, 1, A_NaCl, rho, sig_NaCl, C_NaCl, D_NaCl},
}, /*cutoff=*/12.0);

// Born + Coulomb (Fumi–Tosi salt product path)
cpu::Evaluate eval;
eval.add_potential(std::move(born));
eval.add_potential(cpu::PairCoulCut{12.0});
eval.set_neighbor(cpu::md::NeighborList<>{12.0, 0.3});
eval.build();
eval.run(frame);
// frame.meta.get(kEnergy) [eV]; forces on kAtom_f*
```

- Class `ATV::cpu::PairBorn` inherits `Pair` / `PairPotential`.
- Multi-type: explicit type-pair table (`BornTypeCoeff`); no silent mixing of
  unset pairs (unset → zero interaction). Single-type ctor remains for the
  plain Born foothold.
- Born energy **shift** always on (matches LAMMPS `shift yes` / PairBuck).
- Coulomb is **not** folded into PairBorn — compose with `PairCoulCut` (or later
  `coul/long`). O(N²) Born path withdraws `NeighborConsumer`; Coulomb keeps NL.
- Units: eV, Å on the public surface.

### Physics / system (regression pin)

| Field | Value |
| ----- | ----- |
| Geometry | 4-ion planar square (2 Na + 2 Cl), \(d=2.82\) Å edge, box 40 Å cubic |
| Types / charges | Na type 0, \(q=+1\); Cl type 1, \(q=-1\) |
| Born params | see REFERENCE.md (Na–Na, Cl–Cl, Na–Cl) |
| Cutoffs | Born + Coulomb \(r_c=12\) Å; Born shift yes; Coulomb unshifted |
| Observable | total pair energy \(E\) [eV]; selected forces [eV/Å] |
| Oracle | **LAMMPS** `hybrid/overlay born + coul/cut`, units metal |
| Tolerance | \(\|E-E_{\mathrm{pin}}\| \le 1\times10^{-9}\) eV; \(\|F-F_{\mathrm{pin}}\| \le 1\times10^{-8}\) eV/Å |

### Offline LAMMPS oracle

- Binary: local `lmp` (30 Mar 2026 develop / patch cited in REFERENCE.md).
- Input twin: `tools/metalwalls_oracle/inputs/lammps/in.fumi_tosi_nacl4`.
- ctest **never** launches LAMMPS.

## Files to create or modify

- `.claude/specs/fumi-tosi-nacl/spec.md` (this)
- `.claude/specs/fumi-tosi-nacl/acceptance.md`
- `include/atomiverse/cpu/potentials/pair_born.h`
- `src/cpu/potentials/pair_born.cpp` (GLOB — no source-list edit)
- `tests/cpu/potentials/test_pair_born.cpp` + `tests/CMakeLists.txt`
- `regressions/forcefield/fumi-tosi-nacl/reg_fumi-tosi-nacl.cpp`
- `regressions/forcefield/fumi-tosi-nacl/REFERENCE.md`
- `regressions/CMakeLists.txt` — `atv_add_regression(forcefield fumi-tosi-nacl)`
- `tools/metalwalls_oracle/inputs/lammps/in.fumi_tosi_nacl4`

## Tasks

- [x] Spec + acceptance on `wf/fumi-tosi-nacl` before production code
- [x] Offline LAMMPS pin → REFERENCE.md (`oracle: lammps`, `oracle_version:`)
- [x] Implement multi-type `PairBorn` (shifted BMH; Frame I/O)
- [x] Unit test: analytic single-type \(E\), \(F\); multi-type type-pair lookup
- [x] E2E regression Born+Coul vs offline pin (public API only)
- [x] MD-only build; regression PASS
- [x] Perf gate vs same LAMMPS binary (or `perf_exempt` with reason)
- [x] Commit on `wf/fumi-tosi-nacl` only

## Non-goals

- `born/coul/long` + kspace productization (later; compose `CoulLong`/`KSpacePME`)
- CUDA pair style
- Public Config bags / factories
- Live LAMMPS in ctest
- SCF / PW / GPW / Torch
- Steele wall / wall Born (separate queue slugs)

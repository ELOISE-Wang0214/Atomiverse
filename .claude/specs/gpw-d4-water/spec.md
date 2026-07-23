---
slug: gpw-d4-water
status: approved
created: 2026-07-23
track: cp2k-method-loop
chapter: gpw
grilled: false
---

# gpw-d4-water

## Summary

Ship **periodic Grimme DFT-D4** (Becke–Johnson / rational damping, PBE
parameters, two-body only) as a geometry-only additive correction on
**Γ-point GPW**, gated on a **periodic water-dimer** cell against an offline
**dftd4 3.7.0** pin (CP2K 2026.1 on Arrhenius is built *without* the DFTD4
library — `qs_dispersion_d4` aborts — so the oracle is the same libdftd4
backend CP2K would call when linked). Public surface: `GpwDft::d4 = true`
(and `GPWEvaluator::d4`); `GpwResult::e_disp` [Ha] folded into `e_total`.

## Domain basis

- Caldeweyher et al., *J. Chem. Phys.* **150**, 154122 (2019); periodic
  extension Caldeweyher et al. (2019).
- Roadmap: `docs/roadmap/dft.md` §8 Dispersion (D4 ○ → energy pin).
- Queue: `docs/developer/cp2k-method-queue.md` `gpw-d4-water`.
- Prior pattern: `gpw-pbe-d3bj-ice` (D3(BJ) + `e_disp`); molecular
  `dft-d4-h2o-pbe` (`cpu::scf::D4`, non-periodic).
- Architecture: additive on `GpwDft` / `GPWEvaluator` members — **no** public
  `*Config`, no `make_*`/`build_*`. In-tree tables only (reuse
  `scf/internal/dftd4_tables.h`).

## Design

### Public API

```cpp
cpu::gpw::GpwDft gpw;
gpw.cell[0] = gpw.cell[1] = gpw.cell[2] = 12.0 * kAng2Bohr;
gpw.n_occ = 8;
gpw.xc = cpu::gpw::XcKind::PBE;
gpw.d4 = true;                 // mutually exclusive with d3 / d3bj
gpw.build(atoms);
const auto& res = gpw.run(atoms);
// res.e_disp [Ha]; res.e_total includes e_disp
```

Optional evaluator member: `GPWEvaluator::d4` mirrored into the solver.

### Physics (two-body only)

\[
E_{\mathrm{disp}}^{(2)} = -\sum_{A\le B}\sum_{\mathbf{T}}{}^{\prime}
  \left[
  \frac{s_6 C_6^{AB}}{R^6 + f^6} + \frac{s_8 C_8^{AB}}{R^8 + f^8}
  \right]
\]

with BJ \(f = a_1\sqrt{3 r_4^A r_4^B}+a_2\), lattice vectors \(\mathbf{T}\),
prime excludes \(R=0\). \(C_6\) from EEQ2019 charges + EN-weighted erf CNs
with **periodic** lattice sums (EEQ Coulomb via Ewald as in multicharge /
simple-dftd4). PBE-D4: \(s_6=1\), \(s_8=0.95948085\), \(a_1=0.38574991\),
\(a_2=4.80688534\) Bohr. **ATM / s9 off** (matches `--mbdscale 0` pin).

Default cutoffs (dftd4 defaults): CN 30 Bohr, disp2 60 Bohr.

### Mutual exclusion

`d3`, `d3bj`, and `d4` are pairwise exclusive — any two true throws at run.

## Geometry / oracle

Same water dimer as `gpw-pbe-d3bj-ice` (12 Å cubic, PERIODIC XYZ). Offline
pin from **dftd4 3.7.0** CLI on Arrhenius (`software/.../dftd4/3.7.0-gfbf-2025b`).
See `regressions/gpw/gpw-d4-water/REFERENCE.md`.

## Out of scope

- Three-body ATM / MBD scale
- D4 analytic forces / stress (energy pin only)
- Molecular-only path changes (existing `cpu::scf::D4` stays non-periodic
  Frame API; shared kernel may gain optional PBC for GPW)
- CP2K full QS e_tot parity (dispersion-only hard gate)
- CUDA / MPI industrial scaling

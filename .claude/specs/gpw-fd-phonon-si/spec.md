---
slug: gpw-fd-phonon-si
status: blocked
created: 2026-07-26
domain: response / gpw
queue: cp2k-method-loop
depends_on: [fd-phonon-supercell, gpw-forces-water]
unblocks: []
---

# Spec: `gpw-fd-phonon-si` — FD Γ phonons with GPW analytic forces

**Slug:** `gpw-fd-phonon-si`  
**Domain:** vibrational response (finite displacement) + GPW DFT forces  
**Queue:** CP2K method loop (wave 12)  
**Status:** approved for implement

## Summary

Close the **GPW force-consumer residual** on the existing crystal FD phonon
product: `cpu::vib::Phonon` already accepts a `GPW` kind on
`vib::Phase::force`, but the only product pin is Stillinger–Weber
(`reg_fd-phonon-supercell`). Deliver a textbook **Si₂ diamond-fragment**
Γ-point FD phonon spectrum driven by **GPW analytic forces**, offline-pinned
to CP2K `VIBRATIONAL_ANALYSIS` (Quickstep GPW).

## Domain basis

Central-difference Cartesian force constants on a fixed cell:

\[
C_{i\alpha,j\beta}
  = \frac{F_{j\beta}(R_0 + h\,\mathbf{e}_{i\alpha})
         - F_{j\beta}(R_0 - h\,\mathbf{e}_{i\alpha})}{2h},
\]

with \(h = 10^{-3}\) Å (`kHessianStepAngstrom`). Mass-weighted dynamical
matrix at Γ:

\[
D_{i\alpha,j\beta}
  = \frac{C_{i\alpha,j\beta}}{\sqrt{m_i m_j}},
\quad
\omega = \mathrm{eig}(D)\ \rightarrow\ \mathrm{cm}^{-1}.
\]

All \(3N\) modes kept (no molecular TR projection) —
`PhononDiagonalizer`. Acoustic branches near zero by translational
invariance of \(E\).

Public composition (no new factories / `*Config`):

```cpp
cpu::vib::Phonon phonon(/*na=*/1, /*nb=*/1, /*nc=*/1);
phonon.on(ATV::vib::Phase::force).add(
    cpu::gpw::GPWEvaluator{kinds, cutoff_ry, …, /*compute_forces=*/true});
phonon.build();
cpu::vib::PhononResult r = phonon.run(frame);
```

No `NeighborList` when the only force kind is GPW (no `NeighborConsumer`).

## System twin

| Knob | Value |
| ---- | ----- |
| Structure | Si₂ diamond fragment in cubic \(a = 5.43\) Å, periodic XYZ |
| Positions [Å] | Si (0,0,0), Si (\(a/4,a/4,a/4\)) |
| Basis | DZVP-MOLOPT-GTH |
| PP | GTH-PBE-q4 |
| XC | PADE (LSDA) |
| CUTOFF | 150 Ry |
| BC | **nonperiodic MT** (periodic Γ MOLOPT exceeds ≤60 s e2e budget) |
| SCF | `tol = 1e-7`, `max_iter = 80`, DIIS product path |
| FD step | \(h = 10^{-3}\) Å (CP2K `DX` = \(1.889726\times 10^{-3}\) Bohr) |
| `Phonon` | \(n_a=n_b=n_c=1\) (Γ of the given cell) |

## Acceptance

| ID | Criterion |
| -- | --------- |
| A1 | Mode count \(3N = 6\) |
| A2 | Acoustic residual: ≥3 modes with \(\|\omega\| < 80\) cm⁻¹ (Γ ASR) |
| A3 | Three non-acoustic modes within **50 cm⁻¹** of offline CP2K pin (sorted signed) |
| A4 | All frequencies finite; no NaN |
| A5 | Public API only: `Phonon` + `GPWEvaluator{compute_forces=true}` |
| A6 | Dual gate: \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) or `perf_exempt` with reason |
| A7 | Regression ≤ 60 s wall on epyc login (Release, 4 OMP threads) |

## Out of scope

- Finite-q / dispersion (already product on SW / PW DFPT tracks)
- Molecular `Hessian` TR projection path
- Non-orthorhombic GPW cells
- CUDA / k-mesh phonons

## Files

- Spec: `.claude/specs/gpw-fd-phonon-si/`
- Oracle input: `tools/cp2k_oracle/inputs/gpw-fd-phonon-si.inp`
- Regression: `regressions/gpw/gpw-fd-phonon-si/`
- Unit: extend `tests/cpu/vibration/test_phonon_supercell.cpp` — GPW wiring builds without NeighborList

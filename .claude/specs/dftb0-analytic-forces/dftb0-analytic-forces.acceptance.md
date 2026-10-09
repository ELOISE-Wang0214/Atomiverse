---
title: DFTB0 closed-form nuclear forces — acceptance
slug: dftb0-analytic-forces
oracle: dftb+
---

# Acceptance

Dual gate. Offline pin only. SE-only cmake (`ATV_BUILD_CLI=OFF` so MD
stays off). No unit-test binary: the force is the public `launch`, and a
full non-SCC solve does not belong in `tests/`.

## Regression (`reg_dftb0-analytic-forces`)

DFTB+ 24.1, mio-1-1, `SCC = No`, `PrintForces = Yes`. Transcripts under
`tools/tb_oracle/inputs/dftbplus/dftb0-analytic-forces/`.

| Quantity | Pin | Tol |
| -------- | --- | --- |
| H₂ \(R=0.74\) Å \(E\) | −0.6749509324 Ha, and the published −18.36635052 eV | \(1\times10^{-6}\) eV |
| H₂ \(R=0.74\) Å \(F_x\) | ∓0.001976354621 Ha/Bohr | \(1\times10^{-5}\) eV/Å |
| H₂ \(R=0.90\) Å \(E\) | −0.6630166170 Ha | \(1\times10^{-6}\) eV |
| H₂ \(R=0.90\) Å \(F_x\) | ±0.065975227977 Ha/Bohr | \(1\times10^{-5}\) eV/Å |
| off-eq H₂O \(E\) and \(F\) | `h2o/dftb.out.txt` | same |
| tilted CO \(E\) and \(F\) (pp + transverse) | `co/dftb.out.txt` | same |
| analytic vs energy FD, \(h=10^{-4}\) Å | identity, including a wrapped H₂ | \(1\times10^{-5}\) eV/Å |
| \(\lVert\sum_i F_i\rVert\) | 0 | \(1\times10^{-8}\) eV/Å |
| \(T_\mathrm{atv}/T_\mathrm{dftb+}\) | warm single point, 1 thread | \(\le 10\) |

`REFERENCE.md` states `oracle: dftb+` and `oracle_version: 24.1`.

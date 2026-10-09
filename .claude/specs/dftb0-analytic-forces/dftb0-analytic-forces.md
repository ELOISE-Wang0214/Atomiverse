---
title: DFTB0 closed-form nuclear forces
status: accepted
created: 2026-10-09
slug: dftb0-analytic-forces
oracle: dftb+
---

# DFTB0 closed-form nuclear forces

## Summary

Replace the 6N central-difference forces in `DFTB0::launch` with the
closed-form non-SCC gradient: Hellmann–Feynman + Pulay + pair repulsive.
No SCC shift, no Hubbard γ, no Δq. The H₂ mio-1-1 energy at 0.74 Å
(−18.36635052 eV) must not move.

Oracle is **DFTB+ 24.1** with `SCC = No` and mio-1-1. Not OpenMOPAC, not
xtb, and not an SCC run on a heteronuclear molecule.

## Domain basis

Elstner et al., *Phys. Rev. B* **58**, 7260 (1998). At fixed non-SCC
orbitals of \(H^0 C = S C \varepsilon\),

\[
F_A = -\frac{\partial}{\partial R_A}\Bigl(\sum_i 2\varepsilon_i + E_\mathrm{rep}\Bigr)
    = -\mathrm{Tr}\bigl(P\,\partial_{R_A} H^0\bigr)
      + \mathrm{Tr}\bigl(W\,\partial_{R_A} S\bigr)
      - \partial_{R_A} E_\mathrm{rep},
\]

with \(P_{\mu\nu} = 2\sum_i^\mathrm{occ} C_{\mu i} C_{\nu i}\) and
\(W_{\mu\nu} = 2\sum_i^\mathrm{occ} \varepsilon_i C_{\mu i} C_{\nu i}\).
Only the directed A→B Slater–Koster block is differentiated; the B→A
block contributes the same amount (factor 2). The radial channels and
their bond-vector derivative are the shared `MioSketchPair` already used
by DFTB2/DFTB3. Geometry is one minimum image, matching the energy.

## Design

- Same public type `cpu::semiempirical::DFTB0`. No `*Config`, no
  `make_*` / `build_*`. `with_forces` stays a public member.
- One diagonalization per launch. Forces reuse that density. `launch`
  still does not allocate a second method, and still accumulates energy
  and forces onto the Frame.
- Electronic pair contribution only inside the SK cutoff (where \(H^0\)
  and \(S\) were added). Repulsive derivative follows the spline cutoff
  inside `MioSkTable::d_repulsive_dr`, independent of the SK cutoff, so
  it matches the energy.
- Out of scope: cell stress, k-points, a lattice sum of images, spin,
  elements beyond H/C/N/O, an analytic derivative of the SK interpolant.

## Files

- `include/atomiverse/cpu/semiempirical/dftb0.h`
- `src/cpu/semiempirical/dftb0.cpp`
- `regressions/semiempirical/dftb0-analytic-forces/`
- `tools/tb_oracle/inputs/dftbplus/dftb0-analytic-forces/`
- `docs/roadmap/semiempirical.md` (DFTB0 note: analytic F, not 6N FD)

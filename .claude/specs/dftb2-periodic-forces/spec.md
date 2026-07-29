---
slug: dftb2-periodic-forces
status: approved
created: 2026-07-29
domain: semiempirical / dftb2
queue: tb-method-loop
oracle: dftb+|literature
method: dftb2
---

# Spec: `dftb2-periodic-forces` — lattice-image Γ-PBC nuclear forces

**Slug:** `dftb2-periodic-forces`  
**Domain:** SCC-DFTB2 (MioSketch) under Γ-point PBC  
**Queue:** TB method loop (`docs/developer/tb-method-queue.md`)  
**Oracle family:** **DFTB+** (formula / module boundary); absolute mio-1-1 numbers not required for this deepen slug  
**Unblocks:** lattice-consistent nuclear \(F\) on the shipped Γ energy path

## Summary

Γ-point **energy** already lattice-sums \(H^0\), \(S\), screened \(\gamma\), and the
repulsive pair potential (`build_h0_s` / `gamma_pair` / `repulsive_energy`).
Closed-form nuclear forces (`analytic_forces_restricted` /
`analytic_forces_unrestricted`) still walk **min-image only**. Under a dense
cell the energy and force operators disagree: \(E\) sees multiple lattice
images of a pair while \(F\) differentiates only one.

This slug makes **analytic \(F = -\nabla E\)** use the **same lattice-image
enumeration** as the energy, so a re-SCC central FD of the lattice-summed
energy matches the analytic force on an off-equilibrium 1-D carbon chain cell
(and a polar tight-box water check for \(\Delta q \ne 0\) \(\gamma\) terms).

## Domain basis

SCC-DFTB2 energy (Elstner et al., Phys. Rev. B **58**, 7260 (1998)):

\[
E = \mathrm{Tr}(P H^0) + \tfrac12\sum_{AB}\gamma_{AB}\Delta q_A\Delta q_B
  + E_{\mathrm{rep}},
\]

with Γ-point lattice sums

\[
H^0_{\mu\nu}
  = \sum_{\mathbf{T}} h_{\mu\nu}(\mathbf{R}_B+\mathbf{T}-\mathbf{R}_A),\quad
\gamma_{AB}
  = \sum_{\mathbf{T}} \gamma\bigl(|\mathbf{R}_B+\mathbf{T}-\mathbf{R}_A|\bigr)
\]

(onsite \(\mathbf{T}=0\) special-cased). Fixed-\(\Delta q\) nuclear gradient
(DFTB+ force module; Hourahine et al., JCP **152**, 124101 (2020)):

\[
F_A = -\frac{\partial E}{\partial R_A}
  = -\sum_{\mathbf{T}}
    \Big[
      \mathrm{Tr}\big(P\,\partial_{\mathbf{R}} H^0_{\mathbf{T}}
                     + \tfrac12(V_A+V_B) P\,\partial_{\mathbf{R}} S_{\mathbf{T}}
                     - W\,\partial_{\mathbf{R}} S_{\mathbf{T}}\big)
      + \Delta q_A\Delta q_B\,\partial_r\gamma\,
      + \partial_r E_{\mathrm{rep}}
    \Big]_{\mathbf{R}=\mathbf{R}_B+\mathbf{T}-\mathbf{R}_A}.
\]

Self-image terms (\(A=B\), \(\mathbf{T}\ne 0\)) depend only on lattice vectors
and contribute **zero** nuclear force (they still matter for stress / energy).

## In scope

1. Lattice-image walk in `analytic_forces_restricted` (MioSketch closed shell).
2. Same walk in `analytic_forces_unrestricted` (collinear open shell).
3. Stress pair-virial in the restricted path uses the **same** image bond
   vectors (keeps \(\sigma\) consistent with lattice \(F\) under Γ-PBC).
4. Unit test: off-eq C-chain analytic \(F_x\) vs re-SCC FD of lattice \(E\).
5. Regression `reg_dftb2-periodic-forces` + `REFERENCE.md` + offline oracle note.
6. Public header doc: analytic forces under Γ-PBC are lattice-consistent.

## Out of scope

- k-mesh beyond Γ (`dftb2-kpoints`)
- ThreeOb PBC (ThreeOb path is still molecular)
- Periodic D3 lattice sum (separate; `d3_periodic` if present elsewhere)
- Absolute DFTB+ / mio-1-1 force identity (sketch SK tables)
- DFTB3 / spin noncollinear

## Acceptance

See `acceptance.md`. Dual gate: correctness (analytic vs re-SCC FD under PBC)
+ `perf_exempt` when no DFTB+ binary / fair wall on host.

## API / architecture

- No new public types, no `*Config`, no `make_*` / `build_*`.
- Parameters remain members on `cpu::semiempirical::DFTB2`.
- Frame I/O only; Evaluate + Semiempirical Component.
- SE-only cmake (SEMIEMPIRICAL + EVALUATE + TESTS + REGRESSIONS ON; HF/MD/… OFF).

---
slug: dftb2-periodic-gamma
status: implemented
created: 2026-07-27
domain: cpu/semiempirical
queue: tb-method-loop
method: dftb2
oracle: dftb+|literature
deferred: false
depends_on: []
unblocks: [dftb2-analytic-stress]
grilled: false
---

# Spec: `dftb2-periodic-gamma` — Γ-point periodic SCC-DFTB2 energy

## Summary

Deliver **Γ-point periodic SCC-DFTB2** energy under PBC for
`cpu::semiempirical::DFTB2`: lattice-summed Slater–Koster H⁰/S, screened
SCC γ, and pair repulsive within the static SK cutoff. Frame in / energy
out via **Evaluate + Semiempirical** (SE-only cmake). Family oracle =
**DFTB+**; absolute energy pin is the ATV mio-style Hamiltonian twin (SK
not DFTB.org mio-1-1 bit-identical — same hygiene as `h2o_dftb2_mio`) with
literature method citations. k-point meshes and Ewald long-range SCC are
follow-ons.

## Domain basis

At the Γ-point, Bloch phase \(e^{i\mathbf{k}\cdot\mathbf{T}}=1\). Two-center
matrix elements become real-space lattice sums over cell translations
\(\mathbf{T}\):

\[
H^0_{\mu\nu}
  = \sum_{\mathbf{T}}
    h_{\mu\nu}(\mathbf{R}_\nu+\mathbf{T}-\mathbf{R}_\mu),\qquad
S_{\mu\nu}
  = \sum_{\mathbf{T}}
    s_{\mu\nu}(\mathbf{R}_\nu+\mathbf{T}-\mathbf{R}_\mu),
\]

truncated at the static pair cutoff \(R_\mathrm{cut}\) from
`dftb_mio_tables.h`. Onsite \(\mathbf{T}=0\), same atom: diagonal
\(e_s, e_p\) and \(S=1\). Self-images \(\mathbf{T}\ne 0\) contribute
hopping/overlap on the same atom block.

Screened SCC (Ohno–Klopman-like \(\gamma(R)\)) and exponential repulsive:

\[
\gamma_{AB}^\mathrm{eff}
  = \sum_{\mathbf{T}} \gamma(|\mathbf{R}_B+\mathbf{T}-\mathbf{R}_A|)
  \quad(\mathbf{T}=0,\,A=B \to U_A),
\]

\[
E_\mathrm{rep}
  = \tfrac12\sum_{A,B,\mathbf{T}}{}' A_{AB}\,e^{-B_{AB} R}
  \quad(\text{exclude }A=B,\mathbf{T}=0).
\]

**Not in this slug:** reciprocal-space Ewald for bare \(1/R\) Coulomb (γ is
already short-ranged/screened; real-space cutoff is the documented policy),
k-point sampling, analytic stress.

Molecular limit: large box with \(L \gg R_\mathrm{cut}\) recovers the
pre-existing min-image molecular energy (only \(\mathbf{T}=0\) images in
range).

References:

- Elstner et al., *Phys. Rev. B* **58**, 7260 (1998) — SCC-DFTB2.
- DFTB+ periodic machinery (architecture / oracle family only).
- Frauenheim / Aradi DFTB reviews — Γ-point and k-space TB under PBC.

## Public API (no `*Config`, no factories)

```cpp
cpu::semiempirical::DFTB2 dftb;
dftb.scc_tol = 1e-10;
dftb.with_forces = false;   // energy gate first; FD forces still ok
// Frame.box is the PBC cell; lattice sum is always on (molecular limit
// when L >> r_cut).
eval.on(evaluate::Phase::evaluate).add(std::move(dftb));
```

No new public bag type. Optional future member `k_mesh` stays out of scope.

## Acceptance

| Gate | Criterion |
| ---- | --------- |
| Unit | Molecular H₂ @ 20 Å: energy still finite (molecular limit) |
| Unit | Dense 1-D C dimer crystal (\(L_x \approx 2.8\) Å, \(r_\mathrm{C–C}=1.4\)): energy finite and **differs** from a large-box isolated dimer by more than 0.05 eV (proves lattice images enter H/S) |
| Unit | Unsupported Z still throws |
| Regression | `reg_dftb2-periodic-gamma`: 1-D C chain cell; finite E; hard pin ±0.05 eV vs REFERENCE twin |
| REFERENCE | `oracle: literature` (+ family DFTB+); geometry; pin; units eV; lattice-sum policy |
| Build | SE-only cmake; scoped targets only |
| Perf | `perf_exempt` if no offline DFTB+ binary; else \(T_\mathrm{atv}/T_\mathrm{oracle}\le 10\) |

## Out of scope

- k-point Monkhorst–Pack meshes
- Analytic stress / cell derivatives (`dftb2-analytic-stress`)
- True mio-1-1 / DFTB+ absolute energy identity (needs SK fidelity slug)
- Ewald / reciprocal SCC long-range
- MD/MINIMIZER/HF cmake ON; spin-DFTB; DFTB3

## Files

- `.claude/specs/dftb2-periodic-gamma/spec.md` (this)
- `src/cpu/semiempirical/dftb2.cpp` (lattice-sum H⁰/S/γ/rep)
- `include/.../dftb2.h` (doc: Γ-point PBC)
- `tests/cpu/semiempirical/test_dftb2.cpp` (PBC unit cases)
- `regressions/semiempirical/dftb2-periodic-gamma/{reg_*.cpp,REFERENCE.md}`
- `regressions/CMakeLists.txt` registration
- optional `tools/tb_oracle/inputs/dftbplus/dftb2-periodic-gamma/` twin notes

---
slug: ot-kpoints
status: in_progress
created: 2026-07-26
domain: gpw
queue: cp2k-method-loop
depends_on: [ot-minimizer-water, gpw-kpoints-si]
unblocks: []
---

# Spec: `ot-kpoints` — OT orbital minimizer with Monkhorst–Pack k-mesh

## Goal

Productize **Orbital Transformation (OT)** SCF on the **multi-k closed-shell
RKS** GPW path: `GpwDft::use_ot = true` together with `k_mesh` / explicit
`k_points` must run (no throw), converge to the same gapped fixed point as
multi-k DIIS, and expose the usual k payload (`k_frac`, `eps_k`, `C_occ_k`).

## Background

- Γ OT is product for RKS (`reg_ot-minimizer-water`) and UKS (`reg_ot-uks-o2`).
- Multi-k DIIS is product (`reg_gpw-kpoints-si` and family).
- Today `GpwDft::build` rejects `use_ot` on any k-path
  (`"k-point sampling does not support the OT minimizer"`).

## Scope

**In**

- Complex Bloch OT primitives (parameterize / tangent / gradient /
  Löwdin re-orthonormalize / FULL_SINGLE_INVERSE preconditioner) for
  Hermitian \(H_k,S_k\) and occupied \(C_k\) (\(n_b\times n_{\mathrm{occ}}\)).
- Multi-k RKS outer loop in `run_gpw_scf_kpoints`: density from occupied
  Bloch orbitals, rebuild \(H_k\), one PR+ CG OT step per k with shared
  Armijo line search on the k-summed total energy.
- Public surface: existing members only (`use_ot`, `k_mesh`, `k_points`, …).
- Unit tests (fixed-\(H_k\) band energy vs `diagonalize_k`; multi-k OT vs DIIS
  parity on a compact H₂ toy).
- Regression `regressions/gpw/ot-kpoints` with internal exactness gates +
  offline CP2K pin when available (else frozen first-green + `perf_exempt`).
- Oracle input twin under `tools/cp2k_oracle/inputs/ot-kpoints/`.

**Out**

- Multi-k **UKS** OT (still throw; separate productization).
- OT + metallic smearing (integer occupations only).
- OT + multi-k DFT+U (owned by `ot-dftu` / combo follow-on).
- CUDA OT; forces/stress under OT beyond existing DIIS-path reuse after
  convergence payload store.
- Live CP2K in ctest.

## API

```cpp
cpu::gpw::GpwDft m;
m.n_occ = /* … */;
m.use_ot = true;
m.k_mesh[0] = m.k_mesh[1] = m.k_mesh[2] = 2; // or explicit k_points
m.build(atoms);
const auto& r = m.run(atoms); // r.e_total [Ha]; r.n_kpts, r.eps_k, …
```

No public `*Config`, no `make_*` / `build_*` factories.

## Science

VandeVondele & Hutter, *J. Chem. Phys.* **118**, 4365 (2003); Kühne et al.,
*J. Chem. Phys.* **152**, 194103 (2020) §VII C. Complex generalization of the
in-tree real OT manifold: \(C^\dagger S C = I\),
\(C(x)=C_0\cos U + x\,U^{-1}\sin U\) with \(U=(x^\dagger S x)^{1/2}\).

## Acceptance

See `acceptance.md`.

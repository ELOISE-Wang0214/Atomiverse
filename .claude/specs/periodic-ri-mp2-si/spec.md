---
slug: periodic-ri-mp2-si
status: acceptance_met
created: 2026-07-25
domain: gpw / post-hf
queue: cp2k-method-loop
depends_on: [periodic-mp2-si]
unblocks: [periodic-rpa-si]
---

# Spec: `periodic-ri-mp2-si` — RI-GPW periodic MP2 productization

**Slug:** `periodic-ri-mp2-si`  
**Domain:** GPW wavefunction correlation (CPU)  
**Queue:** CP2K capability loop  
**Dep:** `periodic-mp2-si` ✓ (Γ canonical foothold)  
**Unblocks:** `periodic-rpa-si` (same MO-ERI / B intermediates)

## Summary

Productize **periodic closed-shell MP2** beyond the Γ-point canonical
foothold by shipping the **RI-GPW (plane-wave resolution-of-identity)**
assembly of the occupied–virtual MO repulsion integrals \((ia|jb)\)
(Del Ben, Hutter, VandeVondele, *J. Chem. Theory Comput.* **8**, 4177
(2012); CP2K `&WF_CORRELATION` / `METHOD MP2_GPW` + `&WFC_GPW`).

Foothold v1 (`GpwDft::mp2` with real-space Poisson-per-pair + full
\((ia|jb)\) tensor) stays available as the bit-identity reference path
(`mp2_ri = false`). Product default is `mp2_ri = true`: three-index
\(B_{ia}(\mathbf G)\) intermediates in reciprocal space, \((ia|jb)\) by
Hermitian contraction — same grid Coulomb as the Poisson path (identity
on the density grid), memory \(\mathcal O(n_{\mathrm{ov}} N_G)\) instead
of \(\mathcal O(n_{\mathrm{ov}}^2)\).

**Out of scope this wave:** multi-k MP2 (Bloch \(k_i+k_j=k_a+k_b\)
sum; complex C discarded on the k-path today); Gaussian cFIT/RI_AUX
tables (CP2K `BASIS_RI_MOLOPT` — separate data wave); UMP2; analytic
MP2 gradients/stress.

## Domain basis

Pair density on the SCF grid (Γ lattice-summed MOs):

\[
\rho_{ia}(\mathbf r)=\psi_i(\mathbf r)\psi_a(\mathbf r),\qquad
\int\rho_{ia}=S_{ia}=0.
\]

PW-RI (G=0 dropped, neutral pairs):

\[
\tilde\rho_{ia}(\mathbf G)=\mathrm{FFT}[\rho_{ia}],\qquad
B_{ia}(\mathbf G)=\sqrt{\frac{4\pi}{|\mathbf G|^2}}\,
\tilde\rho_{ia}(\mathbf G)\ \ (\mathbf G\neq 0),
\]

\[
(ia|jb)=\frac{w}{N}\sum_{\mathbf G}
B_{ia}^*(\mathbf G)\,B_{jb}(\mathbf G)
\quad\equiv\quad
w\sum_g V_{ia}(g)\,\rho_{jb}(g)
\]

with \(V_{ia}=\mathrm{Poisson}[\rho_{ia}]\) (existing foothold). Canonical
closed-shell doubles sum unchanged (`mp2_rhf_energy`).

## Public API

```cpp
cpu::gpw::GpwDft m;
// ... cell, cutoff_ry, n_occ, xc=PBE, level_shift, build, run ...
m.mp2_ri = true;           // product default; false → Poisson tensor path
m.mp2_n_virt = 0;          // 0 = all kept virtuals
const cpu::gpw::Mp2Result& r = m.mp2(atoms);
// r.ri == true; r.e_corr (< 0); r.e_total = r.e_scf + r.e_corr [Ha]
```

No public `*Config`, no `make_*` / `build_*`. Kernels private
(`internal/mp2_kernels.h`). Multi-k SCF (`n_kpts > 1` or non-Γ mesh)
throws `invalid_argument` (clear residual).

## Validation

| Tier | Gate |
|------|------|
| Unit | RI G-space \((ia|jb)\) ≡ Poisson tensor to 1e-10 on synthetic MOs; doubles sum identities; multi-k / UKS / OT guards |
| E2E | `regressions/gpw/periodic-ri-mp2-si/`: same sc-Si₂ deck as foothold; SCF converged; E_corr < 0; RI vs canonical |ΔE| ≤ 1e-6 Ha; offline CP2K pin band; OS+SS identity |
| Perf | Fair wall vs CP2K when pin host has binary; else `perf_exempt` |

## Acceptance

- [ ] Spec + offline CP2K input under `tools/cp2k_oracle/inputs/`
- [ ] `mp2_ri` product path + `Mp2Result::ri`
- [ ] Unit tests green (`test_cpu_gpw_mp2`)
- [ ] Regression green (`reg_periodic-ri-mp2-si`, ≤ 60 s)
- [ ] Dual gate: correctness + within_10x or perf_exempt
- [ ] Commit on `wf/periodic-ri-mp2-si`

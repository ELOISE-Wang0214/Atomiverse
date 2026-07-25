---
spec: periodic-mp2-kpoints
created: 2026-07-25
status: in_progress
method: gpw
oracle: CP2K
branch: wf/periodic-mp2-kpoints
deps:
  - periodic-ri-mp2-si
  - gpw-kpoints-si
---

# Spec — `periodic-mp2-kpoints`

## Goal

Lift the multi-k guard on `GpwDft::mp2` and ship **momentum-conserving
closed-shell GPW MP2** on top of a converged multi-k RKS reference
(Monkhorst–Pack / explicit list). Product path remains **RI-GPW** pair-density
Coulomb (Del Ben / CP2K `METHOD MP2_GPW`) with complex Bloch MOs; Γ real path
unchanged.

## Why

Queue residual after Γ product `periodic-ri-mp2-si` ✓ and k-SCF
`gpw-kpoints-si` ✓: `GpwDft::mp2` still throws when `n_kpts > 1` because the
k-path discards complex Bloch `C` (stores `eps_k` only). Roadmap §10 Periodic
MP2 notes multi-k residual.

## Public API (no Config / no factory)

```cpp
cpu::gpw::GpwDft m;
m.k_mesh[0] = 2; m.k_mesh[1] = 1; m.k_mesh[2] = 1;
m.k_time_reversal = false;   // equal weights for the multi-k sum foothold
m.mp2_ri = true;
m.mp2_n_virt = 4;            // regression budget
// ... build(); run(atoms);
const cpu::gpw::Mp2Result& r = m.mp2(atoms);  // no longer throws multi-k
// r.e_corr [Ha]; r.n_kpts mirrored via last_result().n_kpts
```

Members already public: `k_mesh`, `k_points`/`k_point_weights`, `mp2_ri`,
`mp2_n_virt`, `mp2_os_scale`, `mp2_ss_scale`. No new bag types.

## Science

1. **Store** complex Bloch MO coefficients `C_k` (nb × m_k) on the multi-k SCF
   result payload (alongside existing `eps_k` / weights).
2. **Complex MOs on grid** via existing `KpointBasisGrid`:
   \(\psi_p^{k}(g)=\sum_\mu C_{\mu p}^{k}\,\chi_\mu^{k}(g)\).
3. **Pair density** \(\rho_{ia}^{k_i k_a}(r)=\psi_i^{k_i*}(r)\,\psi_a^{k_a}(r)\)
   (complex). Coulomb via **linear** real Poisson on Re/Im:
   \(V=\mathrm{Poisson}(\Re\rho)+i\,\mathrm{Poisson}(\Im\rho)\) (G=0 dropped).
4. **Momentum conservation** \(k_b \equiv k_i + k_j - k_a\) (fractional, wrapped
   into the sampled set). Equal-weight meshes:
   \(E=\frac{1}{N_k}\sum_{k_i k_j k_a}\sum_{ijab}(\cdots)/\Delta\).
5. **Complex closed-shell doubles** (Γ-real reduces to existing OS/SS):
   - \(E_{\mathrm{OS}}=\sum |v_{iajb}|^2/\Delta\)
   - \(E_{\mathrm{SS}}=\sum \Re\!\big[v_{iajb}\,\overline{(v_{iajb}-v_{ibja})}\big]/\Delta\)

## System (regression deck)

| Item | Value |
|------|--------|
| Crystal | artificial sc-Si₂ (same as Γ RI-MP2) |
| Cell | cubic a = 5.43 Å, periodic |
| Basis / PP | DZVP-MOLOPT-GTH / GTH-PBE-q4 |
| XC | PBE, cutoff **30 Ry**, level_shift 0.3 Ha |
| k | **MP 2×1×1**, `k_time_reversal=false` → 2 k, w=½ each |
| MP2 | `mp2_ri=true`, `mp2_n_virt=4` (e2e ≤ 60 s) |

## Offline pin strategy

No CP2K binary on worker hosts (same status as prior waves). Binding gates:

1. **Γ recovery identity:** multi-k machinery with explicit
   `k_points = {(0,0,0)}`, w=1 matches real Γ `mp2()` to **1e-8 Ha** on
   `e_corr` (unit + regression).
2. **Hard multi-k gates:** SCF converged; `e_corr` finite and **strictly
   negative**; OS+SS ≡ e_corr; e_total ≡ e_scf+e_corr.
3. **Frozen first-green band** on the 2-k deck (setup integrity). CP2K input
   `tools/cp2k_oracle/inputs/periodic-mp2-kpoints.inp` ships for offline
   re-pin when a binary exists; do not invent absolute CP2K digits.

## Out of scope

- UMP2 / open-shell k-MP2
- Gaussian cFIT / `RI_AUX` tables
- Analytic MP2 gradients / stress
- TR-folded unequal-weight meshes (product later; foothold uses equal w)
- Periodic RPA / GW (separate slugs)
- Live CP2K in ctest

## Files (expected)

- `include/atomiverse/cpu/gpw/gpw_solver.h` — `GpwResult::C_k`; docs on `mp2`
- `include/atomiverse/cpu/linalg/matrix.h` — optional `CplxMatrix` alias
- `include/.../gpw/internal/mp2_kernels.h` + `src/.../mp2_kernels.cpp` — k kernels
- `src/cpu/gpw/mp2.cpp` — multi-k dispatch
- `src/cpu/gpw/gpw_solver.cpp` — store `C_k` in `store_k_payload`
- `tests/cpu/gpw/test_cpu_gpw_mp2.cpp` — complex OS/SS + Γ-through-k identity
- `regressions/gpw/periodic-mp2-kpoints/` + CMake row
- `tools/cp2k_oracle/inputs/periodic-mp2-kpoints.inp`

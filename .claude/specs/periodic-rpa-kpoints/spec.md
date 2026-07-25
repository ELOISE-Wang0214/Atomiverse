# Spec — periodic-rpa-kpoints

**Slug:** `periodic-rpa-kpoints`  
**Method:** Multi-k RI-GPW direct RPA (dRPA@PBE)  
**Oracle:** offline CP2K (`tools/cp2k_oracle/`; never live in ctest)  
**Parent branch:** `dev`  
**Worker branch:** `wf/periodic-rpa-kpoints`

## Goal

Close the multi-k residual of Periodic RPA: after Γ dRPA@KS product
(`reg_periodic-rpa-si`), `GpwDft::direct_rpa` must accept a multi-k SCF
reference and return a finite, strictly negative \(E_c\) via q-resolved
Furche plasmon formula on complex Bloch orbitals (RI-GPW pair densities).

## Science

- Closed-shell **dRPA** (\(\alpha_{\mathrm{hf}}=0\)) on top of multi-k
  **GPW PBE** KS orbitals.
- **Product multi-k path = k-diagonal (vertical \(q=0\)):** for each sampled
  \(k\), build Casida from same-\(k\) pair densities; BZ-average
  \(E_c=\sum_k w_k\tfrac12(\sum_s\Omega_{k,s}-\mathrm{Tr}\,A_k)\).
- \(A-B=\mathrm{diag}(\Delta\varepsilon)\), \(A+B=\mathrm{diag}(\Delta\varepsilon)+4V\).
- dRPA@KS total: \(E_{\mathrm{SCF}}+E_c\) (no EXX — periodic HF blocked).
- Product path: RI-GPW complex Poisson (same pair-density route as multi-k MP2).
- Recover real Γ plasmon when \(N_k=1\) and MOs are real (unit-tested).
- Full inter-\(k\) / finite-\(q\) dielectric coupling is a follow-on.

## Public API

```cpp
cpu::gpw::GpwDft m;
m.k_mesh[0] = 2; m.k_mesh[1] = 1; m.k_mesh[2] = 1;
m.k_time_reversal = false;
m.mp2_ri = true;
m.mp2_n_virt = 2;
// build(); run(atoms);
const cpu::gpw::RpaResult& r = m.direct_rpa(atoms); // multi-k, no throw
// r.n_kpts == 2; r.ri; r.e_corr < 0
```

No public `*Config`, no `make_*` / `build_*` factories. Parameters are
members on `GpwDft` (reuse `mp2_ri`, `mp2_n_virt`).

## System (regression)

Same deck as `periodic-mp2-kpoints`:

| Item | Value |
|------|--------|
| Crystal | diamond Si, bct 4-atom |
| a | 5.43 Å |
| Basis / PP | DZVP-MOLOPT-GTH / GTH-PBE-q4 |
| XC / cutoff | PBE / 80 Ry |
| Valence | n_occ = 8 |
| k | MP 2×1×1, equal weights |
| RPA | RI-GPW, n_virt = 2 |

## Acceptance

1. Multi-k SCF converges; `C_k` stored; `direct_rpa` does **not** throw for nk>1.
2. `RpaResult`: `n_kpts==2`, `ri==true`, `e_corr` finite and **strictly negative**,
   all \(\Omega>0\), `e_total ≡ e_scf + e_corr`.
3. Unit: \(N_k=1\) multi-k kernel ≡ real `direct_rpa_rhf_energy` on same MOs.
4. First-green absolute band on \(E_c\) (setup integrity; re-pin vs CP2K when
   binary available). Offline twin:
   `tools/cp2k_oracle/inputs/periodic-rpa-kpoints.inp`.
5. Wall ≤ 60 s single process. Perf: \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\le 10\)
   when fair twin exists; else **perf_exempt**.

## Out of scope

- URPA / open-shell RPA
- Frequency-integrated RI-RPA with Gaussian RI_AUX (CP2K default path differs)
- Periodic GW (separate slug)
- CUDA ports

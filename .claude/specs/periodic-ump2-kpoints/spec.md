---
spec: periodic-ump2-kpoints
created: 2026-07-26
status: done
method: gpw
oracle: CP2K
branch: wf/periodic-ump2-kpoints
deps:
  - periodic-ump2-si
  - periodic-mp2-kpoints
  - gpw-kpoints-uks
---

# Spec — `periodic-ump2-kpoints`

## Goal

Lift the multi-k guard on UKS `GpwDft::mp2` and ship **momentum-conserving
unrestricted GPW MP2 (UMP2)** on top of a converged multi-k UKS reference
(Monkhorst–Pack / explicit list). Product path remains **RI-GPW** pair-density
Coulomb (Del Ben / CP2K `METHOD MP2_GPW`) with complex Bloch α/β MOs; Γ UMP2
and multi-k RKS paths stay unchanged.

## Why

Queue residual after Γ UMP2 (`periodic-ump2-si` ✓) and multi-k closed-shell
MP2 (`periodic-mp2-kpoints` ✓): `GpwDft::mp2` still throws
`"multi-k UMP2 is a follow-on (Γ UKS only)"`, and multi-k UKS SCF stores only
occupied Bloch blocks (`C_occ_k` / `C_occ_k_beta`) — full virtuals needed for
canonical MP2 are discarded. Roadmap §10 notes multi-k UMP2 residual.

## Public API (no Config / no factory)

```cpp
cpu::gpw::GpwDft m;
m.n_alpha = 8; m.n_beta = 8;     // UKS closed-shell occupations
m.k_mesh[0] = 2; m.k_mesh[1] = 1; m.k_mesh[2] = 1;
m.k_time_reversal = false;       // equal weights
m.mp2_ri = true;
m.mp2_n_virt = 2;                // regression budget
// build(); run(atoms);
const cpu::gpw::Mp2Result& r = m.mp2(atoms);  // multi-k UMP2, no throw
// r.e_corr [Ha]; r.n_kpts == 2; r.ri == true
```

Members already public: `n_alpha`/`n_beta`, `k_mesh`, `mp2_ri`, `mp2_n_virt`,
`mp2_os_scale`, `mp2_ss_scale`. No new bag types.

## Science

1. **Store** full complex Bloch MO coefficients per spin after multi-k UKS
   SCF: `GpwResult::C_k` (α) and `C_k_beta` (β), alongside existing
   `eps_k` / `eps_k_beta` / occupied force payloads.
2. **Complex MOs on grid** via `mo_on_grid_k` for occ/virt of each spin.
3. **Same-spin αα / ββ** (PySCF UMP2):  
   \(E_{\mathrm{SS}}^\sigma=\frac1{4N_k}\sum_{k_i k_j k_a}\sum_{ijab}
   |v_{iajb}-v_{ibja}|^2/\Delta\) with \(k_b\equiv k_i+k_j-k_a\).
4. **Opposite-spin αβ**:  
   \(E_{\mathrm{OS}}=\frac1{N_k}\sum | (i_\alpha a_\alpha|j_\beta b_\beta) |^2
   /\Delta\) (no exchange), same momentum map.
5. **Identity:** when \(C^\alpha\equiv C^\beta\) and \(n_\alpha=n_\beta\),
   multi-k UMP2 ≡ multi-k RMP2 (dual-gate).

## System (regression deck)

| Item | Value |
|------|--------|
| Crystal | diamond Si, **bct 4-atom** (same as `periodic-mp2-kpoints`) |
| Cell | a_t = a/√2, c = a; a = **5.43 Å** |
| Basis / PP | DZVP-MOLOPT-GTH / GTH-PBE-q4 |
| XC | PBE UKS \(n_\alpha=n_\beta=8\), cutoff **80 Ry** |
| k | **MP 2×1×1**, `k_time_reversal=false` → 2 k, w = ½ |
| MP2 | RI-GPW, **n_virt = 2** per spin (e2e ≤ 60 s) |

## Offline pin strategy

No CP2K binary on worker hosts. Binding gates:

1. **Hard multi-k UMP2:** SCF converged with `C_k`/`C_k_beta`; `e_corr`
   finite and **strictly negative**; OS+SS ≡ e_corr; e_total ≡ e_scf+e_corr;
   `n_kpts` / `ri` flags.
2. **Dual-gate vs RKS multi-k:** same deck occupations ⇒ UMP2 E_corr within
   band of `periodic-mp2-kpoints` first-green (−0.0657637511 Ha).
3. **Kernel unit:** complex UMP2 OS/SS formulas at synthetic tensors; α=β
   recovers RMP2; multi-k Γ recovery of real UMP2 when nk=1 and real phases.
4. **CP2K input** `tools/cp2k_oracle/inputs/periodic-ump2-kpoints.inp` for
   offline re-pin; do not invent absolute CP2K digits.

## Out of scope

- Analytic UMP2 forces / stress
- Gaussian cFIT / `RI_AUX` tables
- TR-folded unequal-weight meshes (foothold uses equal w)
- UGW / URPA (separate slugs)
- Live CP2K in ctest

## Files (expected)

- `include/atomiverse/cpu/gpw/gpw_solver.h` — `GpwResult::C_k_beta`; docs
- `include/.../gpw/internal/mp2_kernels.h` + `src/.../mp2_kernels.cpp` —
  `mp2_uhf_energy_kpoints`
- `src/cpu/gpw/mp2.cpp` — multi-k UKS dispatch
- `src/cpu/gpw/gpw_solver.cpp` — store full `C_k` / `C_k_beta` in UKS payload
- `tests/cpu/gpw/test_cpu_gpw_mp2.cpp` — multi-k UMP2 kernel units
- `regressions/gpw/periodic-ump2-kpoints/` + CMake row
- `tools/cp2k_oracle/inputs/periodic-ump2-kpoints.inp`

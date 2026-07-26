# Spec — `periodic-urpa-kpoints`

## Goal

Close the multi-k unrestricted direct RPA residual: **URPA@UKS** with a
Monkhorst–Pack k-mesh on GPW. Product path = finite-\(q\) inter-\(k\)
spin-conserving dRPA over α→α ⊕ β→β (same dielectric as multi-k UGW
screening / multi-k closed-shell dRPA with spins explicit).

**Footholds (deps, already ✓):**
- `periodic-urpa-si` — Γ URPA@UKS (`direct_rpa_uhf`)
- `periodic-rpa-kpoints` / `periodic-rpa-finite-q` — multi-k RKS dRPA
  (`direct_rpa_rhf_energy_kpoints`)
- `periodic-ugw-kpoints` — multi-k UGW already builds the composite
  α⊕β finite-\(q\) Casida for \(W\); this slug reuses that dielectric for
  the Furche energy formula.

## Public API

```cpp
cpu::gpw::GpwDft m;
m.n_alpha = 8; m.n_beta = 8;          // UKS
m.k_mesh[0] = 2; m.k_mesh[1] = 1; m.k_mesh[2] = 1;
m.k_time_reversal = false;
m.mp2_ri = true;
m.mp2_n_virt = 2;
// build(); run(atoms);
const cpu::gpw::RpaResult& r = m.direct_rpa(atoms);  // multi-k URPA, no throw
// r.n_kpts == 2; r.ri; r.e_corr [Ha]; composite omega
```

**Remove** the throw:
`GpwDft::direct_rpa: multi-k URPA is a follow-on (Γ UKS only)`.

## Science

- Coulomb-only Casida, spins explicit (no closed-shell factor 2):
  \(A=\Delta\varepsilon\,\delta+V\), \(B=V\) on composite OV α then β
  per source \(k\).
- Finite-\(q\): transitions \(i_k\to a_{k+q}\), RI-GPW Poisson at transfer
  \(q\) (same as multi-k RKS dRPA / multi-k UGW).
- Furche:
  \(E_c=\sum_q w_q\,\tfrac12\bigl(\sum_s\Omega_{q,s}-\mathrm{Tr}\,A_q\bigr)\),
  \(w_q=n_{\mathrm{src}}(q)/N_k^2\).
- When \(C^\alpha\approx C^\beta\) and \(n_\alpha=n_\beta\), multi-k URPA
  \(E_c\) dual-gates multi-k RKS dRPA.
- Kernel unit: \(N_k=1\), real-valued complex MOs, α=β recovers real
  Γ URPA ≡ closed-shell plasmon (1e-10 Ha).

## Deck (e2e)

Same crystal as `periodic-rpa-kpoints` / `periodic-ugw-kpoints`:

| Item | Value |
|------|--------|
| Crystal | diamond Si, bct 4-atom |
| Cell | a_t = a/√2, c = a; a = 5.43 Å |
| Basis / PP | DZVP-MOLOPT-GTH / GTH-PBE-q4 |
| XC | PBE **UKS** \(n_\alpha=n_\beta=8\) |
| Cutoff | 80 Ry |
| k | MP 2×1×1, `k_time_reversal=false` → 2 k |
| RPA | RI-GPW, **n_virt = 2 per spin** |

## Acceptance

1. Spec + offline CP2K input twin under `tools/cp2k_oracle/inputs/`.
2. Unit: multi-k UHF \(N_k=1\) α=β dual-gates real RHF / real URPA E_c.
3. Regression `regressions/gpw/periodic-urpa-kpoints/`:
   - SCF converged; `C_k`/`C_k_beta` size == n_kpts;
   - `direct_rpa` no throw; `n_kpts==2`, `ri`;
   - \(E_c\) finite and **< 0**; all \(\Omega>0\);
   - e_total ≡ e_scf + e_corr;
   - n_occ = 16 (8+8), n_virt = 4 (2+2);
   - dual-gate band vs RKS multi-k first-green
     \(E_c=-0.0736273682\) Ha (tol 5e-3 Ha);
   - wall ≤ 60 s.
4. No public `*Config` / `make_*` / `build_*`.
5. `perf_exempt` if no fair CP2K plasmon twin (mirror prior RPA slugs).

## Out of scope

- Open-shell (\(n_\alpha\neq n_\beta\)) multi-k URPA e2e (kernel allows;
  regression uses closed-shell occupations).
- Frequency-integrated RI-RPA / EXX total-energy assembly.
- Molecular post-HF (pyscf-method-loop); classical MW; CUDA.

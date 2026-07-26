# Spec: ot-uks-smearing

**Slug:** `ot-uks-smearing`  
**Queue:** `docs/developer/cp2k-method-queue.md`  
**Branch:** `wf/ot-uks-smearing`  
**Oracle:** CP2K offline (`tools/cp2k_oracle`)  
**Domain:** GPW (collinear multi-k UKS)

## Intent

Productize **Orbital Transformation (OT)** SCF together with **metallic
occupation smearing** (Fermi–Dirac / Methfessel–Paxton) on the GPW multi-k
**collinear UKS** path. Today each feature is product alone:

| Path | Status |
| ---- | ------ |
| multi-k RKS OT + FD/MP | ✓ `ot-smearing` |
| multi-k UKS OT (integer occ) | ✓ `ot-uks-kpoints` |
| multi-k UKS DIIS + FD/MP | ✓ `gpw-uks-smearing` |
| multi-k UKS OT + FD/MP | **throws** at `build` / k-UKS OT |

Target fixed point: same Mermin free energy \(F = E_{\mathrm{KS}} - TS\) as
UKS DIIS + smear on the same system (identity gate). Absolute free-energy pin is
a frozen ATV first-green (CP2K absolute differs by grid/`EPS_DEFAULT` class —
same policy as `gpw-uks-smearing` / `ot-smearing`).

## Dependencies (done)

| Dep | Status |
| --- | ------ |
| `ot-smearing` | ✓ multi-k RKS OT + FD/MP |
| `ot-uks-kpoints` | ✓ multi-k UKS OT, integer occ |
| `gpw-uks-smearing` | ✓ multi-k UKS DIIS + FD/MP |

## In scope

1. **Public API:** keep members `use_ot`, `n_alpha`/`n_beta`, `smearing`,
   `degauss` [Ry], `extra_bands`, `mp_order` — no `*Config`, no factories.
2. **`build()`:** accept collinear UKS `use_ot && smearing != none` (multi-k
   path; smear always routes to k-path including effective Γ). Still reject
   multi-k UKS OT + DFT+U (other slug `ot-uks-dftu`).
3. **Multi-k UKS OT + smear algorithm:**
   - Dual spin OT manifolds expand to full independent AO spectrum
     (`n_active = m_min`) when smearing; integer path keeps `n_α` / `n_β`.
   - Shared \(\mu\) over \(\alpha+\beta\) spin-orbitals (`spin_degen=1`);
     occupations \(f_{\sigma kb}\in[0,1]\); density
     \(P_k^\sigma = C_k^\sigma \mathrm{diag}(f)\, (C_k^\sigma)^\dagger\).
   - Each outer step: rebuild \(H_k^\sigma\), subspace-canonicalize
     \(C^\dagger H C\), fill FD/MP, free energy \(F\) with entropy channel.
   - Projected gradient with non-uniform spinor occupations
     \(g = (I - S C C^\dagger)(2 H C f)\) (same helper as RKS OT+smear;
     \(f_i\in[0,1]\) for UKS).
   - Joint PR+ CG + Armijo on \(F\); soft gnorm + free-energy plateau
     (metals) matching RKS OT+smear policy.
4. **Unit tests:** `test_cpu_gpw_smearing` — UKS OT+smear build **accepted**;
   optional OT≡DIIS free-energy unit if cheap; keep reject for bad degauss.
5. **Regression:** `regressions/gpw/ot-uks-smearing/` — BCC Na twin of
   `gpw-uks-smearing` with `use_ot=true`; OT vs DIIS free-energy identity;
   electron count; frozen pin.
6. **Oracle input:** `tools/cp2k_oracle/inputs/ot-uks-smearing/` (DIIS+LSD+SMEAR
   twin; CP2K 2026.1 rejects OT+KPOINTS — document).

## Out of scope

- Multi-k UKS OT + DFT+U (`ot-uks-dftu`)
- OT + DFT+U + smear (`ot-dftu-smearing`)
- UKS OT + smear forces / stress / BOMD
- Live CP2K in ctest
- Molecular post-HF / MetalWalls / CUDA

## Acceptance

| ID | Gate | Tol |
| -- | ---- | ---: |
| A1 | `build` accepts UKS `use_ot` + FD | — |
| A2 | multi-k UKS OT+FD converges (BCC Na) | — |
| A3 | \(\lvert F_{\mathrm{OT}} - F_{\mathrm{DIIS}}\rvert\) | \(10^{-4}\) Ha |
| A4 | \(N_{\mathrm{occ}} = \sum_{\sigma k b} w f = N_e\) | \(10^{-5}\) |
| A5 | frozen free-energy pin (first green) | \(10^{-3}\) Ha |
| A6 | finite \(E_F\), entropy term | — |
| Perf | \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) or `perf_exempt` | — |

## Public surface (no Config)

```cpp
cpu::gpw::GpwDft m;
m.n_alpha = 10;
m.n_beta  = 8;   // ferro seed; free under smear
m.use_ot = true;
m.smearing = cpu::gpw::Smearing::fermi_dirac;
m.degauss = 0.02;      // Ry
m.extra_bands = 4;
m.k_mesh[0] = m.k_mesh[1] = m.k_mesh[2] = 2;
m.build(atoms);
const auto& r = m.run(atoms);  // r.e_total = F [Ha]
```

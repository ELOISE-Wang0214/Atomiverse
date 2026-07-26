# Spec: ot-smearing

**Slug:** `ot-smearing`  
**Queue:** `docs/developer/cp2k-method-queue.md`  
**Branch:** `wf/ot-smearing`  
**Oracle:** CP2K offline (`tools/cp2k_oracle`)  
**Domain:** GPW (closed-shell multi-k RKS)

## Intent

Productize **Orbital Transformation (OT)** SCF together with **metallic
occupation smearing** (Fermi–Dirac / Methfessel–Paxton) on the GPW multi-k RKS
path. Today each feature is product alone (`use_ot`, `smearing`);
`GpwDft::build` and the multi-k OT block **throw** when both are set.

Target fixed point: same Mermin free energy \(F = E_{\mathrm{KS}} - TS\) as
DIIS + smear on the same system (identity gate). Absolute free-energy pin is a
frozen ATV first-green (CP2K absolute differs by grid/`EPS_DEFAULT` class —
same policy as `gpw-metallic-smearing`).

## Dependencies (done)

| Dep | Status |
| --- | ------ |
| `ot-kpoints` | ✓ multi-k RKS OT, integer occ |
| `gpw-metallic-smearing` | ✓ multi-k RKS FD/MP |
| `ot-uks-o2` | ✓ Γ-UKS OT (no multi-k UKS OT) |

## In scope

1. **Public API:** keep members `use_ot`, `smearing`, `degauss` [Ry],
   `extra_bands`, `mp_order` — no `*Config`, no factories.
2. **`build()`:** accept RKS `use_ot && smearing != none`. Still reject
   multi-k UKS OT and multi-k OT+DFT+U (other slugs).
3. **Multi-k OT + smear algorithm (closed-shell):**
   - OT manifold size `n_active = min(m_min, n_occ + max(extra_bands, 1))`.
   - Each outer step: subspace-canonicalize \(C_k^\dagger H_k C_k\), fill
     FD/MP occupations at shared \(\mu\), density \(P_k = C_k f_k C_k^\dagger\),
     free energy \(F\) via existing `build_k_energy` entropy channel.
   - Projected gradient with **non-uniform** band occupations
     \(g = (I - S C C^\dagger)(2 H C f)\).
   - PR+ CG + Armijo on \(F\).
4. **Smear path routing:** when `smearing != none`, always use the multi-k SCF
   path (including effective Γ via `k_mesh = {1,1,1}`) so occupations are not a
   silent no-op on the real Γ DIIS/OT branch.
5. **Unit tests:** `test_cpu_gpw_smearing` — OT+smear build **accepted** for RKS;
   keep reject for bad degauss.
6. **Regression:** `regressions/gpw/ot-smearing/` — BCC Na twin of
   `gpw-metallic-smearing` with `use_ot=true`; OT vs DIIS free-energy identity;
   electron count; frozen pin.
7. **Oracle input:** `tools/cp2k_oracle/inputs/ot-smearing/` (OT + SMEAR twin).

## Out of scope

- Multi-k UKS OT + smear (`ot-uks-kpoints` residual)
- OT + DFT+U + smear / multi-k OT+DFT+U (`ot-dftu-kpoints`)
- DFPT metallic (`dfpt-metallic-si`)
- Live CP2K in ctest

## Acceptance

| ID | Gate | Tol |
| -- | ---- | ---: |
| A1 | `build` accepts RKS `use_ot` + FD/MP | — |
| A2 | multi-k OT+FD converges (BCC Na) | — |
| A3 | \(\lvert F_{\mathrm{OT}} - F_{\mathrm{DIIS}}\rvert\) | \(10^{-4}\) Ha |
| A4 | \(N_{\mathrm{occ}} = \sum w f = N_e\) | \(10^{-5}\) |
| A5 | frozen free-energy pin (first green) | \(10^{-3}\) Ha |
| A6 | multi-k UKS OT still throws | — |
| Perf | \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) or `perf_exempt` | — |

## Public surface (no Config)

```cpp
cpu::gpw::GpwDft m;
m.use_ot = true;
m.smearing = cpu::gpw::Smearing::fermi_dirac;
m.degauss = 0.02;      // Ry
m.extra_bands = 4;
m.k_mesh[0] = m.k_mesh[1] = m.k_mesh[2] = 2;
m.build(atoms);
const auto& r = m.run(atoms);  // r.e_total = F [Ha]
```

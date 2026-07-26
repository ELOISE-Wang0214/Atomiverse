# Spec — `periodic-urpa-si`

## Goal

Productize **Γ-point unrestricted direct RPA (URPA@UKS)** on top of a
converged **GPW UKS** reference, dual-gated against closed-shell dRPA when
\(n_\alpha = n_\beta\) and orbitals match.

## Science

Spin-conserving Coulomb-only Casida over α→α ⊕ β→β (spins explicit — **no**
closed-shell factor 2):

\[
A = \Delta\varepsilon\,\delta + V,\quad B = V,
\quad
E_c = \tfrac12\Bigl(\sum_s\Omega_s - \mathrm{Tr}\,A\Bigr).
\]

Same-spin \((ia|jb)\) via `mo_eri_ovov{,_ri}`; cross-spin via
`mo_eri_ovov_ab{,_ri}`. Kernel already ships as
`internal::direct_rpa_uhf` (UGW screening). Wire it into
`GpwDft::direct_rpa` for UKS Γ; multi-k URPA remains a follow-on.

When \(C^\alpha \approx C^\beta\) and \(n_\alpha=n_\beta\), Furche \(E_c\)
dual-gates closed-shell `direct_rpa_rhf` / `periodic-rpa-si` first-green.

## Public API

```cpp
cpu::gpw::GpwDft m;
m.n_alpha = 4; m.n_beta = 4;   // UKS
// ... members; build(); run(atoms);
const cpu::gpw::RpaResult& r = m.direct_rpa(atoms);
// r.e_corr, r.e_total [Ha]; r.omega_ha (composite α⊕β OV spectrum)
```

Parameters stay members (no `*Config`). Frame I/O unchanged (method-private
workspace). Public energy on Frame remains eV when evaluated via outer
drivers; internal Ha is the pin currency (cited).

## Scope

| In | Out |
|----|-----|
| Γ UKS dRPA energy (`direct_rpa`) | Multi-k URPA |
| RI-GPW product + Poisson foothold | Analytic URPA forces/stress |
| Unit: kernel dual-gate + unconverged/multi-k guards | Live CP2K in ctest |
| Regression: Si₂ UKS PBE deck | Molecular post-HF (PySCF loop) |

## Acceptance (see `acceptance.md`)

1. `GpwDft::direct_rpa` accepts UKS at Γ; multi-k UKS throws
   “multi-k URPA is a follow-on”.
2. Unit: `direct_rpa_uhf` Ec dual-gates `direct_rpa_rhf` on identical
   α=β tensors (1e-10 Ha); unconverged UKS throws.
3. Regression `reg_periodic-urpa-si`: SCF converged; \(E_c<0\); all
   \(\Omega>0\); \(E_{\mathrm{tot}}=E_{\mathrm{SCF}}+E_c\); dual-gate band
   vs `periodic-rpa-si` first-green \(|E_c-(-0.2658539776)|\le 2\times10^{-3}\)
   Ha; wall ≤ 60 s.
4. Offline oracle input under `tools/cp2k_oracle/inputs/periodic-urpa-si.inp`
   (never launched by ctest). `perf_exempt` if no fair wall twin
   (CP2K frequency RI-RPA ≠ plasmon).
5. No public `*Config` / `make_*` / factories.

## Files (expected)

- `src/cpu/gpw/rpa.cpp` — UKS Γ path
- `include/atomiverse/cpu/gpw/{gpw_solver,rpa}.h` — docs
- `tests/cpu/gpw/test_cpu_gpw_rpa.cpp` — guards + dual-gate surface
- `regressions/gpw/periodic-urpa-si/` — e2e + REFERENCE.md
- `tools/cp2k_oracle/inputs/periodic-urpa-si.inp`
- `regressions/CMakeLists.txt` — register

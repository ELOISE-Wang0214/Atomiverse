# Spec — `gpw-dftu-stress`

## Goal

Ship **analytic cell stress** for Γ-point UKS GPW with **Dudarev DFT+U**
(CP2K `PLUS_U` / Mulliken), closing the residual where `gpw_stress` threw
`"DFT+U not supported in this foothold"` after PLUS_U energy landed
(`gpw-dftu-nio`).

## Method

Hellmann–Feynman fixed-\(P\) product path (same skeleton as
`gpw-analytic-stress` / `gpw-uks-stress`):

1. One UKS+DFT+U SCF at the unstrained cell → \(P_\alpha,P_\beta\) and
   energy-weighted \(W_\sigma\) from KS eigenvalues **including** \(V_U\).
2. Central FD of the non-self-consistent Lagrangian
   \[
   E_L = \mathrm{Tr}((P_\alpha+P_\beta)h)
       + E_{\mathrm{es}} + E_{\mathrm{xc}}
       + E_U[n(P_\sigma,S)]
       - \mathrm{Tr}(W_\alpha S) - \mathrm{Tr}(W_\beta S)
   \]
   under orthorhombic diagonal strain at fixed fractional coordinates.
   \(E_U\) is the Dudarev energy on each spin (Mulliken \(n\) from fixed \(P\)
   and strained \(S\)).
3. \(\sigma_{aa} = -\Omega^{-1}\partial E_L/\partial\varepsilon_{aa}\)
   (CP2K `STRESS_TENSOR` sign; Ha / Bohr³).

Oracle twin: multi-SCF strain-FD (`gpw_stress_fd`) with the same `plus_u`
(copy already flows through `copy_gpw`).

## Scope

| In | Out |
|----|-----|
| Γ UKS + `plus_u` analytic stress | Multi-k DFT+U stress |
| Dual-gate vs multi-SCF FD | Live CP2K in ctest |
| NiO mol twin (SZV / PADE / U=6 eV) | Solid AFM rocksalt NiO |
| Frame `kStress` [eV/Å³] | Shear (off-diag remain 0) |

## Acceptance

1. `gpw_stress` / product `GPWEvaluator::compute_stress` **do not throw** when
   `plus_u` is non-empty (Γ path).
2. Unit: cheap H₂ + PLUS_U on H \(L=0\) — no-throw + dual-gate smoke
   max\(|\sigma^{\mathrm{an}}-\sigma^{\mathrm{fd}}|_{\mathrm{diag}} < 10^{-3}\) Ha/Bohr³
   (tight dual-gate is the NiO regression).
3. Regression `reg_gpw-dftu-stress`:
   - SCF energy soft-twin vs CP2K pin from `gpw-dftu-nio` (\(|\Delta E|\le 5\times10^{-2}\) Ha).
   - Dual-gate analytic vs FD diag max err \(\le 2\times10^{-4}\) Ha/Bohr³.
   - Off-diag product ~0; max\(|\sigma_{\mathrm{diag}}|\) floor.
   - \(e_{\mathrm{dftu}} > 10^{-4}\) Ha (U contribution present).
4. Perf: \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\le 10\) or `perf_exempt` with reason
   (CP2K NUMERICAL stress module unavailable on pin host → same pattern as
   `gpw-uks-stress`).

## Deps

- `gpw-dftu-nio` ✓ (PLUS_U energy)
- `gpw-analytic-stress` ✓ / `gpw-uks-stress` ✓ (fixed-\(P\) stress)

## Files

- `src/cpu/gpw/internal/gpw_stress.cpp` — \(E_U\) in `energy_fixed_P`
- `include/atomiverse/cpu/gpw/internal/dft_plus_u.h` — energy-only helper
- `tests/cpu/gpw/test_cpu_gpw_stress.cpp` — unit dual-gate + no-throw
- `regressions/gpw/gpw-dftu-stress/`
- `tools/cp2k_oracle/inputs/gpw-dftu-stress/`

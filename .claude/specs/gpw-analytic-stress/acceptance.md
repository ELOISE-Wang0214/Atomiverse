# Acceptance — `gpw-analytic-stress`

## Must pass

1. **Unit** (`test_cpu_gpw_stress`, GPW build):
   - `gpw_stress` (analytic) diagonals finite; off-diag exactly 0 (orthorhombic v1).
   - On cheap H₂ fixture: \(\max_a|\sigma_a^{\mathrm{an}}-\sigma_a^{\mathrm{FD}}| < 5\times 10^{-5}\) Ha/Bohr³
     (tight residual expected at `tol=1e-8`, \(h=3\times 10^{-3}\)).
   - `stress_axes` mask: frozen axes report exactly 0 for analytic and FD.
   - `GPWEvaluator` with `compute_stress=true`, `stress_use_fd=false` writes
     finite Frame `kStress` [eV/Å³] consistent with `last_stress_ha_bohr3()`.

2. **Regression** (`reg_gpw-analytic-stress`, ≤ 60 s class on login CPU):
   - Same H₂O / PBE / DZVP-MOLOPT-GTH / MT / 400 Ry twin as `gpw-stress-water`.
   - One Evaluate pass analytic stress + separate free-function FD twin
     (or `stress_use_fd` flip) on the same geometry.
   - \(\max_{\mathrm{diag}}|\sigma^{\mathrm{an}}-\sigma^{\mathrm{FD}}| \le 1\times 10^{-5}\) Ha/Bohr³.
   - \(|E - E_{\mathrm{cp2k}}| \le 1\times 10^{-3}\) Ha (same soft energy pin).
   - Off-diag analytic ≈ 0 within 1e-8 Ha/Bohr³.

3. **Architecture**
   - No public `*Config` / `make_*` / factories.
   - Stress helpers under `cpu/gpw/internal/`; public free functions on
     `gpw_solver.h` mirror PW (`pw_stress` / `pw_stress_fd`).
   - Default product path is analytic (`stress_use_fd=false`).

4. **Docs**
   - REFERENCE.md states included terms (fixed-\(P\) HF strain response) and
     deferred terms honestly.
   - Tutorial methods page + index row.

5. **Build**
   - `test_cpu_gpw_stress` and `reg_gpw-analytic-stress` green under
     `-DATV_BUILD_GPW=ON -DATV_BUILD_HF=ON -DATV_BUILD_PW=ON`.
   - Existing `reg_gpw-stress-water` still green (product default may be
     analytic; CP2K NUMERICAL pin residual within existing 5e-6 Ha/Bohr³
     **or** document a deliberate soft-up of that gate if residual moves
     by a known amount — prefer keep 5e-6 if residual allows).

## Out of acceptance

- CP2K `STRESS_TENSOR ANALYTICAL` absolute pin on MT water.
- Shear components; UKS / k-points / DFT+U / D4 stress.
- Term-by-term closed-form Pulay tables without fixed-\(P\) energy samples.

# Spec — second-gen-cpmd-water

**Slug:** `second-gen-cpmd-water`  
**Queue:** cp2k-method-loop  
**Domain:** dynamics (GPW AIMD)  
**Dep:** `cpmd-water` ✓  
**Branch:** `wf/second-gen-cpmd-water`

## Intent

Ship a **second-generation Car–Parrinello MD (SGCP)** foothold on closed-shell
Γ-point GPW water: Kolafa **ASPC** density predictor-corrector + a **partial
SCF corrector** (not full BOMD SCF, not fictitious-mass CPMD, not Niklasson
XL-BOMD). Pin short NVE energy-drift + structure envelope vs literature /
offline CP2K twin family.

## Literature / oracle

- T. D. Kühne, M. Krack, F. R. Mohamed, M. Parrinello, *Phys. Rev. Lett.*
  **98**, 066401 (2007) — second-generation CP (ASPC + few corrector steps).
- J. Kolafa, *J. Comput. Chem.* **25**, 335 (2004) — ASPC coefficients
  (already unit-tested in `test_accelerators` / `ASPCPredictor`).
- CP2K exercise `exercises:common:sgcp` — `EXTRAPOLATION ASPC` + limited
  SCF + optional Langevin (`NOISY_GAMMA` dissipation). Offline pin sketch under
  `tools/cp2k_oracle/inputs/second_gen_cpmd_water_nve.inp` when no binary.
- Twin geometry family: `gpw-bomd-water`, `xl-bomd-water`, `cpmd-water`.

## Product shape (no public `*Config`, no `make_*`)

Public members on `cpu::gpw::GPWEvaluator` (ctor-args / members only):

| Member | Default | Role |
| ------ | ------- | ---- |
| `second_gen_cp` | `false` | Engage SGCP path |
| `second_gen_cp_k` | `3` | ASPC order \(K\) (Kolafa) |
| `second_gen_cp_max_iter` | `5` | Corrector SCF budget (≪ full BOMD) |

Mutually exclusive with `cpmd` and `xl_bomd`. Closed-shell RKS diagonalization
only (`build()` throws on UKS / `use_ot`). Requires `compute_forces=true`.

**Algorithm (per force launch after history is warm):**

1. ASPC **predict** \(P_{\mathrm{pred}}\) from last \(K+1\) densities.
2. SCF with `guess_P = P_pred` and `max_iter = second_gen_cp_max_iter`
   (partial corrector — may not reach `tol`).
3. Forces at the **post-corrector** density (analytic path at that \(P\)).
4. Kolafa **correct** blend \(P_{\mathrm{corr}}\) and **record** for next step.

**Cold seed (insufficient ASPC history):** full `max_iter` SCF once, record
\(P\), BO forces — same pattern as CPMD/XL first step.

Reuse `cpu::scf::ASPCPredictor` (layout convert to/from `linalg::RealMatrix`).
Do **not** re-implement Lagrange coefficients.

## Non-goals (honest)

- Langevin / `NOISY_GAMMA` thermostat (optional production polish; NVE envelope OK for foothold).
- OT corrector (diagonalization corrector only this slug).
- UKS / k-points / electronic thermostat.
- Live CP2K in ctest.

## Acceptance

See `acceptance.md`. Dual gate: unit tests + `reg_second-gen-cpmd-water`.
`perf_exempt` if no fair CP2K twin wall on pin host.

## Files (expected)

- `include/atomiverse/cpu/gpw/gpw_evaluator.h` + `src/cpu/gpw/gpw_evaluator.cpp`
- `tests/cpu/gpw/test_gpw_second_gen_cpmd.cpp`
- `regressions/dynamics/second-gen-cpmd-water/{reg_*.cpp,REFERENCE.md}`
- `tools/cp2k_oracle/inputs/second_gen_cpmd_water_nve.inp` (offline sketch)
- `regressions/CMakeLists.txt` + test registration

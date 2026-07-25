# Acceptance — `cpmd-water`

## Unit (`tests/cpu/gpw/test_gpw_cpmd.cpp`)

| ID | Gate |
| -- | ---- |
| A1 | `Cpmd` constructs with \(\mu>0\); \(\mu\le 0\) throws |
| A2 | `seed` then `ready()`; `C()` shape preserved; `kinetic_ha()==0` at rest |
| A3 | `advance` with \(H=0\) keeps \(C\) fixed (up to Löwdin) and KE≈0 |
| A4 | `advance` before `seed` throws |
| A5 | `GPWEvaluator` `build()` rejects `cpmd` + UKS (`n_alpha/n_beta`) |
| A6 | `build()` rejects `cpmd` + `use_ot` |
| A7 | `build()` rejects `cpmd` + `xl_bomd` |
| A8 | `build()` rejects `cpmd` with `cpmd_mu≤0` or `cpmd_dt≤0` |
| A9 | Fixed-density path (`max_iter==0` + `guess_P`) returns finite \(E\), filled `H_ks`/`S`, `iters==0` (H₂ or minimal system if wall allows; otherwise matrix-level only via Cpmd) |

## Regression (`regressions/dynamics/cpmd-water`)

| ID | Gate |
| -- | ---- |
| A10 | Public API: `MD` + `VelocityVerlet` + `GPWEvaluator{…, compute_forces, cpmd=true}` |
| A11 | After warm-up: finite PE, PE bound \(<-100\) eV, finite KE_nuc ≥ 0 |
| A12 | Production chunks: finite coords; O–H distances \(< 1.8\) Å |
| A13 | CP conserved energy \(E_{\mathrm{cons}}=KE_{\mathrm{n}}+PE+T_{\mathrm{e}}\): relative drift \(\le 0.10\) (softer than BO 0.05 — Löwdin foothold + short traj) |
| A14 | After seed, force steps do **not** re-run multi-iter SCF (`last_result().iters==0` on production samples) — proves not BOMD rename |
| A15 | Wall ideally ≤60 s; if longer, document soft gate in REFERENCE.md (`perf_exempt`) |

## Docs

| ID | Gate |
| -- | ---- |
| A16 | Tutorial `docs/tutorial/methods/cpmd-water.md` + index row |
| A17 | Roadmap MD §AIMD Car–Parrinello → ✓ (foothold); queue row done on merge |
| A18 | Honest scope: UKS/OT/thermostat/RATTLE out; conserved quantity documented |

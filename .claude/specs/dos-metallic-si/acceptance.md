# Acceptance — dos-metallic-si

## Correctness

| # | Gate | Criterion |
| - | ---- | --------- |
| 1 | No throw | `Dos::compute` accepts `method.smearing = fermi_dirac` (no metallic foothold throw) |
| 2 | SCF | `dr.scf_converged == true` |
| 3 | Fermi level | `dr.e_fermi_ry` finite, \(\lvert E_F\rvert < 10^3\) Ry |
| 4 | Occupied DOS | \(\lvert n_\mathrm{elec\_to\_vbm} - n_e\rvert \le 1.0\) with \(n_e=8\) (integral to \(E_F\)) |
| 5 | Spectrum | DOS mesh yields sorted eigenvalues; `n_bands ≥` valence + empty |
| 6 | Insulator regression | Existing `reg_bands-dos-si` / fixed-occ unit tests still pass (no API break) |
| 7 | Unit | Pure-math test of integrate-to-\(\mu\) or metallic DOS weight path without full SCF e2e |

## Performance

| Item | Criterion |
| ---- | --------- |
| Wall | Regression ≤ 60 s |
| Dual gate | \(T_\mathrm{atv}/T_\mathrm{oracle} \le 10\) **or** `perf_exempt` with reason (no CP2K binary on host) |

## Architecture

- No public `*Config` / `make_*` / `build_*`
- Frame I/O only; parameters are members of `PW` / `Dos`
- Offline pin only in `REFERENCE.md`; never live CP2K in ctest
- Scoped build/test only (`atomiverse`, `test_cpu_pw_bands_dos`, `reg_dos-metallic-si`)

## Done when

All correctness gates green + dual gate recorded + commit on `wf/dos-metallic-si`.

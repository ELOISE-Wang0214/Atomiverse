# Acceptance — second-gen-cpmd-water

## Unit (`test_cpu_gpw_second_gen_cpmd`)

- [ ] `build()` rejects `second_gen_cp` + UKS occupations
- [ ] `build()` rejects `second_gen_cp` + `use_ot`
- [ ] `build()` rejects `second_gen_cp` + `cpmd`
- [ ] `build()` rejects `second_gen_cp` + `xl_bomd`
- [ ] `build()` rejects `second_gen_cp` without `compute_forces`
- [ ] `build()` rejects non-positive / non-finite `second_gen_cp_k` or `second_gen_cp_max_iter`
- [ ] ASPC history seeds after first full SCF; production steps use corrector
      budget (`last_result().iters ≤ second_gen_cp_max_iter` once history warm)
      — **no multi-step water e2e in unit tests** (matrix / guard only if SCF too
      heavy; optional tiny H2 plumbing OK if ms-scale)

## Regression (`reg_second-gen-cpmd-water`)

System: isolated H₂O, MT, PBE / DZVP-MOLOPT-GTH / GTH-PBE, 200 Ry, dt 0.5 fs,
ASPC \(K=3\), corrector `max_iter=5`, 1 warm-up + 4 production (twin family).

| Gate | Criterion |
| ---- | --------- |
| Finite energies | PE, KE finite; KE ≥ 0; PE < −100 eV after warm-up |
| Relative drift | \(\|E - E_0\| / \max(1,\|E_0\|) \le 0.10\) (soft vs BO 0.05 — partial SCF) |
| Structure | O–H < 1.8 Å |
| Method identity | After ASPC is warm, production `iters ≤ second_gen_cp_max_iter` (not full BOMD 60) |
| Wall | ≤ 60 s |

## Dual gate

- Correctness: literature envelope + method-identity gates in REFERENCE.md
- Perf: `perf_exempt` if no CP2K binary; else \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\)

## Definition of done

- [ ] Spec + acceptance committed on `wf/second-gen-cpmd-water`
- [ ] Implementation green under scoped build/test
- [ ] `acceptance_met = true`
- [ ] Single feature commit (message cites slug)

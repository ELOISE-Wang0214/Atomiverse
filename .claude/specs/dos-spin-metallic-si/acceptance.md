# Acceptance — `dos-spin-metallic-si`

## Correctness gates (regression `reg_dos-spin-metallic-si`)

| Gate | Criterion |
| ---- | --------- |
| SCF | `scf_converged` |
| Spin | `nspin == 2`, \(\lvert M\rvert \le 0.1\) (smear seed \(M=0\); free under FD) |
| \(E_F\) | finite, \(\lvert E_F\rvert < 10^3\) Ry |
| \(n_e(E_F)\) | \(\lvert n_{\mathrm{elec\_to\_vbm}} - 8\rvert \le 1.0\) (integral to \(E_F\)) |
| Channels | `dos_alpha` / `dos_beta` same length as total `dos` |
| Gap | `gap_ry > 0` on DOS mesh (gapped Si twin) |
| \(E_F\) placement | within \(\pm 4\,\mathrm{degauss}\) of \([E_\mathrm{vbm},E_\mathrm{cbm}]\) |
| Consistency | \(\lvert\int(D_\alpha+D_\beta)-\int D\rvert\) small on shared grid to \(E_F\) |
| API | `Dos::compute` does **not** throw for spin+smear; `Bands` still rejects |

## Twin settings

Diamond Si, \(a=5.43\) Å, ecut 15 Ry, SCF k \(2\times2\times2\), LSDA, ONCV UPF,
`spin_polarized=true`, `multiplicity=1`, FD `degauss=0.02` Ry, `extra_bands=4`,
DOS k=4, `n_bands=8`, \(\sigma=0.05\) Ry, 200 bins.

## Performance

| Item | Criterion |
| ---- | --------- |
| Fair wall | \(T_{\mathrm{ATV}}/T_{\mathrm{CP2K}}\le 10\) when binary available |
| Else | `perf_exempt=true` with reason (no module / no fair twin) |

## Offline oracle

`tools/cp2k_oracle/inputs/dos-spin-metallic-si.note` — LSD + SMEAR FERMI_DIRAC
+ PDOS on diamond Si. Never live CP2K in ctest.

## Dual gate

1. Unit: `ctest -R test_cpu_pw_bands_dos` (scoped).
2. Regression: `--target run_reg_dos-spin-metallic-si` (or build `reg_dos-spin-metallic-si` + run).

`acceptance_met=true` only when both green and commit lands on `wf/<slug>`.

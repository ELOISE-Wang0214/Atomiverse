# Acceptance — mpi-pw-scaling-4rank

Twin: Si diamond 2-atom fcc primitive cell, a = 5.43 Å, ONCV LDA UPF
(`tests/data/pw/Si_lda_nc.upf`), ecut 15 Ry, `k_mesh = 2` (8 k-points),
tol 1e-6, `max_iter = 50`. External pin: PseudoDojo `Si.djrepo` deltafactor
E_min = −16.916316 Ry/cell (ABINIT-class, same PP; van Setten 2018,
Hamann 2013) — the identical cited band as `regressions/pw/pw-lda-si`.

| # | Criterion | Observable | Tolerance / bound | Verified by |
| - | --------- | ---------- | ----------------- | ----------- |
| A1 | Serial distribution observable | `PWResult::n_kpts`, `n_kpts_local` after a 1-rank run | `n_kpts == 8`, `n_kpts_local == n_kpts` | `test_mpi_pw_kpoint_n1` (new TEST) |
| A2 | Round-robin partition at n = 2, 4 | Allreduce over ranks of `n_kpts_local` | SUM == `n_kpts`; MAX − MIN ≤ 1 | `test_mpi_pw_kpoint_n2` / `_n4` (new TEST) |
| A3 | 4-rank correctness vs external pin | child `energy_ry` at np = 1 and np = 4, both converged | each in [−16.916316, −16.916316 + 1.50] Ry (variational incomplete-basis band) | `run_reg_mpi-pw-scaling-4rank` |
| A4 | Rank-count invariance (parity) | E(np=1) − E(np=4) | abs ≤ 1e-8 Ry (deterministic Davidson seeding) | `run_reg_mpi-pw-scaling-4rank` |
| A5 | Real 4-way distribution | child-reduced `nk_local_min` / `nk_local_max` at np = 4 | min == max == 2 (8 k / 4 ranks) | `run_reg_mpi-pw-scaling-4rank` |
| A6 | Scaling wall ratio | T_run(np=1) / T_run(np=4), equal OMP threads per rank | ≥ 1.2 (weak industrial floor; measured value recorded in REFERENCE.md) | `run_reg_mpi-pw-scaling-4rank` |
| A7 | Budget | total regression wall (both launches) | ≤ 60 s | `run_reg_mpi-pw-scaling-4rank` |
| A8 | Perf vs oracle | T_atv/T_cp2k | ≤ 10, else `perf_exempt` with reason | REFERENCE.md — **perf_exempt: no CP2K binary/module on this host (checked 2026-07-25)** |

Gate commands (scoped):

```
cmake --build build -j16 --target atomiverse test_mpi_pw_kpoint reg_mpi-pw-scaling-4rank
ctest --test-dir build -R '^test_mpi_pw_kpoint' --output-on-failure
cmake --build build --target run_reg_mpi-pw-scaling-4rank
```

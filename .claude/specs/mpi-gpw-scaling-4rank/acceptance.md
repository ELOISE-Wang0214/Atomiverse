# Acceptance: `mpi-gpw-scaling-4rank`

## Correctness

| ID | Gate | Tolerance |
| -- | ---- | --------- |
| A1 | Multi-rank (n=4) SCF converges; finite `e_total` | `converged == true`; `isfinite` |
| A2 | Intra-job rank agreement: every rank holds the same `e_total` after Allreduce | \(\|E_{\min}-E_{\max}\| \le 10^{-12}\) Ha |
| A3 | Rank-count invariance: serial (n=1 library path or `mpiexec -n 1`) vs n=4 | \(\|E_1 - E_4\| \le 10^{-8}\) Ha (SCF tol band) |
| A4 | `n_kpts` matches mesh (MP 2×2×2 TR → 4); `n_kpts_local` partitions under RR | sum of local counts over ranks = `n_kpts`; n=1 ⇒ local = total |
| A5 | Serial path with MPI **not** initialized (or nprocs=1) matches pre-MPI energy of the same twin within SCF tol | \(\le 10^{-8}\) Ha vs frozen first-green or self-consistent n=1 |

## Performance

| ID | Gate | Value |
| -- | ---- | ----- |
| P1 | Fair wall vs CP2K multi-rank | **perf_exempt** — no domain decomp; k-RR only; no CP2K twin required this wave |
| P2 | Record \(T_{n=1}/T_{n=4}\) wall ratio in REFERENCE / stdout | informational only; no hard floor |

## Process

- [ ] Spec + acceptance present under `.claude/specs/mpi-gpw-scaling-4rank/`
- [ ] Library k-RR + Allreduce on GPW multi-k path
- [ ] Regression `reg_mpi-gpw-scaling-4rank` green under MPI build
- [ ] Optional unit `test_mpi_gpw_kpoint` n=1/2/4 green
- [ ] Tutorial page `docs/tutorial/methods/mpi-gpw-scaling-4rank.md`
- [ ] Commit on `wf/mpi-gpw-scaling-4rank` only

# mpi-pw-scaling-4rank — MPI PW k-mesh round-robin scaling gate (4 ranks)

- **Queue:** `docs/developer/cp2k-method-queue.md` (`mpi-pw-scaling-4rank`,
  wave 3 `[~]`, staffed claude wave 2 2026-07-25)
- **Domain:** pw (infra gate)
- **Branch:** `wf/mpi-pw-scaling-4rank` (worktree off `dev`)
- **Status:** done

## Problem

The plane-wave SCF (`cpu::pw::PW`) already distributes k-points round-robin
over `MPI_COMM_WORLD` (density + band-energy Allreduce each iteration), and
`tests/mpi/test_pw_kpoint.cpp` asserts cross-rank *agreement* at n = 1/2/4.
What does **not** exist is the multi-rank **industrial gate**:

1. No end-to-end regression runs the Si NC-UPF twin at 4 ranks and pins the
   converged energy against the external oracle band (rank-count invariance
   as a shipped, repeatable gate).
2. Nothing observes the k-point distribution — a silently serial run (layout
   never active, all k on every rank) would pass every existing test.
3. No wall-ratio scaling evidence (queue metric: "wall ratio") is recorded
   anywhere.

## Deliverable

1. **Observable distribution (library, minimal):** `PWResult` gains two
   plain members set by `PW::run()`:
   - `n_kpts` — total Monkhorst–Pack k-points of the run (`k_mesh`³ reduced
     set as generated, i.e. `prob.kg.kpoints.size()`),
   - `n_kpts_local` — k-points owned by *this* rank under the round-robin
     layout (serial / MPI-inactive: `n_kpts_local == n_kpts`).
   No new types, no config bag — two ints on the existing result struct.
2. **MPI unit assertion** (existing tier `tests/mpi/test_pw_kpoint.cpp`,
   already launched at n = 1/2/4): one new TEST asserting the distribution
   partition — `Allreduce SUM(n_kpts_local) == n_kpts` and
   `max − min ≤ 1` across ranks (single SCF pass, `max_iter = 1`, so the
   assertion stays cheap).
3. **Regression** `regressions/pw/mpi-pw-scaling-4rank/`
   (`reg_mpi-pw-scaling-4rank.cpp`, registered via `atv_add_regression`,
   guarded by `ATV_USE_MPI AND MPIEXEC_EXECUTABLE`): flat `main()`, public
   headers only. The binary self-launches its MPI child runs through the
   configured `mpiexec` (compile-time pin `ATV_MPIEXEC`, same pattern as
   `ATV_PW_UPF_PATH`); the child mode (`--child <out>`) runs the
   `pw-lda-si` Si diamond twin (ONCV LDA UPF, ecut 15 Ry, `k_mesh = 2` →
   8 k-points, tol 1e-6) under `MPI_Init`, times `run()` with `MPI_Wtime`,
   reduces the per-rank distribution, and rank 0 reports
   `energy_ry iters converged n_kpts nk_local_min nk_local_max wall_s`.
   Parent compares np = 1 vs np = 4 (equal `OMP_NUM_THREADS` per rank).
4. **Pin file** `REFERENCE.md`: external correctness pin (PseudoDojo
   ABINIT-class deltafactor E_min for this exact ONCV Si PP — the same
   cited oracle band as `regressions/pw/pw-lda-si`), CP2K comparability
   note per `tools/cp2k_oracle/README.md`, measured wall ratio, perf
   exemption (no CP2K binary/module on this host — verified 2026-07-25).
5. **Docs:** `docs/tutorial/methods/mpi-pw-scaling-4rank.md` product page
   (loop convention).

## Non-goals

- No new parallel decomposition (G-vector / band / domain) — k-point RR is
  the shipped capability; this slug is the gate.
- No CUDA, no GPW MPI (sibling slug `mpi-gpw-scaling-4rank`).
- No live oracle at build/ctest/regression time (offline pins only).

## Architecture conformance

- Concrete method `cpu::pw::PW`, members as parameters, Frame in / result
  out; new fields live on the existing public result struct (not a bag).
- Regression is an out-of-library consumer of public headers (+ `<mpi.h>`,
  which propagates from the public `MPI::MPI_CXX` link).
- No factories, no new public types, no free functions.

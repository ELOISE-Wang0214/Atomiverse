# Spec: `mpi-gpw-scaling-4rank`

**Slug:** `mpi-gpw-scaling-4rank`  
**Domain:** infra / gpw  
**Status:** in progress  
**Queue:** `docs/developer/cp2k-method-queue.md` (industrial gate; `perf_exempt` OK until domain decomp lands)

## Goal

Add **k-point MPI round-robin** to the GPW multi-k RKS SCF path, mirroring the
plane-wave pattern (`cpu::pw` / `internal/scf_loop.cpp` + `pw_solver.cpp`):

- Partition Monkhorst–Pack / explicit k-list across `MPI_COMM_WORLD` ranks
  (`ik = rank; ik < n_kpts; ik += nprocs`).
- Each rank collocates only its local k contribution to \(\rho(\mathbf{r})\),
  builds \(H_k\), diagonalizes, and accumulates band / kinetic / nonlocal
  traces for those k only.
- `MPI_Allreduce` the real-space density accumulation, energy scalars, and
  DIIS \(B\)-matrix partial sums so every rank obtains **identical**
  `e_total` (and the same SCF trajectory as serial, within SCF noise only
  if any; target is rank-count invariance).
- Serial / single-rank path remains **bit-identical** when `nprocs == 1` or
  MPI is not initialized (no Allreduce branch; loop stride 1).

This is a **foothold industrial gate**, not full GPW domain decomposition
(real-space grid / AO distribution). Setup of \(S_k, T_k, V_{\mathrm{nl}}^k\)
and the pre-summed k-basis grid may stay replicated on every rank for this
wave — honest partial scaling.

## Non-goals

- Real-space domain decomposition of the density grid or FFT.
- Multi-rank Γ-only GPW (no k mesh → nothing to round-robin; stays serial work).
- k-point UKS, k-point analytic forces, MPI for OT / MP2 / XL-BOMD paths.
- Fair wall twin vs CP2K (`perf_exempt` until decomp + offline twin exist).

## Public API

- No new `*Config` bags / factories.
- Optional payload on `GpwResult`:
  - `n_kpts` — total sampled k after TR fold / explicit list.
  - `n_kpts_local` — k owned by this rank under the round-robin layout
    (equals `n_kpts` when MPI inactive / `nprocs == 1`).
- Existing members (`k_mesh`, `k_points`, …) unchanged; MPI engages when
  built with `ATV_HAS_MPI`, MPI is initialized, and `nprocs > 1`.

## Library changes

| File | Change |
| ---- | ------ |
| `src/cpu/gpw/gpw_solver.cpp` | `run_gpw_scf_kpoints`: MPI layout, k RR loops, Allreduce \(\rho\), \(e_{\mathrm{kin}}\), \(e_{\mathrm{nl}}\), `err_norm` (MAX); fill `n_kpts` / `n_kpts_local`; rebuild full \(H_k\) from last \(V_{\mathrm{eff}}\) for band payload on all ranks |
| `src/cpu/gpw/internal/kpoints.cpp` | `KDiis::step`: Allreduce the Pulay \(B_{ij}\) block when MPI multi-rank so DIIS coefficients match serial (local \(e^k\) only non-zero per rank) |
| `include/atomiverse/cpu/gpw/gpw_solver.h` | Document + expose `n_kpts`, `n_kpts_local` on `GpwResult` |

## Regression

`regressions/gpw/mpi-gpw-scaling-4rank/`:

- Public headers only; Si bct multi-k twin (lightweight cousin of
  `gpw-kpoints-si`: MP 2×2×2 TR, modest cutoff / tol for ≤ 60 s under
  self-launch).
- **Parent** process (no MPI init): runs serial energy + wall; self-launches
  `mpiexec -n 4` of the same binary in worker mode; pins
  \(|E_{n=1} - E_{n=4}|\) rank-invariance; records wall ratio
  \(T_1 / T_4\) (informational; **perf_exempt**).
- **Worker** mode (`--mpi-worker`): `MPI_Init`, SCF, Allreduce min/max energy
  across ranks (intra-job agreement), print `ENERGY=` / `WALL_S=` for parent.

## Unit (optional)

`tests/mpi/test_gpw_kpoint.cpp` — CTest `mpiexec` n=1,2,4 (same pattern as
`test_pw_kpoint`): all ranks agree on `e_total`; finite + converged.

## Build

```bash
module load buildenv-gcc/2026.03-mpich
cmake -B build -DCMAKE_BUILD_TYPE=Release \
  -DATV_BUILD_TESTS=ON -DATV_BUILD_HF=ON -DATV_BUILD_PW=ON \
  -DATV_BUILD_GPW=ON -DATV_BUILD_CUDA=OFF -DATV_USE_MPI=ON \
  -DATV_BUILD_REGRESSIONS=ON -DBLA_VENDOR=OpenBLAS
```

## Acceptance

See `acceptance.md`.

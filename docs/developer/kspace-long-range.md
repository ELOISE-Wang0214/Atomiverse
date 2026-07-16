# KSpace long-range electrostatics (CPU)

This note is the durable contract for `CoulLong` + `KSpaceEwald` /
`KSpacePME` composition, cell refresh, stress, and MPI scope. It closes
cpu-breadth-03 **ac-009** documentation.

## Pairing (fail closed)

Real-space and reciprocal-space halves are ordinary components that implement
`LongRangeParticipant` and declare:

| Field | Meaning |
|---|---|
| `role` | `RealSpace` (`CoulLong`) or `ReciprocalSpace` (`KSpaceEwald` / `KSpacePME`) |
| `alpha` | Ewald splitting α [1/Å] — **must match** between the pair |
| cutoff / accuracy | cold metadata only |

`CPUDriver` validation at `build()` accepts **exactly one** real-space plus
**exactly one** reciprocal partner with matching `alpha`. Missing, duplicate,
or mismatched partners throw before launch with component names and α values.
There is no short-range or energy-only fallback.

NPT / any component with `requires_system_stress()` further needs the
reciprocal style to **write** complete reciprocal stress into system meta.
PME does so after the Essmann influence-function virial implementation.

## Stress (system virial)

Stress is the **system** 6-component virial in `frame.meta` under key
`"stress"` (`f.meta.get(kStress)`), Voigt order xx,yy,zz,xy,xz,yz, units eV
(LAMMPS metal / force·distance). It is **not** a per-atom column.

- Real-space: `CoulLong` accumulates pair virial `r·F`.
- Reciprocal: Ewald and PME write the influence-function geometric factor
  `V_ab = ug · (δ_ab + vt · k_a k_b)` (m-space form for PME).
- Kinetic: `Virial` observer → `kStressKinetic`.
- Barostat pressure: `P = tr(σ_virial + σ_kinetic) / (3V)`, positive compressive.

## Cell refresh vs capacity rebuild

`KSpace` owns a lifecycle ledger (`lifecycle_stats()`):

| Event | When | Cost |
|---|---|---|
| `builds` | first `build()` | allocate atom/grid scratch, create FFT plans (PME) |
| `cell_refreshes` | Box identity or revision changes, **same** FFT/grid shape | update reciprocal basis / influence; overwrite persistent spline scratch (`fill`, not re-`assign`) |
| `capacity_rebuilds` | atom capacity growth **or** FFT mesh dimensions change | replan FFT, reallocate grid |

Caches key on `(Box*, revision)`. Replica frames with equal revisions but
distinct Box identities still refresh.

Instrumentation gate: `PmeLifecycle.StableGridRefreshIsReplanFreeAndMemoryStable`
and related unit tests.

## MPI ceiling (replicated grid)

Current PME is **replicated-grid**: each rank spreads owned charges, the grid
is `MPI_Allreduce`-summed, every rank runs the full FFT, forces are
interpolated on owned atoms. Reciprocal energy/stress and net-charge terms
are added on rank 0; self-energy partitions over owned atoms; `Reduce`
all-reduces frame energy/stress.

This is correct under domain decomposition. Distributed FFT decomposition is
out of scope for the first charged-NPT acceptance.

Halo still only exchanges immediate neighbour ranks — the charged NPT MPI
fixture uses a box large enough that `cutoff + skin` fits in each subdomain
(see `test_mpi_pme_nacl_parity`, n=6 rock-salt).

## Accepted scope

| In scope | Out of scope (explicit) |
|---|---|
| Isotropic / orthorhombic NPT | Triclinic / anisotropic NPT |
| Neutral ionic crystals | 1-2/1-3/1-4 reciprocal exclusion correction (link 04) |
| Serial + multi-rank replicated grid | Distributed FFT PME |
| Soft core + CoulLong + PME composition | CUDA KSpace (deferred) |

## Public consumers

| Consumer | Path |
|---|---|
| Unit + lifecycle | `tests/cpu/kspace/test_pme_coul.cpp` |
| 1/2-rank parity | `tests/mpi/test_pme_nacl_parity.cpp` |
| Literature regression | `regressions/md/npt_pme_nacl/` |
| LAMMPS validation | `benchmarks/md/npt_pme/validation/bench_pme_nacl_lammps` |

Build with `-DATV_BUILD_KSPACE=ON` (FFTW3). MPI tests also need `-DATV_USE_MPI=ON`.

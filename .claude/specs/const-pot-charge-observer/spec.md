# Spec: `const-pot-charge-observer`

| Field | Value |
| ----- | ----- |
| Slug | `const-pot-charge-observer` |
| Branch | `wf/const-pot-charge-observer` |
| Domain | classical electrochemistry (E1 observe) |
| Method | Electrode total charge observer \(Q_{\mathrm{tot}}(t)\) |
| Chapter | E1 |
| Geometry | Parallel-plate foothold (two electrode type groups) |
| Observable | \(Q_{\mathrm{tot}}\) per electrode [e] |
| Public API | `cpu::ElectrodeChargeObserver` (Observer Component on `md::MD`) |
| Regression | `regressions/md/const-pot-charge-observer/` |
| Oracle | **LAMMPS** (`compute reduce sum` on group charge; ELECTRODE package not required for this sum) |
| Queue row | `docs/developer/metalwalls-method-queue.md` `## in_progress` `[~] const-pot-charge-observer` |
| Map | `.claude/notes/metalwalls-capability-map.md` §1 / §8 electrode charge time series (P0) |

## Goal

Land a **Frame-immutable Observer** that reports total electrode charge
\(Q_{\mathrm{tot},e} = \sum_{i\in e} q_i\) for each electrode group, selected by
atom **type**, for use on `cpu::md::MD` (constant-potential E1 and constant-charge
control). Does **not** solve charges — it only measures `kAtom_charge` after the
charge source (const-pot BO solve, or fixed charges) has written them.

## Non-goals

- Constant-potential BO / CG charge solve (`const-pot-bo-parallel-plate`).
- Freeze / immobile walls (`freeze-electrode-atoms`).
- Charging work energy term (`charging-work-energy`).
- Capacitance / charge-fluctuation analysis (`capacitance-charge-fluctuation`).
- Gaussian electrode kernels, 2D Ewald, Thomas–Fermi, piston.
- New Driver category, public `*Config` bag, or `make_*` / `build_*` factories.
- Welding charge logic into `KSpace*`.
- Live LAMMPS/MetalWalls at ctest time.
- CUDA path.

## Architecture (binding)

- **Kind:** `cpu::Observer` (`FrameImmutable | Deferred`).
- **Bind:** default `md::Phase::end_of_step` (after force / charge update).
- **I/O:** read `kAtom_charge` + `kAtom_type` on owned atoms only
  (`owned_count`); never mutates x/v/f/box/topology/charges.
- **Selection:** one electrode = one atom type id (multi-electrode = vector of
  types, order = electrode index). Parallel-plate foothold: two types.
- **Publish:** member `total_charge(e)` [e]; also write
  `frame.meta` key `kElectrodeCharge` (`MetaF64x3`, first ≤3 electrodes) for
  thermo / write-back.
- **MPI:** local sum over owned atoms, then `MPI_Allreduce(SUM)` when a domain
  communicator is present.

## Tasks

1. Write this spec + `acceptance.md` on `wf/const-pot-charge-observer`.
2. Offline pin via LAMMPS (`tools/metalwalls_oracle/inputs/lammps/const_pot_charge_observer.in`):
   parallel-plate toy — type-1 wall atoms charge \(+q\), type-2 wall \(−q\);
   `compute reduce sum` on each group → pin \(Q_L, Q_R\). Record in
   `regressions/md/const-pot-charge-observer/REFERENCE.md` with
   `oracle: lammps` and `oracle_version`.
3. Implement `cpu::ElectrodeChargeObserver` under `cpu/observe/`; wire CMake
   (MD sources) + unit test under `tests/cpu/observe/`.
4. E2E regression `reg_const-pot-charge-observer` (public headers only): build
   Frame with two electrode types + fixed charges matching the pin geometry;
   run observer (direct launch and/or short MD); assert \(Q\) vs pin within tol.
5. MD-only configure/build (HF/PW/GPW/CUDA/Torch **OFF**); run unit test +
   regression until PASS.
6. Perf gate: fair wall of observer launch vs same LAMMPS `run 0` reduce path if
   practical; else `perf_exempt` with reason (micro-observer, O(N) sum).
7. Commit on `wf/const-pot-charge-observer` only (never `dev`, never push).

## Acceptance criteria

See `acceptance.md` (binding).

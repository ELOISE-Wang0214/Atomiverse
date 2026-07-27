# const-pot-charge-halo-mpi — charge halo correctness after Q update under MPI DD

Status: in-flight (worker, metalwalls loop wave 5)
Branch: wf/const-pot-charge-halo-mpi
Queue row: {"slug":"const-pot-charge-halo-mpi","method":"Charge halo correctness
after Q update under MPI DD","geom":"multi-rank planar walls","observable":
"ghost q match owner after solve","oracle":"literature"}

## Problem

`kAtom_charge` is a `kColComm` column, so the forward halo (`comm::HaloSend`,
`pre_neighbor`) already ships owner charges to ghost rows every step — the
schema/border work is done. But every constant-potential charge solver
(`ConstPotBO`, `ConstPotMaZe`, `ConstPotCG`, `ConstPotentialBO`,
`DipoleElectrodeJoint`) rewrites owned electrode charges at `pre_force @ -50/-500`,
i.e. **after** the forward halo of the same step. Under MPI domain
decomposition, the force stages (`force`: pair Coulomb, kspace) then read
**stale ghost charges** — one solve old, and at step 0 identically zero for
electrode ghosts (initial electrode q = 0). No throw, no crash: forces are
silently wrong exactly at rank boundaries. The capability map (§1/§6, P0)
prescribes: "Prefer solve before halo, or `HaloRefresh` including charge after
solve."

## Design (map-prescribed, house pattern)

Mirror the RATTLE precedent (`comm::HaloRefresh` for positions/velocities +
auto-injection in `MD::install_defaults()` when a Constraint is present):

1. **`comm::HaloRefresh::Field::charges`** — re-send the owner `charge` column
   into the EXISTING ghost slots (reuses `dd->send_lists` / `dd->recv_counts`
   from the last rebuild; idempotent; no PBC shift — a charge is invariant
   under translation). Internals: generalize the refresh column set from a
   hard-wired 3-column triple to an N-column list; charges is a 1-column set.
2. **`Component::updates_owned_charges()`** probe (default `false`) — same
   shape as `requires_system_stress()`. Overridden `true` by the five charge
   updaters above. No name matching, no RTTI (house scan rule).
3. **Auto-inject**: `MD::install_defaults()` under `ATV_HAS_MPI` scans pending
   components; when any `updates_owned_charges()`, inject
   `comm::HaloRefresh{Field::charges}` at `md::Phase::pre_force, +900` — after
   every charge solve (all bind at `pre_force` ≤ -50), before the `force`
   phase. Forgetting a manual refresh line must not produce silently wrong
   rank-boundary forces (same rationale as the constraint refreshes).

Not a new Driver category; not inside KSpace*; Frame stays the sole I/O.

## Out of scope (known, queued separately)

- DD-aware **electrode solve** (each rank's `ConstPotBO` currently collects
  only rank-local electrode/electrolyte atoms → the solved Q on >1 rank is a
  partial-system solve). That is the map's "parallel elastance/CG solve" gap,
  not this slug. This slug guarantees the **transport** invariant: whatever
  charge the owner solved, every ghost copy carries it before forces read it.
- CUDA electrode product, PLUMED, full 2D-PBC stack (loop scope-out).

## Acceptance criteria

1. `comm::HaloRefresh{Field::charges}` exists; `name() == "halo_refresh_q"`;
   single-rank `launch()` is a no-op (unit test, serial).
2. `updates_owned_charges()` returns `true` on `ConstPotBO`, `ConstPotMaZe`,
   `ConstPotCG`, `ConstPotentialBO`, `DipoleElectrodeJoint`; `false` on a
   plain Component (unit test, serial).
3. MPI unit test (`tests/mpi/test_halo_refresh_charges.cpp`, n1 + n2): after
   `HaloSend` forward + owner-charge mutation + `HaloRefresh{charges}`,
   every ghost row's charge equals its owner's current charge (tag-matched,
   exact — the transport is memcpy); ghost positions are untouched by the
   charge refresh; n1 run is a strict no-op.
4. Regression `regressions/md/const-pot-charge-halo-mpi/` (out-of-library,
   public headers only, ≤ 60 s): self-relaunching binary (mpi-pw-scaling
   pattern) runs a 2-rank planar-wall capacitor MD (`ConstPotBO` + PairCoulCut
   + NeighborList + VelocityVerlet) with the auto-injected charge refresh; a
   user-side checker Component (defined in the regression file) verifies at
   the force stage of every step that max |q_ghost − q_owner| over all ghost
   rows and ranks is 0 (≤ 1e-15) AND that electrode ghosts with |q| > 0 exist
   (non-vacuous). Ships `REFERENCE.md` with oracle: literature + citations.
5. Dual gate green: scoped MD-only+MPI build; ctest on the new unit tests;
   `run_reg_const-pot-charge-halo-mpi` passes.
6. perf: no fair external twin binary for an MPI-halo self-consistency
   invariant → `perf_exempt=true`.

## Oracle

`oracle=literature` / self-consistency (queue row). The invariant is the
LAMMPS forward-communication contract for constant-potential electrodes: after
the charge solve, updated charges must be forward-communicated to ghosts
before force evaluation (Thompson et al., Comput. Phys. Commun. 271, 108171
(2022); Ahrens-Iwers et al., J. Chem. Phys. 157, 084801 (2022) — fix electrode
issues a forward_comm of q after every conp solve; Marin-Laflèche et al.,
J. Open Source Softw. 5, 2373 (2020) — MetalWalls). No oracle binary is run:
the pinned number is the exact-transport bound 0 (≤ 1e-15). Pin notes:
`tools/metalwalls_oracle/pin_notes_const-pot-charge-halo-mpi.md`.

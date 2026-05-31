# MPI halo + migration refactor proposal — 2026-05-28

## Why this exists

The current `cpu::md::mpi` subsystem (~2700 LOC across `domain.{h,cpp}`,
`halo_exchange.{h,cpp}`, `atom_migration.{h,cpp}`, `gather.cpp`, `rcb.cpp`,
`reduce.cpp`) shipped working code (commit `d8eb8a9` finally lets MSU PE
run end-to-end at 8 ranks), but the design has accumulated several
band-aids that make every new bug look like the last one. The 2026-05-27
stage-2 conserved-drift incident is the clearest example: the root cause
was *not* a halo bug, *not* a migration bug, but a **handoff bug between
the two** that exists only because they share a mutable `dd->ghost_tags`
log.

This note inventories the deficiencies, contrasts them point-by-point
with LAMMPS' `Comm`/`CommBrick`/`Atom`/`AtomVec` design, and proposes a
refactor that closes the deficiencies without touching the user-facing
component API.

---

## Defects in the current design

### D1. `cpu::Frame` is fixed-size — no `grow_arrays()`

`Frame::Frame(int n)` allocates exactly `n` per-atom rows. `Atoms` has
no `resize()` or `reserve()`. When the halo cascade overflows
(`slot + new_ghosts > frame.atoms.tag.size()`), `halo_exchange.cpp:493`
calls `MPI_Abort`. **The capacity is a user-tuned input.**

That is why the freshly-added `test_mpi_nose_hoover_npt_multirank` blew
up when its 80-atom 50Å box yielded a cascade of ~95 ghosts on one rank;
why the previous `partition_atoms` doc note says "Stale rows in
`[n_local, capacity)` are reused later by HaloExchangeForward as ghost
slots"; and why we hand-rolled a separate `dd->capacity` field on top of
`frame.atoms.size()`.

**Contrast (LAMMPS):** `Atom::nmax` is the allocated capacity and grows
on demand. `AtomVec::unpack_border()` and `unpack_exchange()` invoke
`avec->grow()` when `nlocal + nghost > nmax`. `Comm::borders()` never
pre-flights the count.

### D2. `dd->ghost_tags` log is a workaround for the wrong problem

`tag_to_local` mixes owner and ghost entries. The halo full-rebuild
path keeps a sidecar log (`ghost_tags`) of "tags I previously
inserted as ghosts" so it can `erase` them from `tag_to_local`
before reinserting the current step's ghosts.

This created the 2026-05-27 stage-2 incident: `AtomMigration` rebuilds
`tag_to_local` from current owners but originally did *not* clear
`ghost_tags`. Subsequent halo `for (t : ghost_tags) erase(t)` then
silently deleted owner entries that had migrated in *and* were ghosts
on this rank the previous step. The bonded ownership filter then
dropped ~15 bonds per migration. `dd->ghost_tags.clear()` in
`atom_migration.cpp:334` (commit `d8eb8a9`) is patch tape over a deeper
issue: there should not be a ghost-tags log in the first place.

**Contrast (LAMMPS):** `Atom::map_clear()` is called once at the
*start* of `exchange()`. After `exchange()` + `borders()` finish, a
single `Atom::map_set()` rebuilds the map from `[0, nlocal+nghost)`.
Owners and ghosts are entered exactly once per neighbor cycle.
There is no incremental tracking, no undo log, no "what did I
insert last step" memory. Correct by construction.

### D3. Two generation counters chase each other (`neighbor_generation` vs `halo_generation`)

`dd->neighbor_generation` is bumped by `AtomMigration` and read by
bonded ownership filters and the NL. `dd->halo_generation` is bumped by
`HaloExchangeForward` full-rebuild and read *only* by the NL to know
when ghost slot indices have shuffled.

The comment block on `halo_generation` (domain.h:82–99) explicitly says
bumping `neighbor_generation` from the halo "interacts disastrously
with [`AtomMigration::first_launch`] handling on a Zarr-loaded state."
That comment is the smoking gun: we have two counters because one
component's contract leaked into another's state.

**Contrast (LAMMPS):** One signal — `neighbor->ago == 0` after every
reneighbor (which is what triggers a `borders()` call). `borders()`
itself rebuilds the map, the send-lists, and the ghost rows in one
shot. Downstream caches key off `ago == 0` (a per-step boolean), not
a free-running counter.

### D4. `partition_atoms()` is a duplicate of `exchange()`

`DomainDecomposition::partition_atoms()` is "called once on first bind";
`AtomMigration::launch()` is "called every K steps thereafter". They
do nearly the same work — classify atoms by spatial position, ship
them to the destination rank, rebuild `tag_to_local` — but use
different code paths, different MPI patterns (`partition_atoms` uses
local-only filtering on the global frame; `AtomMigration` uses
`MPI_Alltoallv`), and different post-state contracts (`partition_atoms`
sets `dd->capacity`; `AtomMigration` does not).

The `backend.set_frame()` carve-out (`backend.cpp:127`) — "first call
runs partition_atoms; subsequent calls only refresh sublo/subhi" —
exists because we couldn't unify the two.

**Contrast (LAMMPS):** `Comm::exchange()` is called every time the
neighbor list is rebuilt. The very first call (after `Atom::data_atoms`
has read the data file on rank 0 and broadcast/scattered to all ranks)
is structurally identical to step-N migration. There is no
"partition_atoms" duplicate.

### D5. `MPI_Abort` on capacity overflow is wrong response

`halo_exchange.cpp:493–500` aborts the world on cascade overflow. The
right response is **grow the arrays**. We don't, because of D1.

### D6. Send-list rebuild scans frame positions 3 times per axis

`rebuild_send_lists_for_axis` is called 3 times (once per axis dim) in
the cascade. Each call iterates `[0, sources_end)` and writes to
`dd->send_lists[2*dim+0/1]`. The pattern of scanning a growing prefix
is correct, but the function takes `sources_end` as a *separate*
parameter from `dd->n_local + ghosts_so_far`. The caller
(`HaloExchangeForward::launch`) tracks the running `slot` variable and
passes it down; the function itself is unaware of which axis stage it
is. That coupling makes the axis-cascade invariant invisible from
either side.

**Contrast (LAMMPS):** `for dim in [0,1,2]: for ineed in [0,1]:` —
the cascade is a 2-level for-loop directly inside `borders()`.
`nlast = atom->nlocal + atom->nghost` is recomputed inline at the
boundary of each swap. No helper function with a hidden contract.

### D7. Atom payload field list duplicated in 3 places

`atom_migration.cpp:42–49` hardcodes `kPayloadDoubles=10`,
`kPayloadInt32=1`, `kPayloadInt64=1`, `kPayloadInt8=3`. The `move_row`
helper a few lines below lists the 14 field assignments by name. The
halo `pack_positions_for_face` lists the 4-double layout (x, y, z, tag)
inline. Adding a new per-atom field (e.g. per-atom charge, per-atom
stress tensor for stress profiles, per-atom virial — there's a TODO
about this in halo_exchange.h:191–199) requires editing all three
places, plus the new ghost-slot rewrite in halo, plus the new
exchange pack/unpack. No central abstraction.

**Contrast (LAMMPS):** `AtomVec::pack_border_vel`, `pack_exchange`,
`pack_comm` are *member functions* of the atom style. Each atom style
(atomic, full, molecular, charge, dipole, …) defines what its
ghost/exchange/comm payload is in one place. New fields → extend the
one class.

### D8. `MD::install_infra()` auto-injects halo with `1.1 * pair_cutoff`

The 10% skin is a global constant baked into `md.h:177`. The comment
there acknowledges callers with a custom NL skin > 10% have to
register `HaloExchangeForward{halo_cut}` *explicitly* before
`set_pair()`. There's no skin-detection: the halo cutoff is fixed
relative to the pair cutoff, ignoring the NL skin actually configured
later via `set_neighbor`.

**Contrast (LAMMPS):** `comm_modify cutoff` (or its automatic
counterpart `Neighbor::oneatom * skin` propagated to comm) sets the
ghost cutoff to **`cutneigh = max(force_cutoff, bonded_cutoff) +
neighbor_skin`**. The skin is the *neighbor list*'s skin, not a
separately invented halo-only constant.

### D9. The `every=10` reneighbor cadence is per-halo-component, not per-driver

`HaloExchangeForward::every_ = 10` (constructor default, halo_exchange.h:93).
`NeighborList::every` is separate (defaults to whatever `set_neighbor`
takes). These should be the same number by definition — every halo full
rebuild *is* a reneighbor event. We let them drift independently.

### D10. Bonded ownership filter (`select_owned_tuples`) is a per-step
filter on top of LAMMPS-style storage

Every bonded force kernel (Bond / Angle / Dihedral) calls
`select_owned_tuples` which iterates *all* bonds, resolves each atom
tag through `tag_to_local`, and keeps the bond only when
`min(tag_i, tag_j)` resolves to a slot `< n_local`. This works but
duplicates topology work each call. LAMMPS stores
`atom->num_bond[i]` and `atom->bond_atom[i][m]` *per-atom*, so the
owner is implicit (atom `i` owns its bonds with index `m`). No
filtering, no map lookup in the hot loop.

This is more invasive to refactor (touches Frame's topology layout),
so I'd defer it to a Phase 2 of the refactor. The fix from this note
removes the *coupling* between bonded filter correctness and ghost-tag
log correctness, which is the immediate pain point.

---

## Proposed refactor (Phase 1 — minimal viable)

Goal: replace D1–D6 with LAMMPS-shaped equivalents while keeping the
public component API (`HaloExchangeForward`, `HaloExchangeReverse`,
`AtomMigration`, `DomainDecomposition`) intact. Bonded refactor (D10)
deferred to a Phase 2 because it changes the data model.

### Step 1. Make `cpu::Frame` genuinely growable

Add `Frame::Atoms::resize(int new_nmax)` that calls `xt::xtensor::resize`
on every per-atom array. Add `Frame::Atoms::nmax()` that returns the
current allocated capacity (distinct from `size()` which now returns
`n_local`).

The semantic split becomes:
- `frame.atoms.size()` — **count of live owners** (= `dd->n_local`).
  Loops over "do this for every owned atom" use this.
- `frame.atoms.nmax()` — **allocated capacity** (≥ `n_local + n_ghost`).
  Halo and migration loops over ghost slots use this when bounds-checking.

`frame.atoms.tag.size()` (xtensor's intrinsic size) becomes equal to
`nmax`, which is what `halo_exchange.cpp:570` already reads. The check
at line 492 stays, but instead of `MPI_Abort` it calls
`frame.atoms.resize(round_up(slot + new_ghosts, growth_factor))`.

Closes D1, D5.

### Step 2. Drop `dd->ghost_tags` and rebuild `tag_to_local` lazily

Remove the `ghost_tags` vector from `DomainDecomposition`. The
`tag_to_local` map becomes *owners only* after `AtomMigration` and
*owners + ghosts* after `HaloExchangeForward` full-rebuild.

Rebuild contract:
- `AtomMigration::launch` ends with `rebuild_tag_to_local(dd, frame,
  [0, n_local))` (owners only).
- `HaloExchangeForward::launch` full-rebuild ends with `dd->n_ghost = slot - dd->n_local; rebuild_tag_to_local(dd, frame, [0, n_local + n_ghost))` (owners + ghosts).
- `HaloExchangeForward::launch` position-only-refresh does **not** touch
  `tag_to_local` (ghost slot indices are stable between rebuilds — same
  contract as LAMMPS `forward_comm`).

`rebuild_tag_to_local` is a single helper:
```cpp
void rebuild_tag_to_local(DomainDecomposition& dd, const Frame& frame,
                          int begin, int end) {
    dd.tag_to_local.clear();
    dd.tag_to_local.reserve(end - begin);
    for (int i = begin; i < end; ++i)
        dd.tag_to_local.emplace(frame.atoms.tag[i], i);
}
```

Closes D2, D3 (since `halo_generation` was only there to compensate for
`tag_to_local` ghost-slot churn; now the map is fully rebuilt at every
full halo, so consumers can just read `dd->neighbor_generation`).

### Step 3. Collapse `partition_atoms` into `AtomMigration::first_launch`

Replace `DomainDecomposition::partition_atoms` with a single
`AtomMigration::launch_initial(ctx)` path that does the same logic as
the steady-state path but enters with `dd->n_local == 0` and the *full
global atom set* on every rank (the existing input contract for the
multi-rank entry).

`backend.set_frame` then just sets `frame_` pointer and calls
`compute_sublo_subhi`; the migration component takes care of the
initial partition on its first launch.

The `dd->n_local == 0` first-launch branch is a 5-line classifier-loop
identical to what's in `partition_atoms.cpp` today; the rest of
`AtomMigration::launch` already handles the "rebuild compacted owner
prefix, append received atoms, rebuild tag_to_local" steps.

Closes D4. Also removes the awkward "subsequent set_frame skips
partition" comment block in backend.cpp.

### Step 4. Unify the cascade into a single function with no helper-arg coupling

Move the body of `rebuild_send_lists_for_axis` and
`forward_axis_stage` inline into `HaloExchangeForward::launch`, written
as one nested for-loop that exactly mirrors LAMMPS:

```cpp
int nfirst = 0, nlast = dd->n_local, slot = dd->n_local;
for (int dim = 0; dim < 3; ++dim) {
    for (int direction = 0; direction < 2; ++direction) {  // -dim, +dim
        // Rebuild send_lists[2*dim + direction] by scanning [nfirst, nlast)
        // MPI Isend/Irecv counts then positions+tags
        // Unpack into [slot, slot + recv) and grow frame if needed
        // slot += recv
    }
    nfirst = 0;
    nlast = slot;  // include all ghosts so far for next axis
}
dd->n_ghost = slot - dd->n_local;
rebuild_tag_to_local(*dd, frame, 0, slot);
```

Closes D6. Also makes the cascade invariant ("next axis scans all
prior ghosts") visible without indirection.

### Step 5. Pull halo cutoff from NL skin instead of a hardcoded 1.1 factor

Drop the 10% multiplier in `MD::install_infra`. Replace with:
- Read the registered `NeighborList`'s `(cutoff, skin)` after
  `scan_pending_infra_` (`s.neighbor_skin = ...`).
- Auto-injected halo cutoff = `s.pair_evaluator_cutoff + s.neighbor_skin`.
- If no NL is registered, fall through to the throw-on-missing-NL path
  that already exists.

Closes D8 and D9 (now there's one skin, owned by the NL).

### Step 6. Drop the `ago_ / every_` cadence from `HaloExchangeForward`

Replace with a "consumed from NL" model: the NL decides each step
whether a reneighbor is needed (it already does — `is_stale()` based on
max displacement vs skin/2). On reneighbor-needed steps, the NL bumps
`dd->reneigh_request_step`. The halo reads that flag at `launch` time
and runs full-rebuild iff it's been bumped this step.

Halo and NL now share the same trigger; no per-component every_
parameter to drift apart. The signal lives on `dd` (or
`ctx.runtime.outer_status`), not on either component's private state.

Closes D9 properly (D8 was already closed in Step 5).

---

## Out of scope for this refactor

- **D7 (atom payload duplication)** — needs an `AtomVec`-style
  abstraction. Worth doing later; not on the path to "MSU runs."
- **D10 (bonded ownership filter via per-atom topology)** — needs to
  reshape `Frame::Bonds` from `{atomi[], atomj[], type[]}` to a
  per-atom-indexed `{num_bond[i], bond_atom[i][m], bond_type[i][m]}`
  layout. Substantial.
- **CommTiled / RCB unification** — `rcb.cpp` is already a separate
  component; leaving it alone.

---

## Validation plan

After the refactor lands:

1. All 7 functional MPI tests must still pass under
   `unset OMP_*; srun -n N -c 1 build/tests/...`. (The pre-existing
   `FaceBufsAreInstanceOwned` test is incidental — its setup probes
   `dd->capacity` re-binding, which the Step 1 (grow_arrays) change
   makes the question moot. Plan to delete that test rather than fix
   it.)
2. MSU PE stage-2 conserved drift must stay < 1% over 2000 steps at
   8 ranks (currently -1.35% with the patch).
3. MSU PE full 5-stage pipeline must complete cleanly at 8 ranks at
   scale=0.01 (currently does — 940s).
4. Add a stress test: 100k atoms in a box where halo zone is ~80% of
   cell volume (forces frequent autogrow). Confirm no MPI_Abort.

## File-by-file edit estimate

- `include/atomiverse/cpu/frame.h`: +30 LOC (Atoms::resize, Atoms::nmax,
  separation of `size()` vs `nmax()`).
- `src/cpu/frame.cpp` (if it exists, else inline): +20.
- `include/atomiverse/cpu/md/mpi/domain.h`: −20 (remove `ghost_tags`,
  `halo_generation`, `capacity`; capacity comes from `frame.atoms.nmax()`).
- `src/cpu/md/mpi/domain.cpp`: −50 (delete `partition_atoms`).
- `src/cpu/md/mpi/halo_exchange.cpp`: net −150 (inline the cascade,
  delete the helper, delete the position-only-only path — fold into
  one body with a `full_rebuild` branch only where the payload differs).
- `src/cpu/md/mpi/atom_migration.cpp`: +30 (add first-launch branch),
  −15 (remove `ghost_tags.clear` since the field is gone).
- `src/cpu/backend.cpp`: −20 (delete the partition-vs-rebind carve-out).
- `include/atomiverse/cpu/md/md.h`: +5 −5 (Step 5).

Net: ~−200 LOC, more contracts implicit in structure rather than spelled
out in comments.

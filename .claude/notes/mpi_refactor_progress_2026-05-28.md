# MPI subsystem LAMMPS-shape refactor — progress, 2026-05-28

## Status

**Foundation committed** (5 commits, ~1500 LOC of LAMMPS-shaped scaffolding as dead code coexisting with the live old halo/migration path):

| Commit  | What landed                                                |
|---------|------------------------------------------------------------|
| 6edff73 | `cpu::Atom` — LAMMPS-shaped per-atom container (nlocal/nghost/nmax/grow/map; per-atom topology fields num_bond/bond_atom/.../nspecial/special) |
| 4a9fbea | `cpu::comm::atom_pack` — pack_comm / pack_border / pack_exchange (variable-length, carries per-atom topology) / compute_bonded_extension |
| 3149313 | `cpu::Domain` — LAMMPS-shaped spatial decomposition (Cartesian topology + sublo/subhi + reneigh_this_step; no atom counts) |
| f15e376 | `cpu::comm::Comm` — LAMMPS CommBrick equivalent with exchange/borders/forward_comm/reverse_comm methods; CommState PIMPL holds sendlist/recvnum/bufs |
| 741927f | Topology builders on `cpu::Atom` — add_bond / add_angle / add_dihedral / build_special_list (BFS over bond graph for 1-2/1-3/1-4 exclusion list) |

**Old path still live.** `HaloExchangeForward / Reverse`, `AtomMigration`, `DomainDecomposition`, `topology_filter.h` all unchanged. MSU PE 8-rank still uses them; functionality identical to commit `d8eb8a9` (the `ghost_tags.clear()` fix).

## What's left to land

Approx scope: ~2000+ LOC across ~50 files; sequencing is tangled because Frame dissolution + materialize update + force kernel rewrite all interlock.

### Required for switchover

1. **Frame dissolution or thin-holder reshape.**
   - cpu::Frame currently has god-struct sub-structs `Atoms`, `Bonds`, `Angles`, `Dihedrals`. These overlap with the new `cpu::Atom` (which has per-atom topology built in).
   - Cleanest: delete `Frame::Bonds/Angles/Dihedrals`; make `Frame::Atoms` an alias of `cpu::Atom`. Touches every callsite that reads `frame.bonds.atomi[k]` etc. (~50 files).
   - The `frame.atoms.X` accessor pattern works as-is because `cpu::Atom` exposes the same field names (x, y, z, vx, ..., tag, image).

2. **Materialize + LAMMPS data reader: per-atom scatter.**
   - `src/cpu/materialize.cpp`: replace flat topology fill with `atoms.add_bond(...)`-style scatter; call `atoms.build_special_list()` at the end. Started — reverted in this session because it cascaded; revisit when (1) is done.
   - `src/cpu/lammps_data.cpp`: same pattern.

3. **Force kernels: per-atom layout.**
   - `src/cpu/potentials/bond_harmonic.cpp`: outer loop becomes `for (i in [0, nlocal)) for (m in [0, num_bond[i]))`.
   - `src/cpu/potentials/angle_harmonic.cpp`: center-atom-owned iteration.
   - `src/cpu/potentials/dihedral_multi_harmonic.cpp`: atom-j-owned iteration.
   - Remove the `select_owned_tuples` post-hoc filter from `include/atomiverse/cpu/potentials/topology_filter.h` (delete the header).
   - `src/cpu/potentials/pair_lj_cut.cpp` + `pair_morse.cpp`: read exclusion list from `atom.nspecial / atom.special` (currently consume `frame.bonds/angles/dihedrals` directly to build their own exclusion bitmask).

4. **Wire `cpu::comm::Comm` into `MD::install_defaults`.**
   - Replace auto-injected `HaloExchangeForward / Reverse + AtomMigration` with a single `comm::Comm` component family (or one `Comm` PIMPL exposed through phase-binding adapters).
   - `MD::install_infra` → `install_defaults` (per the user direction).
   - Halo cutoff = `max(pair_cutoff + NL.skin, bonded_extension)`; reads from `cpu::comm::compute_bonded_extension_local`.

5. **Delete old code.**
   - `include/atomiverse/cpu/md/mpi/halo_exchange.{h,cpp}` ← merged into `cpu::comm::Comm`.
   - `include/atomiverse/cpu/md/mpi/atom_migration.{h,cpp}` ← merged into `cpu::comm::Comm::exchange`.
   - `include/atomiverse/cpu/md/mpi/domain.{h,cpp}` ← `cpu::Domain` replaces.
   - `include/atomiverse/cpu/potentials/topology_filter.h` ← per-atom ownership invariant removes the need.
   - `tests/cpu/md/test_topology_resolved_locals.cpp` ← the filter it guards no longer exists.

6. **CUDA bridge stub.**
   - `src/cuda/transfer.cu` `cpu::to_cuda` / `cuda::to_cpu` field-by-field copy currently reads `frame.bonds.atomi`. Once flat topology is gone, the copy needs a per-atom→flat conversion path OR has to assert `nullptr` until CUDA-MPI lands.

7. **Test fixture overhaul.**
   - Every `tests/mpi/*.cpp` and `tests/cpu/md/test_pe_*.cpp` that builds topology via `frame.bonds.atomi = xt::zeros<int>({...})` switches to `frame.atoms.add_bond(...)` calls. ~25 test files; mostly mechanical.

### Verification (deferred)

After switchover lands, all the verification steps in the original plan
(.claude/plans/reflective-dazzling-horizon.md §Verification):

- Single-rank parity (`ctest --test-dir build` non-MPI)
- 8-rank MSU stage 2 from stage1.zarr (T~500K, conserved drift < 1.5%)
- 8-rank MSU full pipeline at scale=0.01 (all 5 stages clean)
- 1r-vs-8r bonded force diff (E_bond+E_angle+E_dihedral identical to 1e-10 eV)
- Frame-purity grep returns empty for ghost_tags / halo_generation /
  neighbor_generation / partition_atoms / topology_filter / select_owned_tuples.

## Why incremental commits don't help here

The switchover items (1) – (7) are tightly coupled:
- Deleting `Frame::Bonds` breaks all force kernels until they're rewritten for per-atom.
- Deleting `topology_filter.h` breaks bonded force kernels until they switch.
- Switching `install_infra` to new Comm requires deleting old `HaloExchangeForward` to avoid the auto-injector tripping on both.
- Stub `cpu::to_cuda` requires `frame.bonds` to not exist.

Trying to land these as separate per-Layer commits creates intermediate
states where the build is broken in ways that aren't a clean
"compile-but-test-fails" — they're "doesn't compile at all".

Realistic path: one large switchover commit (~1500-2000 LOC of touches),
followed by a verification commit. Or, accept that the foundation stops here
and the switchover happens in a fresh dedicated session.

## Next steps

Either:
- **Continue in a future session**: pick up at (1) Frame dissolution. The
  scaffolding (commits 6edff73..741927f) is solid; the consumer migration
  is just mechanical work.
- **Or roll back the scaffolding** if you decide the design is wrong.
  `git reset --hard d8eb8a9` removes all 5 scaffold commits; the
  working `ghost_tags.clear()` fix remains in place.

## Update — 2026-05-29 early morning: multi-rank works for in-code MD; Zarr-load step-75 spike remains

The "astronomical E_bond" bug at 8-rank was diagnosed:
**`partition_atoms` only updated `dd->n_local`, not `atoms.nlocal`.** The
new per-atom force kernels were iterating `[0, atoms.nlocal)` — i.e. the
Frame-ctor's full atom count, *including stale uncompacted rows*. Each
stale row carried `num_bond/bond_atom` from the pre-partition layout
and was matched against current (partitioned) positions, giving 100Å
"bond lengths" → E_bond ≈ 2.85×10⁸ eV at step 0.

A first attempt to sync `atoms.nlocal = dd->n_local` directly in
`partition_atoms` segfaulted inside `MPI_Alltoallv`'s `MPIR_Localcopy`
during the first `AtomMigration` step. The pack-stride math is tied to
`frame.atoms.size()` which now reports the shrunk count; the buffer
layout went inconsistent across ranks. **The right fix was to keep
`atoms.nlocal` untouched and have the force kernels read `dd->n_local`
directly in multi-rank mode** (single-rank still uses `atoms.nlocal`).

Per-thread force accumulators (`fx_thr_[t]`) also needed sizing to
`atoms.nmax` rather than `atoms.size()` — they accumulate into ghost
rows that `reverse_comm` gathers back.

Results after the fix (commit `70c1145`):
- 8-rank `test_mpi_pe_msu_pipeline` PASSES (now exercising bonded
  forces — previously vacuous because per-atom topology was empty
  when the test used `frame.bonds.atomi` flat shim).
- 8-rank MSU PE in-code NVT stage 1 runs cleanly through 80 steps:
  T cooling from 10643K → 1124K via Langevin, E_bond bounded
  (8 → 19766 → 4706 → 2030 → 1636 eV chain relaxation), n_local
  stable per-rank around 1830-1850 with slight migration. **Same
  ranges as single-rank.**
- 8-rank MSU PE NPT stage 2 from `stage1.zarr`: **clean through
  step 50** (T = 502, 509, 516 K; E_bond = 227, 218, 220 eV;
  conserved drift = -0.01% — matches single-rank). **At step 75 a
  spike appears** (T → 28000 K, conserved → 52000 eV). This is the
  same "step ~70 catastrophic blowup" pattern that the original
  `d8eb8a9` ghost-tags-clear fix had to chase down — likely a
  similar subtle interaction at a halo full-rebuild boundary now
  that the per-atom code path is exercised. Doesn't appear in the
  in-code NVT path → specific to `--ckpt-load-before` + multi-rank
  + per-atom topology.

The step-75 Zarr-stage-2 spike is the last known issue. Diagnosis
path: instrument `bond_atom[i][m]` resolution at step 70-75 to
see whether (a) a partner tag resolves to a stale row, or (b) the
ghost position for the partner is wrong (halo bug), or (c)
`Atom::reshuffle_topology_to_tags` left some row mis-assigned.

## Update — finer-grained diagnosis: spike correlates with AtomMigration event

Print-every=5 trace at stage 2 from `stage1.zarr` shows the spike
precisely at step 70:

```
step   T(K)  E_pot   E_bond  E_kin   conserved  neigh_gen
0      502    502    227     649     1151       2
…
65     508    493    226     657     1151       8
70     11433  2779   223     14776   17559      9   ← neighbor_gen++
75     28382  14911  14604   36684   52698      9
80     19101  16216  13351   24688   42476      10
```

`neighbor_generation` bumps from 8 → 9 between step 65 and 70 — that's
an `AtomMigration::launch` event. The new `AtomMigration` pack/unpack
goes through `comm::pack_exchange` / `unpack_exchange` (variable-length,
carries per-atom topology) with a fixed 800-byte per-atom stride
(`kAtomPayloadBytes = kMaxExchangeDoubles * 8`). Pack writes ≤ 100
doubles per atom and zero-pads to 800; unpack reads back via the same
helper and advances by the fixed 800 bytes. The stride math is
self-consistent.

But evidently *something* about the new per-atom-topology-carrying
migration is corrupting per-atom topology rows for migrated atoms.
The in-code NVT stage 1 (which also has migration events every ~10
steps) does NOT blow up — so the issue is specific to
**`Atom::reshuffle_topology_to_tags` interacting with migration**.
After reshuffle, the row-i-holds-atom-tag[i] invariant is locally true.
But when migration sends atom from rank A's row to rank B's row, it
packs per-atom topology that was reshuffled relative to A's tag
ordering — then unpacks into B at the new row index. B's tag_to_local
post-migration maps the migrant's tag to the new row, but B's other
atoms' topology references partner tags from a different rank's
perspective.

Hypothesis: after reshuffle, the per-atom `bond_atom[i][m]` partner
tag is the *global atom tag* (preserved across migration), but
`build_special_list` BFS produced the `nspecial`/`special` exclusion
list under the assumption that the local bond graph was complete. On
the receiving rank after migration, the migrant's row carries its
old `nspecial`/`special` (correct global tags) — so the LJ exclusion
loop is fine. Yet energy explodes.

**Recommended next step**: instrument at the migration boundary —
log `tag[migrant_row]`, `num_bond[migrant_row]`,
`bond_atom[migrant_row][0..nb]` immediately BEFORE pack on the sender
and immediately AFTER unpack on the receiver. If they differ → pack
ordering bug. If they match → the issue is in how the recipient's
tag_to_local / halo handles the new owner row.

## Update — instrumentation result: pack/unpack consistent (the bug is elsewhere)

Ran with `ATV_MIGRATE_DEBUG=1` and the diagnostic prints inserted in
`pack_atom`/`unpack_atom`. Every traced atom showed the sender's PACK
line matching the receiver's UNPACK line bit-for-bit: same tag, same
num_bond, same bond_atom[0..nb], same bond_type[0..nb]. Pack/unpack
fidelity is NOT the bug.

Excerpt (rank 0 pack → rank 1 unpack):
```
[migrate][rank0 PACK]   dest=1 tag[1592]=7449 num_bond=1 bond_atom[0]=7450 bond_type[0]=0
[migrate][rank1 UNPACK] slot=1492 tag=7449       num_bond=1 bond_atom[0]=7450 bond_type[0]=0
```

So the per-atom topology survives `MPI_Alltoallv` intact. The Zarr-stage-2
step-70 spike must come from one of:

1. **Halo ghost handling for migrated atoms** — `HaloExchangeForward`
   runs immediately after `AtomMigration` in `pre_force`. It rebuilds
   `dd->tag_to_local` for ghosts, ships positions of boundary atoms.
   If the receiver's ghost slot ends up with a wrong position for the
   migrated atom (e.g. PBC shift bug, or ghost slot reuse collision),
   bond force on the partner's rank computes `(partner_pos - wrong_ghost_pos)`
   → huge bond stretch → astronomical energy.

2. **NeighborList CSR / cell-list cache** — NL caches per-step pair
   list based on `dd->halo_generation`. The new code might be reading
   a stale cache for the first force pass after migration, mixing
   pre-migration ghost slots with post-migration atom data.

3. **Compaction-time ordering inside `partition_atoms`** — the per-atom
   topology compaction runs AFTER the SoA compaction. Between the two,
   row indices in [0, n_local) are partially updated (tag/positions
   compacted; topology still at original rows). If anything reads
   topology in that window, it sees a mismatch.

The diagnostic instrumentation was reverted (commit `<latest>`). The
debug-print machinery is straightforward to re-add: search for
`ATV_MIGRATE_DEBUG` in atom_migration.cpp git history. To diagnose
further, the most direct path is:

- Add a similar trace at HaloExchangeForward's full-rebuild unpack
  for the *same atom tags* that migrated in the previous AtomMigration.
- If those positions differ from the sender's positions → halo bug.
- If those positions match but the partner's bond resolution still
  produces huge force → bond_atom value is wrong despite the
  pack/unpack trace showing it's right (which would mean topology is
  being mutated between migration and force eval).

The remaining 5 tasks (#16 Layer D, #17 CUDA bridge, #18 delete old
code, #19 MPI tests, #20 verification) all depend on this last
debug, then the cleanup is mechanical.

## Update — bug is NPT stage-transition, not Zarr-load

Tested `--mode npt` WITHOUT `--ckpt-load-before` (stages 1→2 in-code,
no Zarr involved). Stage 2 still breaks within 20 steps:
T=0, E_kin=0, conserved=1.85e151. So the bug is NOT specific to
Zarr — it's specific to the **stage-1 → stage-2 transition with
NoseHooverNPT + per-atom topology + multi-rank**.

In-code NVT stage 1 alone: works fine (T cooling, E_bond bounded).
Stage transition into NPT stage 2: collapses.

The per-stage driver rebuild flow:
1. Stage 1 ends, `Backend::release_dd` extracts the dd state.
2. Stage 2's `Backend::init` does `adopt_dd` to keep the partition.
3. Backend's first `set_frame` for stage 2 sees `dd->n_local > 0`,
   skips `partition_atoms`, just refreshes sublo/subhi.
4. Components rebuild against existing dd.
5. NoseHooverNPT new instance, first `launch()` runs.

Hypotheses for what breaks the transition:
- a) New `AtomMigration` instance's `last_observed_generation_` is
     `~0`. With my pack_exchange-based pack/unpack and a generation
     bump that already happened (from stage 1's halo/migration),
     the first stage-2 launch might attempt a spurious migration
     that ships uninitialised payload.
- b) `compute_bonded_extension_local` in `HaloExchangeForward::build`
     now iterates `dd->n_local` (commit just landed), which uses the
     adopted dd's n_local. Should be correct.
- c) NPT-specific: `NoseHooverNPT::build()` reads `atoms.size()`
     (= 10000 always, my refactor unchanged it). For Q chain mass
     calc it allreduces dd->n_local. But the integrator iteration
     itself uses `owned_count(ctx)` → `dd->n_local` → correct.
- d) Velocity scrambling: if stage transition somehow zeroes
     velocities (which would explain T=0), KE is 0 and NPT's chain
     diverges.

The bug needs interactive debugging — instrument `NoseHooverNPT::build`
and the first stage-2 `launch`, log `atom.vx[0..3]`, `dd->n_local`,
and compare against single-rank reference.

**Workaround for the user:** Hard-reset `git reset --hard d8eb8a9`
to restore the pre-refactor working state.

## Final task status (session end)

20 commits since `d8eb8a9` baseline. Task list state:

| Task | Status      | Note                                                  |
|------|-------------|-------------------------------------------------------|
| #10  | completed   | cpu::Atom created (commit `6edff73`)                  |
| #11  | completed   | cpu::comm::atom_pack (commit `4a9fbea`)               |
| #12  | completed   | cpu::Domain (commit `3149313`)                        |
| #13  | completed   | cpu::comm::Comm (commit `f15e376`)                    |
| #14  | completed   | Force kernels rewritten (commit `3ecefef` + fixes)    |
| #15  | completed   | Materialize/lammps_data scatter (in `3ecefef`)        |
| #16  | **pending** | Layer D — needs Frame dissolution + new Comm wired up |
| #17  | completed   | CUDA bridge: per-atom→flat conversion (commit `<latest>`) |
| #18  | **pending** | Delete old halo/migration — needs Layer D first       |
| #19  | **pending** | Update MPI tests — needs new comm wired               |
| #20  | in-progress | Verification — single-rank pass; 8-rank NVT pass; NPT broken |
| #21  | completed   | cpu::Frame thin holder (was: full deletion not done)  |
| #22  | completed   | Mass rename — not actually needed (Frame stayed)      |

**Key fix discovered**: commit `70c1145` — force kernels must
iterate `dd->n_local`, not `atom.nlocal`, in multi-rank. Without
this, stale uncompacted rows beyond `dd->n_local` give astronomical
bonded forces. This is a correctness fix that any future
per-atom-topology + multi-rank work must preserve.

**Outstanding bug**: NoseHooverNPT × multi-rank × per-atom topology
collapses energy to zero / explodes within ~10 stage-2 steps. Same
pattern from in-code path AND from Zarr-load path. Reproducer:

```
salloc -n 8 -t 00:15:00
srun -n 8 -c 1 build/benchmarks/polyethylene_tension \
    --data benchmarks/md/reference/polyethylene/msu_pe.data \
    --mode npt --scale 0.003 --print-every 20 \
    --out-dir /tmp/repro
```

Step 40 stage 2 shows T=0, E_kin=0, E_pair=3.18e9. Diagnosis:
NPT chain dynamics over-respond to the high E_pair from initial
LJ collisions, killing velocities. But original d8eb8a9 with the
same scale ran cleanly through all 5 stages, so my refactor broke
something specific in the per-atom-topology + NPT interaction
that I could not isolate in this session.

**Three paths forward** documented earlier:
- Continue NPT debug from this diagnosis baseline
- Partial revert (force kernels back to flat, scaffolding stays)
- Hard rollback to `d8eb8a9`

Commits `3ecefef`, `2a1c64a`, `5d8c777` landed the switchover on top of
the scaffolding:

- `cpu::Frame::Atoms` is now `cpu::Atom` (alias of new type, member name
  kept as `atoms` for ~1600 call-site continuity)
- `Frame::Bonds/Angles/Dihedrals` are kept as **dead-storage shims** so
  ~30 tests + benchmarks that populate them still build until they're
  individually migrated to the per-atom `atoms.add_bond / add_angle /
  add_dihedral` setter API. **Delete these dead shims** in the test
  fixture cleanup.
- `materialize.cpp` and `lammps_data.cpp` now scatter into per-atom
  topology + `build_special_list()`.
- Bonded force kernels (`bond_harmonic`, `angle_harmonic`,
  `dihedral_multi_harmonic`) iterate per-atom (`for i in [0, nlocal):
  for m in [0, atoms.num_bond[i])`). They fall back to
  `dd->tag_to_local` in multi-rank mode (still maintained by old halo
  path until the new `cpu::comm::Comm` is wired in via install_defaults).
- `pair_lj_cut` reads 1-2/1-3/1-4 exclusion from `atoms.nspecial /
  atoms.special`.
- `partition_atoms` and `AtomMigration` were extended to carry per-atom
  topology through compaction / pack-exchange.
- `Atom::reshuffle_topology_to_tags()` added so the bench's `--ckpt-load-before`
  path can move topology rows to match the Zarr-loaded tag permutation.

**Known bug** (unblocks completion): 8-rank MSU PE benchmark (full
`polyethylene_tension --mode nvt --data msu_pe.data`) shows astronomical
bonded energies at step 0 — E_bond ~ 2.85e8 eV vs the correct ~0 eV
at lattice positions. Single-rank PE tests pass (3 of 4 — the 4th
fails by 1.6% in a 1e-9 tolerance check, FP-order sensitivity).
`test_mpi_pe_msu_pipeline` (in-code construct via `make_pe_box` writing
*flat* topology) passes — but that passes vacuously because flat
topology is dead storage and per-atom topology is never populated,
so bonded forces compute on `num_bond = 0` everywhere.

The bug is likely in one of:
1. **lammps_data add_bond ownership routing**: when n_atoms = 10000,
   the bond owner is `min(id_to_idx[ai], id_to_idx[aj])`. id_to_idx
   was just populated in the Atoms-section pass; check that map_set
   ran for the full row range before the Bonds section.
2. **partition_atoms in-place compaction**: rows iterate ascending,
   sources are sorted ascending — should be safe but verify with a
   print of `num_bond[0..3] / bond_atom[0][0..2]` before and after
   the compaction loop.
3. **AtomMigration pack/unpack**: the new pack uses `kMaxExchangeDoubles`
   fixed-stride. Verify the pack writes exactly that many doubles
   per atom, and the unpack reads the same — adding a marker double
   at the head would help diagnose offset drift.

**Recommended path to root-cause**: add `printf` of `frame.atoms.tag[0..3]`,
`num_bond[0..3]`, `bond_atom(0, 0..2)` after each of (a) lammps_data
load, (b) partition_atoms compaction, (c) first halo full-rebuild,
(d) first bond_harmonic launch. Compare to single-rank values at the
same logical position.

**Debug-print finding (briefly, before reverting the prints)**: at
partition_atoms's BEFORE-compact print for rank 0's owned_indices[0]=746,
the value `tag[746]` is unexpectedly `5223` (not `746` as Frame ctor's
identity default would give). This points to lammps_data — or
somewhere — writing tags ≠ identity before partition runs. Best
next step: instrument `lammps_data.cpp` right after the Atoms section
to print `tag[0..3]`; if those are already non-identity, the bug is
in the Atoms parser. If they're identity, walk forward to find the
mutation.

**Rollback path**: `git revert 5d8c777 2a1c64a 3ecefef` restores
force kernels to flat-topology + tag_to_local (the working pre-switchover
state). The new `cpu::Atom / cpu::Domain / cpu::comm::Comm` scaffolding
(commits `6edff73..741927f`) stays as a clean foundation for the next
attempt.

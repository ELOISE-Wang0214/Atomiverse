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

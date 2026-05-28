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

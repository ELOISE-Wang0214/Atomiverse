# MSU PE 10K-atom multi-rank NPT collapse — unresolved

**Date**: 2026-05-29
**Status**: Open. Reproduces consistently. Distinct from the smaller "chains-per-axis=4" NPT bug (that one is fixed).

## Symptom

8-rank NPT with the real MSU PE deck (`benchmarks/md/reference/polyethylene/msu_pe.data`,
10000 united-atom CH2 sites, real chain topology) collapses at step ~6
(when `++neighbor_generation` is bumped every NPT step in `nose_hoover_npt.cpp:670`)
or at step ~10 (when not bumped, default halo `every_=10` cadence).

Step-by-step (with NPT bump on, starting from fresh `stage1.zarr` checkpoint
captured by current binary):

```
step   T_K       P_bar      E_bond   E_kin
  0    519       -8946       2162    670
  1    543       -8809       2131    702
  2    590       -8502       2069    762
  3    658       -8036       1980    850
  4    744       -7427       1868    961
  5    844       -6697       1736    1091
  6    9935      +20290      1590    12841   ← spike
  7    35306     +93082      1962    45633
 13    1.7e6     +4e6        33648   2.2e6   ← total blowup
 14    7.7e-159  +1.3e9      117982  1e-158  ← velocities zeroed
 15+   0         varying     105000  0       ← atoms frozen, P drifts
```

T and KE go to ~0 (denormal floats) after step 13, atoms freeze at high-PE state,
NPT integrator's per-step box scale clamp (`kEpsDtMaxFraction=0.005`) prevents box
runaway but conserved Hamiltonian keeps drifting.

## What was tried

A/B at 8 rank with fresh `stage1.zarr`:

1. **All forces enabled, NPT bump on** → crash at step 6 (first migration step)
2. **All forces enabled, NPT bump off** (every_=10) → crash at step 10 (halo full_rebuild)
3. **bonded forces disabled (LJ only), NPT bump off** → crash at step 10 (same)
4. **all forces disabled** → driver refuses to build (Capability::Force required)

→ The crash trigger is at the **halo full_rebuild step**, regardless of whether
bonded force kernels are active. Not a bonded-force-kernel bug.

The shared element across crash steps:
- Migration runs (clears `tag_to_local`, compacts owners, unpacks new owners,
  bumps `neighbor_generation`)
- HaloExchangeForward sees `generation_changed=true`, does `full_rebuild` —
  scans `[0, n_local)` for face send-list, cascades through 3 axes,
  rebuilds `tag_to_local` for `[0, n_local + n_ghost)`

Both pack/unpack of per-atom topology (`comm::pack_exchange` /
`unpack_exchange` in `src/cpu/comm/atom_pack.cpp`) was audited and looks
symmetric. `move_row` in `atom_migration.cpp:127` was audited and copies
all SoA + per-atom topology fields correctly.

## Why the 1152-atom config works

`--chains-per-axis 4` (1152 atoms, 8 chains × 144 atoms = 18 atoms per chain × 8)
deliberately places each chain inside a single rank's slab. Atoms drift only within
their slab; no migration ever fires. HaloExchangeForward full_rebuild runs every
step (NPT bump) but with **stable ghost set composition** (no migration → ghost slots
stay assigned to the same atoms across rebuilds).

The MSU PE deck has real ~Å-bonded chains that span multiple slabs at 8-rank
(2×2×2 brick with 80Å box → 40Å slabs). Chain motion crosses slab boundaries
in early NPT, triggering migrations. The full_rebuild + migration interaction
is where the bug lives.

## Hypotheses for next round

1. **send_list rebuild uses raw `q[i]` without PBC consideration**
   (`halo_exchange.cpp:274` — `if (q[i] - lo < cutoff)`). NPT integrator does
   `wrap_orthorhombic(f, n)` for owners, so theoretically `q[i] ∈ [0, L)` for
   all owners pre-rebuild. But if migration appends a newly-arrived owner
   whose pre-migration position was wrapped on the SENDER side and the wrap
   was applied to the sender's coordinate system, the post-unpack position
   on the receiver might be outside `[0, L)`. Add an assertion for
   `wrap_orthorhombic` consistency post-migration before halo rebuild.

2. **ghost slot index reshuffle between rebuilds breaks LJ NL CSR**.
   `halo_generation` bump after full_rebuild signals NL to rebuild CSR
   (halo_exchange.cpp:636, neighbor_list.h:144 `halo_stale` check). But
   between when migration writes `dd->tag_to_local` (line 302 of
   atom_migration.cpp) and halo rebuild calls `forward_axis_stage`,
   `tag_to_local` only has owners. If any in-flight component reads
   `tag_to_local` between these two points expecting ghosts to be there,
   garbage results. Add invariant: only NL and force kernels read
   `tag_to_local`, both after halo rebuild.

3. **per-atom topology partner tag falls outside halo ghost set at corner/edge
   slabs**. With 2×2×2 brick, corner ranks have 7 cartesian neighbors via
   cascading 6-face halo. If a bond crosses a corner (e.g., atom in rank 0
   bonded to atom in rank 7), the halo cascade for axis x → y → z must
   eventually ship rank 7's atom to rank 0. Verify `forward_axis_stage`'s
   3-axis cascade correctly reaches corner neighbors. Session memory
   already documents one cpu_mpi_halo_edge_corner_bug — this may be related
   or a sibling.

4. **Repro on a single rank**: change `ATV_FORCE_DIMS=8,1,1` (1D
   decomposition) and rerun. If 1D doesn't crash but 2×2×2 does → corner/edge
   halo bug. If both crash → migration-pack bug. If neither → ghost-slot
   reshuffle in 3D cascade.

## What to instrument next round

- In `halo_exchange.cpp` `rebuild_send_lists_for_axis` (line 261), add
  per-rank atomic counter for `q[i] < 0 || q[i] >= box.L[dim]` cases — should
  be zero post-`wrap_orthorhombic`.
- In `atom_migration.cpp` after Alltoallv unpack, dump first 3 atoms'
  `(x,y,z, tag, num_bond, bond_atom[0..nb])` per rank to stderr at
  step==6 (with NPT bump).
- In `neighbor_list.h::launch`, log when `halo_stale` triggers rebuild
  vs disp_rebuild — confirm halo bumps generation correctly.

## Reference

Same bug as documented in the previous session's compaction summary:
"Unsolved: NPT × multi-rank × per-atom topology collapse — T=0, E_kin=0,
conserved=1.85e151 at stage-2 step 40". My `++neighbor_generation` fix at
`nose_hoover_npt.cpp:670` shifted the crash from step 40 to step ~6-10 by
forcing per-step halo rebuilds, but did **not** address the underlying
halo full_rebuild + migration interaction bug. The 1152-atom config never
exposed it because no migrations fire.

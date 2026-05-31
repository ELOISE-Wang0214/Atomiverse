# CPU MPI halo: face-only exchange misses edge/corner neighbors

**Date:** 2026-05-27
**Status:** Root cause identified, fix not yet implemented
**Bench symptom:** MSU PE 10k atoms at 8 ranks (2×2×2 PBC) — stage 2 NPT T drifts from 500→570 K by step ~400, runs away to 1000+ K by step ~500, catastrophic explosion around step ~700. Each previously-committed fix (ghost-rows, partition, min-image multi-wrap) moved the crash later but didn't eliminate it.

## Root cause

`src/cpu/md/mpi/halo_exchange.cpp::HaloExchangeForward::launch` exchanges with only **6 face neighbors** (`±x, ±y, ±z`) in a single MPI round:

```cpp
for (int f = 0; f < 6; ++f) {
    // build send_list[f] for face f, send to dd->face_neighbors[f]
}
```

In a true 3D Cartesian decomposition with PBC (e.g. 2×2×2), each rank has up to 26 neighbors:
- 6 face neighbors
- 12 edge neighbors (differ in 2 axes)
- 8 corner neighbors (differ in 3 axes)

Many collapse via PBC for small grids — but they're still distinct ranks that need to receive ghosts.

For PE-style polymer chains, a bond `(i, j)` whose owner-tag atom is on rank A and partner is on rank A's edge or corner neighbor B is silently dropped on A:
- `detail::refresh_resolved_locals<K>` looks up both atom tags in `dd->tag_to_local`.
- Partner's tag is missing (not in halo on A) → `continue` (skip).
- A is the bond's owner per `min(tag_i, tag_j)` rule, but A cannot see the partner.
- Per Newton's 3rd law force pairs, A's owner atom doesn't get this bond's force; partner's owner rank also doesn't compute it (it's not the owner per the rule).

## Confirmed empirically

`[partition]` trace + `[tuple_skip_bond]` diagnostic at stage-2 step 0 (load from Zarr checkpoint at end of stage 1):

```
[partition] tag=272 pos=(65.580,41.187,40.077) coord=(1,1,1) → rank 7
[partition] tag=273 pos=(67.431,38.781,37.217) coord=(1,0,0) → rank 4
[partition] tag=274 pos=(66.865,38.219,38.541) coord=(1,0,0) → rank 4

[tuple_skip_bond] rank=7 step=0 tags=(272,273) missing_side=j
  owner_pos=(65.580,41.190,40.079) box=(80.059,...)
```

Rank 4 = `(1,0,0)` differs from rank 7 = `(1,1,1)` in **y AND z** → edge neighbor, not face. 273 is never sent as ghost to rank 7. Bond force on atom 272 missing every step.

Same pattern on rank 1 (tag 9596 missing partner), rank 5 (tag 749), etc. — many bonds, angles, dihedrals dropped.

## Fix options

1. **LAMMPS-style multi-stage forward exchange** (proper fix, matches `comm_brick::forward_comm`):
   - 3 sequential axis stages (or 6 face stages).
   - Each stage's pack includes owners **plus** ghosts received in earlier stages.
   - After all stages, ghosts reach all 26 neighbors via cascade.
2. **Direct 26-way exchange**:
   - Add 12 edge + 8 corner faces to `dd->send_lists`.
   - More straightforward but more bookkeeping.

Effort estimate: ~1-2 days for #1 (refactor `HaloExchangeForward::launch` + `rebuild_send_lists` + buffer management).

## Workaround (no code change)

1D decomposition is safe — no edge/corner. For 8 ranks: `dims = {8, 1, 1}`. But:
- Doesn't scale to 32/64/128 ranks with reasonable slab depth.
- `MPI_Dims_create` currently picks 3D automatically.
- Would need an override env var or CLI flag.

## Why hidden until now

- `test_mpi_nose_hoover_npt_multirank` (8 ranks, 80 atoms): too few atoms for chain bonds to cross 3D corners. Passes.
- LAMMPS reference comparison: LAMMPS *does* multi-stage halo, so it works. Atomiverse looks correct at small scale.
- All previously-committed multi-rank fixes (ghost-rows, partition, min-image) were *necessary* but didn't expose this layer.

## Diagnostic tools (in git history, can be reapplied)

- `[partition]` trace: `src/cpu/md/mpi/domain.cpp::partition_atoms` — log a specific tag's destination rank.
- `[tuple_skip]`/`[tuple_skip_bond]`: `src/cpu/potentials/topology_filter.h::refresh_resolved_locals` — report any bonded tuple whose partners aren't in `tag_to_local`.

## Reproducer

```bash
# Stage 1 dump (one-time, ~17 min):
sbatch /tmp/atv_pe_dump1.sbatch
# Stage 2 from checkpoint with print-every=1 (~3 min iteration):
sbatch /tmp/atv_pe_load2.sbatch
# Crash at step 704 (T>2000 detector) — same trajectory regardless of code state.
```

## Already-fixed companion bugs (all needed in production but not sufficient)

- `00d83f5` ghost-rows scaling in NPT/Deform/Berendsen
- `cc21bdf` skip `partition_atoms` on `adopt_dd` (atom count preservation)
- `66eb9cb` multi-wrap `Box::min_image` (was the catastrophic step-194 bug)
- `9b1f70a` multi-wrap min-image in `cell_list` + `pair_lj_cut`
- `2eeef52` Zarr checkpoint dump+load + per-step mid-stage hook (fast debug iteration)
- molrs `4ff6121`: `molrec_read_zarr_first_frame` FFI (needed by atomiverse `ZarrReader`)

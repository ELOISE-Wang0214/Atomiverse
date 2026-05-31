# MSU PE 10K-atom multi-rank NPT collapse — RESOLVED

**Date opened**: 2026-05-29
**Date resolved**: 2026-05-29
**Status**: FIXED. Root cause was a stale special-bonds exclusion cache in
`cpu::PairLJCut`, NOT a halo / migration topology bug.

## Root cause

`cpu::PairLJCut` (`src/cpu/potentials/pair_lj_cut.cpp`) built its
tag-keyed special-bonds exclusion map `excluded_by_tag_` **once, in
`build()`**, from the atoms that rank owned *at build time*. The map is
keyed by global tag (migration-stable in principle), but it was never
refreshed afterwards.

When `AtomMigration` moved an atom into a rank whose build-time snapshot
never contained that atom's tag (e.g. an atom that started >3 bonds away
from the slab boundary and later drifted across), `rebuild_filtered_list`
looked up `excluded_by_tag_[tag_i]`, found nothing, and applied **no
exclusion**. The migrated atom's 1-2/1-3/1-4 partners (a ~1.53 Å bonded
pair) were then handed to the LJ kernel → `(σ/r)¹²` with `σ/r ≈ 2.6` →
thousands of eV per pair → catastrophic repulsive spike.

This explains every observed symptom:
- **E_pair / P explode while E_bond stays low** at the spike step — it is
  LJ on a bonded pair, not the bond kernel.
- **"LJ-only still crashes"** — the exclusion uses topology, independent
  of whether the bond *force* kernel is active.
- **Config-dependent onset** (step ~6 with one checkpoint, step ~68 with
  another) — it triggers the first time a *non-boundary* atom migrates
  into a rank whose build map lacks it AND a bonded partner is within
  cutoff there.
- **1152-atom config never reproduced** — each chain stays inside one
  rank's slab, so no atom ever migrates; the build-time map is never
  stale.

The per-atom `nspecial`/`special` lists themselves were always correct —
they travel with the atom via `pack_exchange`/`unpack_exchange`/
`move_row`. Only the cached `PairLJCut` map went stale.

## Fix

`src/cpu/potentials/pair_lj_cut.cpp` + `.h`:
- Extracted the map construction into
  `refresh_exclusion_map(frame, n_owners)`, which rebuilds
  `excluded_by_tag_` from the **live** per-atom topology of the current
  owner rows.
- Call it at the top of `rebuild_filtered_list`, which is already gated on
  `pairs.version` and therefore re-runs on every neighbor-list rebuild —
  i.e. after every migration (the halo full-rebuild bumps
  `halo_generation` → NL `halo_stale` rebuild → `pairs.version++`).
- Gate the filtered path on a new `any_excl_` flag instead of
  `!excluded_by_tag_.empty()`, so the refresh stays active even on a step
  where this rank momentarily owns no excluded pair.
- `build()` still seeds the map once (for single-`compute()` callers like
  `Evaluate` that never reneighbor).

## Validation

- **Real deck**: 8-rank stage-2 NPT from `stage1.zarr` previously
  collapsed at step 68 (T 508→11879→0, conserved 1151→1.5e111). After the
  fix it runs 500 stable steps (T 500–516 K, conserved 1150.7→1145.3,
  box growing smoothly toward P=0).
- **Regression test**:
  `tests/mpi/test_nose_hoover_npt_multirank.cpp ::
  MigrationKeepsSpecialBondsExclusion` (8-rank). Chains laid entirely in
  the x>box/2 half with a −x drift march across the slab boundary under
  NVE; asserts `|E_pair|` stays bounded. With the fix reverted it FAILS
  at step 60 (E_pair = 8.5e9 eV); with the fix it holds at −0.39 eV
  through 150 steps.

## Follow-ups (not blocking)

- `cuda::PairLJCut` may have the same build-once exclusion pattern; the
  CUDA path does not currently exercise CPU `AtomMigration`, so it was not
  reproduced here. Audit if CUDA multi-rank migration is added.
- The NL ghost-pair dedup in `cell_list.cpp:413` uses strict `<`
  (`tag_ptr[j] < tag_ptr[i]`); it would not exclude a genuine
  same-tag duplicate ghost. No duplicate is produced today (verified the
  collapse was the exclusion bug, not a duplicate), but `<=` would be a
  cheap defense-in-depth if a future halo change can ever duplicate a tag.

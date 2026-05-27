# MSU stage-2 NPT — extended Hamiltonian drift at every halo full rebuild

**Date:** 2026-05-27
**Status:** Root cause area localized, mechanism not yet pinned, two fix attempts failed.
**Bench symptom:** MSU PE 10k atoms at 8 ranks, fresh stage1.zarr loaded — T slowly drifts 500 → 900 K over ~900 steps, then catastrophic ε̇ runaway → T → 2.5×10⁶ → NaN cascade by step ~1100.

## Smoking-gun signal

After landing the multi-stage halo cascade fix (commit `3a9f1a0`), the residual T-runaway is **not** from missing corner ghosts — it's a slow accumulation of energy that overwhelms the thermostat.

Adding `f.scalars.conserved` (MTK extended Hamiltonian, the quantity NoseHoover-NPT preserves to machine precision in LAMMPS) to the bench's `thermo.csv` showed:

```
step  cons      delta
0     1150.67   —
1     1150.60   −0.07     ← FP noise
2..9  1150.61   ±0.01     ← FP noise
10    1147.18   −3.43     ← JUMP at halo full rebuild (every_=10)
11..19 1147.14  ±0.01
20    1148.25   +1.11     ← JUMP
30    1150.34   +2.05     ← JUMP
40    1151.46   +1.13     ← JUMP
50    1153.86   +2.43     ← JUMP
```

`conserved` jumps by 1-3 eV at every halo full-rebuild step (`every_=10`); flat to FP noise between rebuilds. Cumulative drift → 5000 eV by step 1000 → ε̇ runaway → catastrophic blowup.

## Confirmed pre-existing

`77a72ce`'s commit message already references the same runaway under the buggy 6-face halo: "stage-2 T-runaway (now ~step 1500 in 1D, was step 700 in 3D)". The halo cascade fix moves the crash from step ~700 to step ~1000 — same drift mechanism, different timing.

## Fixes tried, both **failed**

### Attempt 1: halo cutoff += NL skin (1.3× pair_cutoff)

Hypothesis: ghosts flicker in/out of the cutoff zone as atoms drift across the 10.5 Å boundary; LAMMPS exchanges ghosts within `cutoff + skin` so the band absorbs drift smoothly.

Implementation: `MD::install_infra()` set halo cutoff to `1.3 × pair_evaluator_cutoff` (matches NL default auto-skin `0.3 × cutoff`).

Result: **No effect on the drift pattern at MSU scale** (same -3.43, +1.11, +2.05 jumps). **Capacity overflow** on small-box tests (`test_pe_msu_pipeline` 144 atoms in 48 Å box, slab 24 Å < 2 × 13.65 = 27.3 → atoms enter both ±d send_lists → ghost duplication → frame.atoms.size() overflow → MPI_Abort).

Reverted.

### Attempt 2: bump `dd->neighbor_generation` at halo full rebuild

Hypothesis: halo full rebuild can reshuffle ghost slot indices (atoms entering/leaving cutoff, or just re-binning across axes); NL's CSR pair list stores `j` as raw slot index; without a generation bump the CSR keeps pointing at slots whose *atoms* have changed → small per-rebuild error in LJ + bonded virial.

Implementation: `++dd->neighbor_generation` at the end of halo full rebuild path.

Result: **CATASTROPHIC** — `conserved` jumps to 1×10²⁸ at step 1 (immediate), T=310000 K at step 1. The bump apparently interacts disastrously with `atom_migration`'s first_launch path on Zarr-loaded state. Reverted.

## Other hypotheses (untested)

- Bonded virial accumulating across rebuild boundaries (e.g., `bond_harmonic`'s `owned` cache returning a stale list after a halo rebuild reshuffled slot indices)
- NoseHooverNPT MTK Trotter symmetry violation specific to the Atomiverse phase split (P1-P6 in `initial_integrate`, P7 in `final_integrate`, P8-P11 in `end_of_step`) vs LAMMPS's all-in-final_integrate
- ε̇ update at P3 uses *previous step's* stress (Atomiverse intentional — `ReduceEnergy` runs in end_of_step). LAMMPS recomputes pressure fresh at initial_integrate, breaking that 1-step lag

## How to reproduce + iterate

1. Load `.debug/halo_corner/stage1.zarr` (fresh, fixed-binary 10k-atom MSU end-of-stage-1 frame).
2. `srun -n 8 build-main/benchmarks/polyethylene_tension --data <msu_pe.data> --mode npt --ckpt-load-before 2:stage1.zarr --scale 0.04 --print-every 1 --out-dir <out>`
3. Inspect `<out>/thermo.csv`'s `conserved_eV` column. Expect jumps at every multiple of 10.

`.debug/halo_corner/verify_stage2.sbatch` is the canonical reproduction script.

## Diagnostic infrastructure that landed

- `benchmarks/md/polyethylene_tension.cpp`: thermo CSV gains a `conserved_eV` last column (frame.scalars.conserved is written by NoseHooverNPT::launch at end_of_step). Downstream tools (`compare_msu.py`, `plot_msu_compare.py`) are column-name-keyed so they're unaffected.
- `include/atomiverse/cpu/md/cell_list.h`: `CellList::skin()` accessor.
- `include/atomiverse/cpu/md/neighbor_list.h`: `NeighborList<Strategy>::skin()` accessor.
- `tests/mpi/test_nose_hoover_npt_multirank.cpp`: `BondedKineticNonZero` + `BondsOnlyKineticNonZero` now accept 1-rank as well as 8-rank (`world_size != 8 && world_size != 1 → SKIP`), so the same test can localize whether a particular failure is multi-rank-specific.

## Next-session pickup

The conserved-jump pattern at halo full-rebuild boundaries is the smoking gun. Two productive next steps:

1. **Halo-only counter:** instead of bumping `dd->neighbor_generation` (which triggers `atom_migration`), add a separate `dd->halo_generation` counter that only NL subscribes to. Avoids the gen-counter coupling that broke attempt 2.
2. **Per-phase H trace:** log `extended_hamiltonian` after each of P1-P11 inside NoseHooverNPT for one rebuild-boundary step (e.g., step 10), and the same for the previous non-rebuild step (step 9). The phase where H jumps localizes the leak.

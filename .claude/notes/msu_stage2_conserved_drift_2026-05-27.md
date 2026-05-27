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

## Update (later same day) — partial fixes landed in 7b9263d

Two follow-on attempts:

3. **Halo-only `halo_generation` counter** — implemented in 7b9263d. Decoupled from `neighbor_generation` so the NL CSR rebuilds at every halo full rebuild without triggering `atom_migration`'s first_launch cascade. **Neutral effect** on the conserved trajectory (drift pattern identical) — confirms the bug is NOT in NL stale-CSR.
4. **10 % halo skin (1.1× pair_cutoff)** — also landed in 7b9263d. 1.3× = NL auto-skin overflows the 144-atom-in-48 Å test slabs (`2*halo > slab`); 1.1× is the largest factor that still keeps that test under capacity. **Partial benefit**: the eventual catastrophe at step ~1100 is softened (T = 99k K instead of 807k K), but the per-rebuild +1-3 eV drift up to step 1000 is identical to baseline. So skin extension helps at high T (when atoms move fast enough to actually drift past the bare-cutoff halo edge) but doesn't address the underlying mechanism.

## Mechanism is STILL unclear

The boundary-step thermo breakdown nails the leak to **bonded forces**, not LJ:

```
step 9:  E_pair=-347.056  E_bond=230.901  E_angle=227.355  E_dih=397.511  E_kin=641.882  cons=1150.61
step 10: E_pair=-347.055  E_bond=228.456  E_angle=226.738  E_dih=395.648  E_kin=643.376  cons=1147.18
                                ↑Δ=-2.45        ↑Δ=-0.62        ↑Δ=-1.86       ↑Δ=+1.49        ↑Δ=-3.43
```

ΔE_pair is FP noise (-0.001), ΔE_bond/angle/dih are 5 eV combined drop, ΔE_kin only catches 1.5 of it → 3.43 eV escapes the atom subsystem. NoseHoover-NPT's thermostat / barostat chains do not absorb it.

Hypothesis "atoms drift past halo cutoff → cross-rank bonded tuple skipped" doesn't survive the 1.1× skin null result. The bonded mass loss is real but its mechanism is not "ghost set shrinking".

## Next-session pickup (still open)

- **Per-phase Hamiltonian trace.** Log `extended_hamiltonian` inside NoseHooverNPT::launch after each of P1-P11 for the step-9→step-10 boundary and step-19→step-20 boundary. The phase that introduces the jump localizes the leak. (Probably need a debug-only flag — the conditional reads + logging adds ~100 ns/step which is fine for a 2000-step verify.)
- **Bond-by-bond audit.** Dump every owner-resolved bond's `r` at step 9 and step 10, identify which bonds changed by more than typical thermal drift. If a specific subset is mis-handled, that points at a specific topology corner case.
- **Cross-check vs LAMMPS at single-rank.** Single-rank Atomiverse currently can't reload stage1.zarr (cross-rank-count tag/index bug noted earlier), but a fresh single-rank stage-1 + stage-2 run with `econserve` thermo would tell us whether the drift exists at single-rank too — i.e., whether this is an MPI-specific bug or a NoseHooverNPT integrator bug that's always present and only exposed by the larger 8-rank energy scale.

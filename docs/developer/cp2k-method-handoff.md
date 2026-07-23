# CP2K capability loop — handoff

Last harness pass: assign wave 1 (2026-07-23) — 16 items → `in_progress`.

## Wave summary

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) |
| ---- | ------ | ------ | ---------- | --------------------------------------------- |
| gpw-water-pbe-mt | gpw | [~] wave 1 | — | — |
| gpw-forces-water | gpw | [~] wave 1 | — | — |
| gpw-uks-o2-gamma | gpw | [~] wave 1 | — | — |
| gpw-stress-water | gpw | [~] wave 1 (after energy/forces) | — | — |
| gpw-pbe-d3bj-ice | gpw | [~] wave 1 | — | — |
| pw-lda-si | pw | [~] wave 1 | — | — |
| pw-pbe-gga-si | pw | [~] wave 1 | — | — |
| pw-smearing-al | pw | [~] wave 1 | — | — |
| pw-spin-polarized-o2 | pw | [~] wave 1 | — | — |
| pw-analytic-stress-si | pw | [~] wave 1 (after pw-lda-si) | — | — |
| gpw-bomd-water | dynamics | [~] wave 1 (after forces) | — | — |
| pw-bomd-h | dynamics | [~] wave 1 | — | — |
| fd-phonon-supercell | response | [~] wave 1 | — | — |
| ot-minimizer-water | gpw | [~] wave 1 | — | — |
| gpw-dftu-nio | gpw | [~] wave 1 | — | — |
| bands-dos-si | pw | [~] wave 1 | — | — |

## Oracle

Offline CP2K pins live in each regression `REFERENCE.md`.
Pin path: `tools/cp2k_oracle/` (see README). **Never** call live CP2K from ctest.

## Next assign candidates (pending head)

1. `mpi-gpw-scaling-4rank` (industrial / `perf_exempt` OK)
2. `berry-polarization` (after k-path / bands)
3. `wannier-centers` (after `bands-dos-si`)

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

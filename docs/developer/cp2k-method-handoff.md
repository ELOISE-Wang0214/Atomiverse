# CP2K capability loop — handoff

Last harness pass: **2026-07-24 wave 2** on `dev`
(`chore(harness): cp2k wf wave 2`).
Prior: **2026-07-23 wave 1** + assign wave 2 (staff 10).

**Workflow order (binding):** Recover incomplete first → Census →
Assign (`incomplete ≫ pending`) → Worktrees (**spec-first**) → Merge → Harness →
Tutorial. See `docs/developer/cp2k-method-loop.md` and
`.grok/workflows/cp2k-method-loop.rhai`.

## Wave 2 summary (harness 2026-07-24)

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) |
| ---- | ------ | ------ | ---------- | --------------------------------------------- |
| `gpw-stress-water` | gpw | **merged → done** | `regressions/gpw/gpw-stress-water` | **4.75** |
| `gpw-bomd-water` | dynamics | **merged → done** | `regressions/gpw/gpw-bomd-water` | **perf_exempt** |
| `pw-analytic-stress-si` | pw | **merged → done** | `regressions/pw/pw-analytic-stress-si` | **perf_exempt** |
| `fd-phonon-supercell` | response | **merged → done** | `regressions/pw/fd-phonon-supercell` | **perf_exempt** |
| `gpw-dftu-nio` | gpw | **merged → done** | `regressions/gpw/gpw-dftu-nio` | **0.21** |
| `bands-dos-si` | pw | **merged → done** | `regressions/pw/bands-dos-si` | **perf_exempt** |
| `pw-uks-kmesh-fe` | pw | **merged → done** | `regressions/pw/pw-uks-kmesh-fe` | **perf_exempt** |
| `gpw-pbesol-si` | gpw | **merged → done** | `regressions/gpw/gpw-pbesol-si` | **perf_exempt** |
| `gpw-d4-water` | gpw | **merged → done** | `regressions/gpw/gpw-d4-water` | **perf_exempt** |
| `gpw-kpoints-si` | gpw | **incomplete** | — | — (assign leftover; no merge) |

- **Blocked this wave:** none new (hybrid/GAPW/PAW/CUDA-GPW/meta-GGA stand)
- **Merge-failed this wave:** none
- **Incomplete remaining (1):** `gpw-kpoints-si` — resume first next assign
- **Queue after harness:** pending = 11; in_progress = 1; done = 19
  (wave 1 + wave 2); blocked = 5
- **Roadmap:** CPU ✓ for **GPW stress (strain-FD)**, **GPW BOMD**, **PW
  Nielsen–Martin analytic stress**, **FD supercell phonons**, **periodic
  DFT+U**, **bands/DOS**, **k-point UKS PW**, **PBEsol**, **periodic D4**.
  No ⚡ (ratios / `perf_exempt` only).
- **Specs closed:** all 9 wave-2 merged slugs (acceptance met on `dev`);
  do not re-open without superseding intent.

**Resume first:** `gpw-kpoints-si` (only open incomplete). Then pending head
after deps cleared by wave 2.

## Wave 1 summary (historical)

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) |
| ---- | ------ | ------ | ---------- | --------------------------------------------- |
| `gpw-water-pbe-mt` | gpw | **merged → done** | `regressions/gpw/gpw-water-pbe-mt` | **perf_exempt** |
| `gpw-forces-water` | gpw | **merged → done** | `regressions/gpw/gpw-forces-water` | **3.56** |
| `gpw-uks-o2-gamma` | gpw | **merged → done** | `regressions/gpw/gpw-uks-o2-gamma` | **0.98** |
| `gpw-pbe-d3bj-ice` | gpw | **merged → done** | `regressions/gpw/gpw-pbe-d3bj-ice` | **perf_exempt** |
| `pw-lda-si` | pw | **merged → done** | `regressions/pw/pw-lda-si` | **perf_exempt** |
| `pw-pbe-gga-si` | pw | **merged → done** | `regressions/pw/pw-pbe-gga-si` | **perf_exempt** |
| `pw-smearing-al` | pw | **merged → done** | `regressions/pw/pw-smearing-al` | **0.008** |
| `pw-spin-polarized-o2` | pw | **merged → done** | `regressions/pw/pw-spin-polarized-o2` | **perf_exempt** |
| `pw-bomd-h` | dynamics | **merged → done** | `regressions/pw/pw-bomd-h` | **perf_exempt** |
| `ot-minimizer-water` | gpw | **merged → done** | `regressions/gpw/ot-minimizer-water` | **~9.0** |

## Next assign candidates (pending head, after incomplete clear)

1. `gpw-kpoints-si` (**incomplete** — must clear first)
2. `gpw-variable-cell-ice` (**dep:** `gpw-stress-water` ✓)
3. `dfpt-phonon-si` (**dep:** `fd-phonon-supercell` ✓) / `xl-bomd-water` (**dep:** BOMD ✓)
4. `berry-polarization` / `wannier-centers` (**dep:** `bands-dos-si` ✓)
5. `periodic-mp2-si` / MPI scaling industrial

## Oracle

Offline CP2K pins live in each regression `REFERENCE.md`.
Pin path: `tools/cp2k_oracle/` (see README). **Never** call live CP2K from ctest.

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

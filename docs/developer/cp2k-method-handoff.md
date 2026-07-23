# CP2K capability loop — handoff

Last assign pass: **2026-07-24 wave 3** on `dev`
(`chore(queue): cp2k assign wave 3`).
Prior harness: **2026-07-24 wave 2** (`chore(harness): cp2k wf wave 2`).

**Workflow order (binding):** Recover incomplete first → Census →
Assign (`incomplete ≫ pending`) → Worktrees (**spec-first**) → Merge → Harness →
Tutorial. See `docs/developer/cp2k-method-loop.md` and
`.grok/workflows/cp2k-method-loop.rhai`.

## Assign wave 3 (2026-07-24) — staff 10

| Priority | Slug | Chapter | is_resume | Notes |
| -------- | ---- | ------- | --------- | ----- |
| 1 | `gpw-kpoints-si` | gpw | **true** | only incomplete; no open spec — write-first |
| 2 | `gpw-variable-cell-ice` | gpw | false | dep stress ✓ |
| 3 | `dfpt-phonon-si` | response | false | dep FD phonons ✓ |
| 4 | `xl-bomd-water` | dynamics | false | dep GPW BOMD ✓ |
| 5 | `berry-polarization` | response | false | dep bands ✓ |
| 6 | `wannier-centers` | response | false | dep bands ✓ |
| 7 | `cpmd-water` | dynamics | false | after solid BOMD ✓ |
| 8 | `periodic-mp2-si` | gpw | false | solid Γ GPW energy ✓ |
| 9 | `mpi-gpw-scaling-4rank` | infra | false | industrial; `perf_exempt` OK |
| 10 | `mpi-pw-scaling-4rank` | infra | false | k-point RR; industrial |

- **Held pending (2):** `periodic-rpa-si` (after MP2), `periodic-gw-si` (staff full)
- **Queue after assign:** pending = 2; in_progress = 10; done = 19; blocked = 5
- **Workers:** each NEW slug writes `.claude/specs/<slug>/` first; resume
  `gpw-kpoints-si` also needs a live spec (INDEX: no open CP2K spec dir)

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
| `gpw-kpoints-si` | gpw | **→ wave 3 resume** | — | incomplete leftover from assign wave 2 |

- **Blocked this wave:** none new (hybrid/GAPW/PAW/CUDA-GPW/meta-GGA stand)
- **Merge-failed this wave:** none
- **Specs closed:** all 9 wave-2 merged slugs (acceptance met on `dev`);
  do not re-open without superseding intent.

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

## Next after wave 3 (held pending)

1. `periodic-rpa-si` (**after** `periodic-mp2-si` foothold)
2. `periodic-gw-si` (deps clear; staff-held)

## Oracle

Offline CP2K pins live in each regression `REFERENCE.md`.
Pin path: `tools/cp2k_oracle/` (see README). **Never** call live CP2K from ctest.

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

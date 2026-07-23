# CP2K capability loop — handoff

Last harness pass: **2026-07-23 wave 1** on `dev`
(`chore(harness): cp2k wf wave 1`).
Last assign: **2026-07-23 wave 2** — 6 resumes + 4 new pending → `## in_progress`
(staff 10; see queue).

**Workflow order (binding):** Recover incomplete first → Census →
Assign (`incomplete ≫ pending`) → Worktrees (**spec-first**) → Merge → Harness →
Tutorial. See `docs/developer/cp2k-method-loop.md` and
`.grok/workflows/cp2k-method-loop.rhai`.

## Wave 1 summary

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

- **Blocked this wave:** none
- **Merge-failed this wave:** none
- **Assign wave 2 (2026-07-23):** staff 10 = **6 incomplete resumes** + **4 new**
  pending promoted. No open CP2K/GPW/PW specs in `.claude/specs/` (Live is
  electrochem-only). No unmerged `wf/<slug>` for these 10 slugs.
- **In progress (10):** resumes — `gpw-stress-water`, `gpw-bomd-water`,
  `pw-analytic-stress-si`, `fd-phonon-supercell`, `gpw-dftu-nio`, `bands-dos-si`;
  new — `gpw-kpoints-si`, `pw-uks-kmesh-fe`, `gpw-pbesol-si`, `gpw-d4-water`.
- **Queue after assign:** pending = 11; in_progress = 10; done = 10 (wave 1);
  blocked = 5 (hybrid/GAPW/PAW/CUDA-GPW/meta-GGA)
- **Roadmap:** CPU ✓ for **Γ-UKS GPW**, **PW true GGA PBE**, **PW collinear
  spin**, **metallic smearing**, **GPW D3(BJ)**; OT product switch + water pin;
  PW BOMD foothold already ✓ on MD roadmap. No ⚡ (no standing CP2K perf suite
  in ctest — ratios / `perf_exempt` only).
- **Specs closed:** all 10 wave-1 slugs (acceptance met on `dev`); wave-2 workers
  write `.claude/specs/<slug>/` on `wf/<slug>` **before** code (spec-first).

**Resume first (deps ✓):** `gpw-stress-water`, `gpw-bomd-water`,
`pw-analytic-stress-si`, then `fd-phonon-supercell`, `gpw-dftu-nio`,
`bands-dos-si`. Do **not** demote incomplete to pending.

## Next assign candidates (pending head, after wave-2 clear)

1. `gpw-variable-cell-ice` (**dep:** `gpw-stress-water`)
2. `periodic-mp2-si` (Γ GPW energy ✓) / MPI scaling industrial
3. `dfpt-phonon-si` (**dep:** `fd-phonon-supercell`) / `xl-bomd-water` (**dep:** BOMD)
4. `berry-polarization` / `wannier-centers` (**dep:** `bands-dos-si`)

## Oracle

Offline CP2K pins live in each regression `REFERENCE.md`.
Pin path: `tools/cp2k_oracle/` (see README). **Never** call live CP2K from ctest.

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

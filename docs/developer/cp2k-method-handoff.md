# CP2K capability loop — handoff

Last harness pass: **2026-07-23 wave 1** on `dev`
(`chore(harness): cp2k wf wave 1`).

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
- **Incomplete remaining (queue `[~]`, 6):** `gpw-stress-water` (deps ✓ — energy+forces),
  `pw-analytic-stress-si` (deps ✓ — `pw-lda-si`), `gpw-bomd-water` (deps ✓ — forces),
  `fd-phonon-supercell`, `gpw-dftu-nio`, `bands-dos-si`
- **Queue:** pending = 15; in_progress = 6; done = 10 (wave 1); blocked = 5
  (hybrid/GAPW/PAW/CUDA-GPW/meta-GGA)
- **Roadmap:** CPU ✓ for **Γ-UKS GPW**, **PW true GGA PBE**, **PW collinear
  spin**, **metallic smearing**, **GPW D3(BJ)**; OT product switch + water pin;
  PW BOMD foothold already ✓ on MD roadmap. No ⚡ (no standing CP2K perf suite
  in ctest — ratios / `perf_exempt` only).
- **Specs closed:** all 10 wave-1 slugs (acceptance met on `dev`); INDEX Live
  empty for this loop.

**Resume batch (≤6 remaining incomplete, dep-safe first):** `gpw-stress-water`,
`gpw-bomd-water`, `pw-analytic-stress-si`, then `fd-phonon-supercell`,
`gpw-dftu-nio`, `bands-dos-si`. Do **not** demote incomplete to pending.

## Next assign candidates (pending head, after incomplete clear)

1. `gpw-kpoints-si` (deps energy/forces ✓)
2. `pw-uks-kmesh-fe` (deps spin + smearing ✓)
3. `gpw-d4-water` / `gpw-pbesol-si` (GGA/dispersion ladder)
4. `mpi-gpw-scaling-4rank` / `mpi-pw-scaling-4rank` (industrial / `perf_exempt` OK)

## Oracle

Offline CP2K pins live in each regression `REFERENCE.md`.
Pin path: `tools/cp2k_oracle/` (see README). **Never** call live CP2K from ctest.

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

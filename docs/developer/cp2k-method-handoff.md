# CP2K capability loop — handoff

Last recover pass: 2026-07-23 (reconfirm) — inventory of incomplete CP2K/GPW/PW
work on parent `dev`.

**Workflow order (binding, 2026-07-23):** Recover incomplete first → Census →
Assign (`incomplete ≫ pending`) → Worktrees (**spec-first**) → Merge → Harness →
Tutorial. See `docs/developer/cp2k-method-loop.md` and
`.grok/workflows/cp2k-method-loop.rhai`.

**Recover inventory (this shift):** **16** unfinished wave-1 `[~]` rows in
`docs/developer/cp2k-method-queue.md` `## in_progress` (table below). No open
`.claude/specs/<slug>/` for GPW/PW/periodic (INDEX Live empty). No unmerged
`wf/<slug>` for this loop (local `wf/*` are PySCF molecular only). No
half-written `regressions/{gpw,pw}/<queue-slug>/` (product PW footholds live as
`regressions/pw/h_*` only; **no** `regressions/gpw/`). Queue already normalized —
do **not** demote incomplete to pending. **Resume priority:** open spec >
unmerged `wf/` > queue `[~]` > handoff → all 16 are queue `[~]` only; staff up
to **10** dep-safe first (skip stress/BOMD until energy/forces, analytic PW
stress until `pw-lda-si`). Spec is source of truth: workers write
`.claude/specs/<slug>/` **first** on `wf/<slug>`, then implement + offline CP2K
pin (`tools/cp2k_oracle/`).

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

**Resume batch (≤10, dep-safe, staff first):** `gpw-water-pbe-mt`,
`gpw-forces-water`, `gpw-uks-o2-gamma`, `gpw-pbe-d3bj-ice`, `pw-lda-si`,
`pw-pbe-gga-si`, `pw-smearing-al`, `pw-spin-polarized-o2`, `pw-bomd-h`,
`ot-minimizer-water`. **Held (still incomplete, honor deps):**
`gpw-stress-water` (after energy+forces), `pw-analytic-stress-si` (after
`pw-lda-si`), `gpw-bomd-water` (after forces), plus `fd-phonon-supercell`,
`gpw-dftu-nio`, `bands-dos-si` for next staff slots.

## Oracle

Offline CP2K pins live in each regression `REFERENCE.md`.
Pin path: `tools/cp2k_oracle/` (see README). **Never** call live CP2K from ctest.

## Next assign candidates (pending head)

1. `mpi-gpw-scaling-4rank` (industrial / `perf_exempt` OK)
2. `berry-polarization` (after k-path / bands)
3. `wannier-centers` (after `bands-dos-si`)

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

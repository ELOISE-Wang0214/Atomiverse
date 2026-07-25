# CP2K method loop — worktree → merge → dev

Abstract pattern: [`oracle-capability-loop.md`](oracle-capability-loop.md).  
Scaffold skill: `.agents/skills/create-oracle-loop/`.

## Continuous discovery (binding)

**Every wave starts with GapCensus:** re-diff CP2K GPW/PW/periodic surface +
`docs/roadmap/dft.md` vs ATV on `target_branch`, refill `## pending`, update
roadmap. Foothold ≠ capability closed. Stop `queue_empty` only when
`gaps_remaining=0` — **not** after one seed batch.

## Pipeline

```text
parent (dev)
  recover incomplete
       │
  for each wave:
    GapCensus → Assign → Worktrees → Merge? → Clearup → Harness? → Tutorial?
```

Instance = template + CONFIG/`extra_law` only
(`.grok/workflows/cp2k-method-loop.rhai`).

| Role | `isolation_worktree` | Writes |
| ---- | -------------------- | ------ |
| Recover / GapCensus / assign / merge / harness | **false** (parent) | `dev` only |
| Feature worker | **true** | private worktree only until merge |

Workers must report `worktree_path`, `branch` (`wf/<slug>`), `head_sha`.

**Recover incomplete first** (`## in_progress` `[~]`, open GPW/PW specs,
unmerged `wf/<slug>`). Spec is source of truth. Offline CP2K pins only
(`tools/cp2k_oracle/`); never live CP2K in ctest.

**Scope:** GPW / PW / periodic DFT. **Not** molecular post-HF (PySCF loop) or
classical MD electrochemistry (MetalWalls loop).

## Why not edit dev in parallel?

Parallel agents on one tree thrash shared GPW/PW headers and leave a half-broken
`dev`. Isolation = each feature lands as a branch; merge is ordered and
abortable per slug.

## Args

| Arg | Default | Meaning |
| --- | ------- | ------- |
| `staff_count` / `parallelism` | 10 | Concurrent worktrees |
| `max_waves` | 20 | Assignment rounds this shift |
| `target_branch` | `dev` | Merge destination |
| `do_perf` | true | Fair wall vs CP2K when available |
| `mode` | `exhaust` | Until gaps gone / budget |

## Launch

```text
agent_budget: 1024
args: { "mode":"exhaust", "staff_count":10, "max_waves":20, "target_branch":"dev", "do_perf":true }
```

## Queue

`docs/developer/cp2k-method-queue.md`

## Run log

(autonomous harness appends lines here)

- **2026-07-25 harness wave 2 (parent `dev` @ `f5739506`):** 10 merged →
  done: `gpw-kpoint-forces`, `gpw-revpbe-si`, `gpw-kpoints-uks`,
  `gpw-uks-forces`, `gpw-kpoint-stress`, `gpw-mbd-ice`,
  `periodic-mp2-kpoints`, `periodic-rpa-si`, `dfpt-finite-q-si`,
  `gpw-variable-cell-molopt`. Incomplete remaining **0**; pending residual
  **3** (`gpw-rpbe-water`, `gpw-vdw-df-ice`, `periodic-gw-si`); blocked **5**.
  Specs closed (deleted) for all 10. Roadmap: revPBE ✓, free-atom MBD ✓,
  Periodic RPA ✓ (Γ dRPA@KS), finite-q DFPT productized; multi-k MP2 residual
  closed; GPW k-UKS / k-forces / k-stress / Γ-UKS forces / MOLOPT CELL_OPT
  prose. Clearup already ran (no isolation worktrees; no builds). **Queue not
  exhausted** — next wave GapCensus. No push.
- **2026-07-25 GapCensus wave 2 (parent `dev` @ `507bb5f8`):** recover found
  unmerged `wf/gpw-kpoint-forces` (ahead 1) → `## in_progress` `[~]`. Prior
  pending (−k-forces) re-verified true gaps; **+4** follow-ons
  (`dfpt-finite-q-si`, `periodic-mp2-kpoints`, `gpw-kpoint-stress`,
  `gpw-uks-forces`) → pending **12**; in_progress **1**; blocked **5**. No
  ○/✓ table flips; stale DFPT/k-point intro prose refreshed. Next: Assign
  (prefer resume k-forces + surface-disjoint pending).
- **2026-07-25 harness wave 1 (post GapCensus staff; parent `dev` @
  `83a100df`):** 9 merged → done: `gpw-blyp-water`, `gpw-uks-pbe`,
  `pw-uks-fe-smearing`, `pw-variable-cell-si`, `gpw-ts-ice`,
  `phonon-dispersion-si`, `periodic-wannier-si`, `second-gen-cpmd-water`,
  `periodic-ri-mp2-si`. Incomplete remaining **0**; pending residual **9**;
  blocked **5**. Specs closed (deleted) for all 9. Roadmap: TS ✓, Periodic
  MP2 ✓ (Γ RI product; multi-k residual), MLWF ✓, 2nd-gen CP ✓; phonon
  IFC-Fourier productized; GPW BLYP / UKS-PBE / Fe smear / PW CELL_OPT prose.
  Clearup already ran (no isolation worktrees; no builds). **Queue not
  exhausted** — next wave GapCensus. No push.
- **2026-07-25 harness wave 1 (claude, earlier staff batch):** 5 merged into
  `dev` (head `bd39507a`) → done: `gpw-kpoints-si` (MP k-point GPW RKS),
  `gpw-variable-cell-ice` (vs CP2K CELL_OPT), `dfpt-phonon-si` (Γ DFPT, PW),
  `berry-polarization` (KSV, PW), `periodic-mp2-si` (Γ GPW MP2 **foothold** —
  later closed by `periodic-ri-mp2-si`). Incomplete then: `xl-bomd-water` WIP
  (later merged). Merge gate: cp2k domain suite green. Clearup: 71 merged
  `wf/*` pruned. Specs closed for all 5; roadmap Berry ✓, k-points prose;
  MP2 ○ at that time.
- **2026-07-25 bootstrap (parent `dev` @ `5bdae27d`):** GapCensus recover
  clean (incomplete=0, CP2K specs=0, wf/*=0, builds=none). Pending **14**
  re-verified true gaps; done **30**; blocked **5**. No roadmap ○/✓ flips.
  Ready for Assign. Suggested staff: `gpw-blyp-water`, `pw-uks-fe-smearing`,
  `pw-variable-cell-si`, `gpw-ts-ice`, `phonon-dispersion-si`,
  `periodic-ri-mp2-si`, `second-gen-cpmd-water`, `gpw-variable-cell-molopt`
  (hold revPBE/RPBE / MBD/vdW-DF / RPA/GW for surface order).
- **2026-07-25 GapCensus wave 1 (parent `dev` @ `f6838676`):** recover clean.
  Prior pending **14** re-verified; **+4** follow-on slugs (`gpw-uks-pbe`,
  `gpw-kpoints-uks`, `gpw-kpoint-forces`, `periodic-wannier-si`) → pending
  **18**; done **30**; blocked **5**. No table ○/✓ flips; fixed stale
  `docs/roadmap/index.md` GPW blurb. Ready for Assign.

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

- **2026-07-25 harness wave 1 (claude):** 5 merged into `dev` (head
  `bd39507a`) → done: `gpw-kpoints-si` (MP k-point GPW RKS),
  `gpw-variable-cell-ice` (vs CP2K CELL_OPT), `dfpt-phonon-si` (Γ DFPT, PW),
  `berry-polarization` (KSV, PW), `periodic-mp2-si` (Γ GPW MP2 **foothold** —
  capability row stays ○). Incomplete: `xl-bomd-water` WIP checkpoint
  `61d235f4` on `wf/xl-bomd-water` (untested — build never completed). Not
  staffed, still `[~]`: `wannier-centers`, `cpmd-water`,
  `mpi-gpw-scaling-4rank`, `mpi-pw-scaling-4rank`. Merge gate: cp2k domain
  suite green (`ctest -R '^test_cpu_(gpw|pw|scf|xc)_|^test_xc_'`). Clearup:
  71 merged `wf/*` pruned. Specs closed for all 5 merged slugs; roadmap
  synced (Berry ✓, GPW k-points prose, variable-cell productization; MP2 ○).
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

# CP2K method loop — worktree → merge → dev

## Pipeline

```text
parent (dev)
  census + assign
       │
       ├─ worker₁ isolation_worktree  ──┐
       ├─ worker₂ isolation_worktree  ──┼─ pin + impl + verify + commit on wf/<slug>
       └─ workerₙ isolation_worktree  ──┘
       │
  merge agent (parent): git merge each green HEAD → dev
       │
  harness agent (parent): queue + roadmap + notes; commit chore(harness)
```

| Role | `isolation_worktree` | Writes |
| ---- | -------------------- | ------ |
| Census / foreman / merge / harness | **false** (parent) | `dev` only |
| Feature worker | **true** | private worktree only until merge |

Workers must report `worktree_path`, `branch` (`wf/<slug>`), `head_sha`.

**Recover incomplete first** (`## in_progress` `[~]`, open GPW/PW specs,
unmerged `wf/<slug>`). Spec is source of truth. Offline CP2K pins only
(`tools/cp2k_oracle/`); never live CP2K in ctest.

## Why not edit dev in parallel?

Parallel agents on one tree thrash shared GPW/PW headers and leave a half-broken
`dev`. Isolation = each feature lands as a branch; merge is ordered and
abortable per slug.

## Args

| Arg | Default | Meaning |
| --- | ------- | ------- |
| `staff_count` / `parallelism` | 10 | Concurrent worktrees (wave 1 staff batch) |
| `max_waves` | 20 | Assignment rounds this shift |
| `target_branch` | `dev` | Merge destination |
| `do_perf` | false | Optional fair wall vs CP2K |
| `mode` | `exhaust` | Until queue empty / budget |

## Launch

```text
agent_budget: 256
args: { "mode":"exhaust", "staff_count":10, "max_waves":20, "target_branch":"dev", "do_perf":false }
```

`agent_budget` = LLM-turn payroll for the shift, **not** headcount.

## Queue

`docs/developer/cp2k-method-queue.md`

## Run log

(autonomous harness appends lines here)

- 2026-07-23 recover: 16 incomplete `[~]` inventory on `dev`; no open specs; no unmerged CP2K `wf/*`; handoff resume batch of 10 dep-safe
- 2026-07-23 assign wave 1: gpw-water-pbe-mt, gpw-forces-water, gpw-uks-o2-gamma, gpw-pbe-d3bj-ice, pw-lda-si, pw-pbe-gga-si, pw-smearing-al, pw-spin-polarized-o2, pw-bomd-h, ot-minimizer-water → in_progress staff; held: gpw-stress-water, pw-analytic-stress-si, gpw-bomd-water, fd-phonon-supercell, gpw-dftu-nio, bands-dos-si
- 2026-07-23 wave 1 merge: all ten green HEADs → dev (gpw-water-pbe-mt, gpw-forces-water, gpw-uks-o2-gamma, gpw-pbe-d3bj-ice, pw-lda-si, pw-pbe-gga-si, pw-smearing-al, pw-spin-polarized-o2, pw-bomd-h, ot-minimizer-water); conflict repair commits on PW/GPW unions; blocked=[]; merge-failed=[]
- 2026-07-23 harness wave 1: queue merged→done (10); incomplete remaining=6; roadmaps ○→✓ for Γ-UKS GPW, PW true GGA PBE, PW collinear spin, metallic smearing, GPW D3(BJ); specs closed (10); ratios forces=3.56 uks=0.98 smearing=0.008 ot≈9.0 others perf_exempt; queue_remaining=21 (pending 15 + incomplete 6); handoff rewritten
- 2026-07-23 assign wave 2: 6 resumes + 4 new pending → in_progress (staff 10): gpw-stress-water, gpw-bomd-water, pw-analytic-stress-si, fd-phonon-supercell, gpw-dftu-nio, bands-dos-si, gpw-kpoints-si, pw-uks-kmesh-fe, gpw-pbesol-si, gpw-d4-water
- 2026-07-24 wave 2 merge: nine green HEADs on dev (all harness list except gpw-kpoints-si); blocked=[]; merge-failed=[]
- 2026-07-24 harness wave 2: queue 9→done; incomplete remaining=1 (`gpw-kpoints-si`); roadmaps ○→✓ for GPW stress, GPW BOMD, PW analytic stress, FD phonons, periodic DFT+U, bands/DOS, k-UKS PW, PBEsol, periodic D4; specs closed (9); ratios stress=4.75 dftu=0.21 others perf_exempt; queue_remaining=12 (pending 11 + incomplete 1); handoff rewritten
- 2026-07-24 assign wave 3: 1 resume (`gpw-kpoints-si`) + 9 new pending → in_progress staff 10: gpw-kpoints-si, gpw-variable-cell-ice, dfpt-phonon-si, xl-bomd-water, berry-polarization, wannier-centers, cpmd-water, periodic-mp2-si, mpi-gpw-scaling-4rank, mpi-pw-scaling-4rank; held: periodic-rpa-si (after MP2), periodic-gw-si
- 2026-07-24 recover: incomplete_count=10 (all wave-3 `[~]`); open GPW/PW/CP2K specs=0; no unmerged CP2K `wf/*` ahead of `dev` (wave-1/2 local branches fully merged leftovers); no half-written wave-3 regressions; pending remains 2 (RPA after MP2, GW staff-held); blocked 5 unchanged; resume batch = all 10 write-spec-first

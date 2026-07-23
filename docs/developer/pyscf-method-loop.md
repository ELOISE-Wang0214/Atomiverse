# PySCF method loop — worktree → merge → dev

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

## Why not edit dev in parallel?

Parallel agents on one tree thrash shared SCF headers and leave a half-broken `dev`.  
Isolation = each feature lands as a branch; merge is ordered and abortable per slug.

## Args

| Arg | Default | Meaning |
| --- | ------- | ------- |
| `staff_count` / `parallelism` | 5 | Concurrent worktrees |
| `max_waves` | 20 | Assignment rounds this shift |
| `target_branch` | `dev` | Merge destination |
| `do_perf` | false | Optional |
| `mode` | `exhaust` | Until queue empty / budget |

## Launch

```text
agent_budget: 256
args: { "mode":"exhaust", "staff_count":5, "max_waves":20, "target_branch":"dev", "do_perf":false }
```

`agent_budget` = LLM-turn payroll for the shift, **not** headcount.

## Queue

`docs/developer/pyscf-method-queue.md`

## Run log

(autonomous harness appends lines here)

- 2026-07-23 assign wave: frozen-core-mp2-h2o (post-hf), cisd-h2co-sto3g (excited-states), blyp-water-sto3g (dft), ghf-grad-h3-sto3g (hartree-fock), pcm-grad-h2o (multiscale) → in_progress; skip=[]
- 2026-07-23 wave 1 merge: all five green HEADs → dev (blyp-water-sto3g, frozen-core-mp2-h2o, cisd-h2co-sto3g, ghf-grad-h3-sto3g, pcm-grad-h2o); blocked=[]; merge-failed=[]
- 2026-07-23 harness wave 1: queue merged→done (5); in_progress cleared; roadmaps ○→✓ for landed methods; handoff written; queue_remaining=21

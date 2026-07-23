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
- 2026-07-23 assign wave 2: scs-mp2-h2o-sto3g (post-hf), full-lr-tddft-h2-pbe (excited-states), roks-o2-sto3g (dft), static-polariz-h2o (hartree-fock), ief-pcm-h2o (multiscale) → in_progress; skip=[]
- 2026-07-23 wave 2 merge: all five green HEADs → dev (scs-mp2-h2o-sto3g, full-lr-tddft-h2-pbe, roks-o2-sto3g, static-polariz-h2o, ief-pcm-h2o); blocked=[]; merge-failed=[]
- 2026-07-23 harness wave 2: queue merged→done (5); in_progress cleared; roadmaps ○→✓ for SCS-MP2, full LR-TDDFT, ROKS energy, RHF polarizability prose, IEF-PCM already ✓; handoff written; queue_remaining=16
- 2026-07-23 assign wave 3: sos-mp2-h2o-sto3g (post-hf), mp2-grad-h2o-sto3g (post-hf), uks-tda-o2-lsda (excited-states), d3bj-molecular-pbe (dft), uhf-hessian-o2 (hartree-fock) → in_progress; skip=[]
- 2026-07-23 wave 3 merge: all five green HEADs → dev (sos-mp2-h2o-sto3g, mp2-grad-h2o-sto3g, uks-tda-o2-lsda, d3bj-molecular-pbe, uhf-hessian-o2); blocked=[]; merge-failed=[]
- 2026-07-23 harness wave 3: queue merged→done (5); in_progress cleared; roadmaps ○→✓ for RMP2 analytic gradient + RHF analytic Hessian inventory fix; prose for UKS-TDA, UHF Hessian, D3BJ/SOS already ✓; handoff written; queue_remaining=11
- 2026-07-23 assign wave 4: romp2-o2-sto3g (post-hf), roks-grad-o2-sto3g (dft), h2-b3lyp-rpa-tddft (excited-states), spin-flip-tda-o2 (excited-states), rohf-hessian-o2 (hartree-fock) → in_progress; skip=[]
- 2026-07-23 wave 4 merge: all five green HEADs → dev (romp2-o2-sto3g, roks-grad-o2-sto3g, h2-b3lyp-rpa-tddft, spin-flip-tda-o2, rohf-hessian-o2); blocked=[]; merge-failed=[]
- 2026-07-23 harness wave 4: queue merged→done (5); in_progress cleared; roadmaps ○→✓ for ROMP2, ROKS grad, hybrid full-LR B3LYP, ROHF Hessian; SF-TDA already ✓; ratios romp2=1.37 roks-grad=1.117 rpa=0.430 sf-tda=9.56 rohf-hess=perf_exempt; no ⚡ (no standing perf suite); handoff written; queue_remaining=8
- 2026-07-23 tutorial: methods/romp2-o2-sto3g.md
- 2026-07-23 tutorial: methods/roks-grad-o2-sto3g.md
- 2026-07-23 tutorial: methods/h2-b3lyp-rpa-tddft.md
- 2026-07-23 tutorial: methods/spin-flip-tda-o2.md
- 2026-07-23 tutorial: methods/rohf-hessian-o2.md
- 2026-07-23 assign wave 5: adc2-h2co-sto3g (excited-states), mp3-h2-sto3g (post-hf), ccsd-h2-sto3g (post-hf), g0w0-h2-sto3g (excited-states), b2plyp-h2o (dft) → in_progress; skip ccsd-t (needs CCSD), casci, delta-scf; pending_left=3
- 2026-07-23 wave 5 merge: all five green HEADs → dev (adc2-h2co-sto3g, mp3-h2-sto3g, ccsd-h2-sto3g, g0w0-h2-sto3g, b2plyp-h2o); blocked=[]; merge-failed=[]
- 2026-07-23 harness wave 5: queue merged→done (5); in_progress cleared; eom-ee-ccsd unblocked→pending (CCSD energy ships); roadmaps ○→✓ for MP3, CCSD, G0W0 (ADC(2)+B2PLYP already ✓); ratios adc2=0.267 mp3=0.774 ccsd=0.317 g0w0=0.109 b2plyp=1.75; no ⚡ (no standing perf suite); handoff written; queue_remaining=4
- 2026-07-23 tutorial: methods/adc2-h2co-sto3g.md
- 2026-07-23 tutorial: methods/mp3-h2-sto3g.md
- 2026-07-23 tutorial: methods/ccsd-h2-sto3g.md
- 2026-07-23 tutorial: methods/g0w0-h2-sto3g.md
- 2026-07-23 tutorial: methods/b2plyp-h2o.md
- 2026-07-23 assign wave 6: ccsd-t-h2o-sto3g (post-hf), casci-h2-sto3g (post-hf), delta-scf-n2 (excited-states/hf), eom-ee-ccsd (excited-states) → in_progress; skip=[]; pending_left=0 after this wave
- 2026-07-23 wave 6 merge: all four green HEADs → dev (ccsd-t-h2o-sto3g, casci-h2-sto3g, delta-scf-n2, eom-ee-ccsd); blocked=[]; merge-failed=[]
- 2026-07-23 harness wave 6: queue merged→done (4); in_progress cleared; pending empty; roadmaps ○→✓ for CCSD(T), CASCI, EOM-EE-CCSD (ΔSCF already ✓); ratios ccsd-t=0.333 casci=0.076 delta-scf=7.83 eom-ee=0.193; no ⚡ (no standing perf suite); handoff written; queue_remaining=0
- 2026-07-23 tutorial: methods/ccsd-t-h2o-sto3g.md
- 2026-07-23 tutorial: methods/casci-h2-sto3g.md
- 2026-07-23 tutorial: methods/delta-scf-n2.md
- 2026-07-23 tutorial: methods/eom-ee-ccsd.md
- 2026-07-23 assign wave 1 (post re-census): eom-ip-ccsd (excited-states), uks-full-lr-o2-lsda (excited-states), casscf-h2-sto3g (post-hf), ccsd-grad-h2o-sto3g (post-hf), rks-hessian-h2o-sto3g (dft), ump2-grad-o2-sto3g (post-hf) → in_progress; skip=[nevpt2 needs CASSCF, eom-cc2 needs CC2]; pending_left=15

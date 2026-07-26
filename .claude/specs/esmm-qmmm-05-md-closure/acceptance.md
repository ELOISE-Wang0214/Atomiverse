---
spec: esmm-qmmm-05-md-closure
created: 2026-07-26
criteria:
  - id: ac-001
    summary: QM/MM NVE conserves momentum exactly
    type: runtime
    pass_when: |
      regressions/multiscale/qmmm_md_conservation runs 200 NVE steps of
      QM(RHF/STO-3G water, mask column) + charged-LJ MM atoms with
      embed_cutoff ≥ half the box diagonal and |Σp| < 1e-8 amu·Å/fs at every
      sampled step (no drift trend).
    status: pending
  - id: ac-002
    summary: Energy conservation with Verlet scaling
    type: runtime
    pass_when: |
      |ΔE_total|/N < 1e-4 eV/atom over the run, and halving dt reduces the
      drift by 4 ± 1 (O(dt²)); REFERENCE.md records engine tol and the
      tolerance provenance (conservation identities, not snapshots), with
      Frenkel & Smit §4.3 and Senn & Thiel 2009 §3 cited.
    status: pending
  - id: ac-003
    summary: Warm start engaged inside QM/MM MD
    type: runtime
    pass_when: |
      The case asserts production-step SCF iteration counts drop after step 1
      (predictor(3) engaged through the persistent engine).
    status: pending
  - id: ac-004
    summary: Case discipline
    type: code
    pass_when: |
      One flat main(); public headers only (no internal/); ≤ 60 s
      single-threaded; gated ATV_BUILD_HF AND ATV_BUILD_MD; never in default
      ctest; perf_exempt reason recorded.
    status: pending
  - id: ac-005
    summary: Roadmap and tutorials tell the new story
    type: docs
    pass_when: |
      docs/roadmap/multiscale.md marks "QM/MM molecular dynamics" ✓ and the
      Textbook-QM/MM paragraph describes engine slot + mask column + MM-side
      ownership; ch19/ch20 show the mask-column + set_engine + ordinary MM
      styles example; the MC posture (reset_history_each_launch + ASPC
      guard) is documented.
    status: pending
out_of_scope:
  - QM/MM-MC sampling pin; QM/MM × CPMD; PME-QM/MM; GEEP; ONIOM; CUDA
---

# Acceptance — `esmm-qmmm-05-md-closure`

**Done** means the roadmap's "QM/MM molecular dynamics" row is honestly ✓: an
NVE case proves momentum exactness (the third-law fix), energy conservation
with the right dt scaling, and live warm start — and the narrative docs teach
the redesigned surface.

## AC-001 — momentum
The single most direct witness of the Wave-2 physics fix.

## AC-002 — energy + scaling
Separates integrator error from force error; pins the conservation story.

## AC-003 — warm start
The engine-slot payoff: no more cold SCF every step.

## AC-004 — regression discipline
House rules for the e2e tier.

## AC-005 — docs
A redesign the tutorials don't teach didn't happen.

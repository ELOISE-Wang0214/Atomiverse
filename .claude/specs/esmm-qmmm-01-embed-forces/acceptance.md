---
spec: esmm-qmmm-01-embed-forces
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Analytic QM-center embedding leg matches FD
    type: runtime
    pass_when: |
      tests/cpu/scf/test_embedding.cpp pins analytic electron–MM gradient
      w.r.t. QM nuclei against central FD (h = 1e-4 Bohr) ≤ 1e-6 Ha/Bohr,
      and the FD block formerly at gradient.cpp:663-722 is gone (no
      core_matrix rebuild loop in the embedding gradient path).
    status: pending
  - id: ac-002
    summary: Forces on MM sites exist and satisfy Newton's third law
    type: runtime
    pass_when: |
      GradientResult::embedding_forces is filled [n_charges,3] eV/Å;
      analytic site gradient matches FD of E w.r.t. a charge position
      ≤ 1e-6 Ha/Bohr; with embed_self_energy=false,
      |Σ_QM F + Σ_MM F| ≤ 1e-9 eV/Å; embedding_forces is empty when
      point_charges is empty.
    status: pending
  - id: ac-003
    summary: embed_h2o_mm_forces regression matches offline PySCF pin
    type: runtime
    pass_when: |
      regressions/scf/embed_h2o_mm_forces builds under ATV_BUILD_HF, never
      launches PySCF, and the MM-side gradient matches the REFERENCE.md pin
      (PySCF 2.14, examples/qmmm/30-force_on_mm_particles.py recipe)
      element-wise ≤ 1e-5 Ha/Bohr; ≤ 60 s.
    status: pending
  - id: ac-004
    summary: Existing embedding pins unchanged
    type: runtime
    pass_when: |
      regressions/scf/embed_h2o_grad and regressions/scf/h2o_embedding_field
      pass unmodified in tolerance (h2o_embedding_field keeps
      E = −74.98418931920534 Ha with the default embed_self_energy = true).
    status: pending
  - id: ac-005
    summary: One public EmbeddingSite type; nested/duplicate structs deleted
    type: code
    pass_when: |
      grep finds no `PointChargeSite` and no nested `struct PointCharge`
      under include/ or src/; cpu::EmbeddingSite/EmbeddingForce live in
      cpu/component/qm_engine.h with Doxygen units; point_charges is
      std::vector<EmbeddingSite> on MeanFieldMethod.
    status: pending
  - id: ac-006
    summary: ECP+embedding combination is loud, not silently wrong
    type: runtime
    pass_when: |
      If ecp + point_charges is reachable, the analytic path throws a typed
      exception with a guard test; otherwise the spec records
      unreachability evidence in the implementation summary.
    status: pending
out_of_scope:
  - QmEngine/EmbeddingEngine protocol classes (esmm-qmmm-02)
  - Coupling rewrite, mask column, embedding cutoff (esmm-qmmm-03)
  - GEEP / PME-QM/MM / periodic sums; CUDA
---

# Acceptance — `esmm-qmmm-01-embed-forces`

**Done** means point-charge embedding has analytic gradients on both sides of
the QM/MM seam — QM nuclei AND MM sites — with Σ F = 0, routed through
`GradientResult::embedding_forces`, one public `EmbeddingSite` type, an
`embed_self_energy` gate defaulting to the pinned standalone convention, and a
PySCF-pinned regression for the MM-side forces.

## AC-001 — analytic QM-center leg
Kills the O(6N×n_charges) FD sweep; parity gate against the old numerics.

## AC-002 — third law
The physics this stage exists to deliver; momentum-exact QM/MM MD depends on it.

## AC-003 — e2e pin
Offline-oracle discipline: numbers pinned in REFERENCE.md, never live at ctest.

## AC-004 — non-regression
The two shipped embedding pins are the safety net for the rewrite.

## AC-005 — type unification
The template-nested PointCharge hazard (RHF::PointCharge ≠ UKS::PointCharge)
must be gone before esmm-qmmm-02 can erase engines.

## AC-006 — loud unsupported paths
House rule: no silent fallbacks.

---
spec: esmm-qmmm-03-electrostatic-rewrite
created: 2026-07-26
criteria:
  - id: ac-001
    summary: QM region comes from the qm_z mask column, not a ctor list
    type: code
    pass_when: |
      cpu/multiscale/qm_mask.h publishes kAtomQmZ ("atoms","qm_z", int,
      dynamic wire id) with Doxygen how-to; QmmmElectrostatic's ctor takes no
      QM-atom/tag argument (grep); build() throws with guidance when the
      column is missing or all-zero; mask semantics 0=MM, >0=QM with Z=value
      are unit-tested (resolve → rows + Z).
    status: pending
  - id: ac-002
    summary: Engine slot works and rejects non-embeddable engines loudly
    type: runtime
    pass_when: |
      set_engine/set_engine_ref (rvalue owns / lvalue borrows; plain lvalue
      set_engine deleted) accept any QmEngine; build() with AM1 throws a
      message containing the engine name; RHF↔RKS swap unit test runs both;
      the engine object persists across launches (warm-start iteration
      counts drop after step 1 with predictor enabled).
    status: pending
  - id: ac-003
    summary: Exclusion invariant enforced; MM–MM owned once
    type: runtime
    pass_when: |
      build() throws listing tags when any QM atom has nonzero kAtom_charge;
      the coupling sets site self-energy false (MM–MM Coulomb excluded from
      the engine total — unit-tested against the closed-form MM–MM term);
      QmMmVdw touches only QM–MM pairs (unit test).
    status: pending
  - id: ac-004
    summary: Embedding sum has cutoff + minimum image + per-atom charges
    type: runtime
    pass_when: |
      collect_sites honors embed_cutoff under minimum image, reads charges
      from kAtom_charge (the uniform mm_charge ctor arg is gone — grep),
      and omits link MM partners; unit tests pin each.
    status: pending
  - id: ac-005
    summary: Persistent sub-frame; per-run-static region; serial v1
    type: runtime
    pass_when: |
      The private QM Frame is allocated in build() and its address is stable
      across launches (no per-step allocation); refresh() throws on QM-count
      change; build() throws when mpi_size > 1.
    status: pending
  - id: ac-006
    summary: qmmm_es_h2o_pyscf regression pins the full coupling path
    type: runtime
    pass_when: |
      regressions/multiscale/qmmm_es_h2o_pyscf (ATV_BUILD_HF gate, ≤ 60 s,
      public API only, offline pin in REFERENCE.md) matches the PySCF
      qmmm.mm_charge pins for energy, QM forces, AND MM-site forces, and
      |Σ F| ≤ 1e-9 eV/Å; regressions/multiscale/qmmm_es_h2 is deleted with
      superseding intent recorded.
    status: pending
  - id: ac-007
    summary: Contract citizenship (+= / kEnergyElectronic / forces_requested)
    type: runtime
    pass_when: |
      Launch accumulates (pre-seeded kEnergy preserved — unit test); the QM
      share lands in kEnergyElectronic; forces_requested=false skips force
      assembly and site-force scatter; multiscale builds with ATV_BUILD_HF
      alone (CMake re-gating verified by a build config in CI or a recorded
      configure log).
    status: pending
out_of_scope:
  - Mechanical/Polarizable port (esmm-qmmm-04); MD closure (esmm-qmmm-05)
  - PME-QM/MM, switching functions, GEEP, ONIOM, adaptive QM/MM, MPI mode
---

# Acceptance — `esmm-qmmm-03-electrostatic-rewrite`

**Done** means electrostatic QM/MM is the connection layer: any
embedding-capable engine, mask-column region declaration, per-atom MM charges,
cutoff+minimum-image embedding, third-law site forces scattered back, MM–MM
owned solely by the user's force field — composing into MD/Minimize/MC/
Hessian as an ordinary `Potential`.

## AC-001 — mask column
The operator's binding decision; membership + element in one column.

## AC-002 — engine slot
The 自由组合 core: engines are configured by the user and handed over.

## AC-003 — exclusion without new machinery
Zero-charge invariant + self-energy-off = no double counting, provably.

## AC-004 — embedding sum hygiene
Cutoff, images, and real charges — the three toy shortcuts removed.

## AC-005 — lifecycle
No per-step allocation; loud limits (count change, MPI).

## AC-006 — e2e pin
The coupling reproduces the bare-SCF PySCF pins exactly.

## AC-007 — pipeline citizenship
The Wave-1 contract is what makes the coupling co-bindable.

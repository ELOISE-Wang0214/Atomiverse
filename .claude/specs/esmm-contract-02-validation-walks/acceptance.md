---
spec: esmm-contract-02-validation-walks
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Nested HF node satisfies MC energy validation
    type: runtime
    pass_when: |
      tests/cpu/test_mc_scaffold.cpp: an MC whose only energy source is a
      nested Driver<scf::Phase>+HF stub builds successfully; an MC with no
      source still throws.
    status: pending
  - id: ac-002
    summary: Nested HF node satisfies PathDriver force validation
    type: runtime
    pass_when: |
      PathDriver unit test: a band whose only force source is a nested HF
      stub builds; a band with neither component-kind nor HF node still
      throws.
    status: pending
  - id: ac-003
    summary: Both walks are node-aware for_each_pending_kind
    type: code
    pass_when: |
      mc_base.h and path_driver.h validate_components call
      for_each_pending_kind with a dynamic_cast<const HF*> driver-node
      lambda (md.h pattern); no name()-based matching introduced.
    status: pending
  - id: ac-004
    summary: No suite regressions
    type: runtime
    pass_when: |
      Full ctest green; existing MC and path regressions pass unchanged.
    status: pending
out_of_scope:
  - run_slot_ nested dispatch (esmm-contract-03)
  - forces_requested plumbing (esmm-contract-04)
---

# Acceptance — `esmm-contract-02-validation-walks`

**Done** means MC and PathDriver recognize nested Gaussian method drivers as
energy/force sources exactly the way MD/Evaluate/Minimize/Hessian already do,
without weakening the no-source guard.

## AC-001 / AC-002 — recognition parity

The free-combination matrix rows "RHF × MC" and "RHF × NEB" cannot even build
today; after this spec they build (correct execution arrives with
esmm-contract-03/04).

## AC-003 — house pattern

dynamic_cast against kind classes on the node-aware walk — no classifier
functions, no string matching.

## AC-004 — guard intact

The change adds visibility, never permissiveness.

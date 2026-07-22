---
spec: backend-mirror-05-dissolve
created: 2026-07-22
criteria:
  - id: ac-001
    summary: public Backend god class removed
    type: code
    pass_when: |
      No public class named Backend in include/atomiverse/cpu/driver or
      include/atomiverse/cuda/driver; product code does not construct ATV::*::Backend.
    status: pending
  - id: ac-002
    summary: DD lifecycle owned under communication
    type: code
    pass_when: |
      DomainDecomposition ownership/bind used by drivers lives under
      communication/ headers/namespaces, not a driver/backend god header.
    status: pending
  - id: ac-003
    summary: device owns streams/threads setup surface
    type: code
    pass_when: |
      CUDA stream / CPU thread configuration used at run start is reached via
      device/ APIs (or Driver calling device/), not Backend::init god method.
    status: pending
  - id: ac-004
    summary: no_backend regression runs a preset
    type: runtime
    pass_when: |
      regressions/layout/no_backend/ builds and runs a public preset (Evaluate or
      dynamics MD) without including or naming Backend.
    status: pending
  - id: ac-005
    summary: full suite green after dissolve
    type: runtime
    pass_when: |
      Full check + ctest pass.
    status: pending
out_of_scope:
  - Path renames; Criterion rename (done in 04); docs
---

# Acceptance — backend-mirror-05-dissolve

Done means Backend is gone, responsibilities live in communication/device/asynchronous/driver, tests rewritten.

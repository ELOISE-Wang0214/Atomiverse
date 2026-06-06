---
slug: async-analysis-snapshot-ring-01-uniformity
criteria:
  - id: ac-001
    summary: Delivered snapshot slot is internally consistent for its step
    type: scientific
    pass_when: |
      test_snapshot_publisher's freeze-then-dispatch regression (every-step
      capture + slowed on_snapshot over analytic trajectory x_i(step)=f(tag_i,step))
      asserts, for every delivered slot, slot.x[i]==f(slot.tag[i], slot.step) for
      all i; test passes.
    status: pending
  - id: ac-002
    summary: BackpressurePolicy enum and Config::policy are removed
    type: code
    pass_when: |
      grep finds no `enum class BackpressurePolicy` in core/types.h and no
      `policy` member in AnalysisExecutor::Config; src/cpu/backend.cpp has no
      cfg.policy assignment.
    status: verified
    last_checked: 2026-06-06
  - id: ac-003
    summary: acquire_slot is pure drop_oldest with constant ring depth K
    type: code
    pass_when: |
      grep finds no grow_ring_locked, wait_for_free_slot_locked, cv_slot_free_,
      ring_grew_, ring_size_grew, or policy_ symbols in analysis_executor.{h,cpp};
      acquire_slot_locked has a single drop_oldest path; ring_capacity stays
      constant after init.
    status: verified
    last_checked: 2026-06-06
  - id: ac-004
    summary: Component-level backpressure override removed
    type: code
    pass_when: |
      grep finds no `backpressure()` declaration in
      include/atomiverse/cpu/async_component_base.h.
    status: verified
    last_checked: 2026-06-06
  - id: ac-005
    summary: Capture mask computed once, no per-step recompute
    type: code
    pass_when: |
      grep finds no mask_union_firing symbol in src/cpu/snapshot_publisher.cpp;
      Publish caches mask_union_all and reuses it across launch() calls.
    status: verified
    last_checked: 2026-06-06
  - id: ac-006
    summary: Block/Grow backpressure tests removed; drop_oldest + global-K asserted
    type: code
    pass_when: |
      tests/cpu/test_analysis_executor.cpp contains no BlockBackpressure or
      GrowBackpressure test and no cfg.policy assignment; round-robin +
      n_dropped + global ring-depth K assertions remain.
    status: verified
    last_checked: 2026-06-06
  - id: ac-007
    summary: Config::ring_size documents global K semantics
    type: docs
    pass_when: |
      Doxygen on Config::ring_size states it is the single global ring depth K
      (larger K = fewer drops, more memory).
    status: pending
  - id: ac-008
    summary: Full CPU test suite passes
    type: runtime
    pass_when: |
      `cmake -B build -DATV_BUILD_TESTS=ON && cmake --build build -j &&
      ctest --test-dir build --output-on-failure` (ATV_BUILD_CUDA=OFF, login
      node) exits 0.
    status: verified
    last_checked: 2026-06-06
---

# Acceptance criteria

- **ac-001** binds the Domain-basis invariant: the freeze-then-dispatch RED regression is the falsifiable scientific check; it must be written failing first (analytic trajectory, slowed worker) and pass after the uniform path lands.
- **ac-002 … ac-005** lock the structural collapse: enum gone, executor single-path with constant K, no component override, mask computed once.
- **ac-006** ensures the test suite reflects the single policy (removed cases double as compile-fences against reintroducing block/grow).
- **ac-007** covers the K-semantics docstring required for link 02 consumers.
- **ac-008** is the integration gate on the login node CPU build.

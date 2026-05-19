---
slug: parallel-md-finish-06-close-out
criteria:
  - id: ac-001
    summary: Full ctest cpu+mpi labels exit green with no regressions
    type: runtime
    pass_when: |
      `cmake -B build -DCMAKE_BUILD_TYPE=Release
      -DATV_BUILD_TESTS=ON -DATV_BUILD_BENCHMARKS=ON
      -DATV_BUILD_RHF=ON -DATV_USE_MPI=ON && cmake --build
      build -j && ctest --test-dir build -L "cpu|mpi"
      --output-on-failure` exits 0. Among the tests listed,
      both `test_nose_hoover_t_ramp` and `test_md_no_alloc`
      appear (the latter SKIPPED on Darwin is acceptable).
      The total `tests passed` count is ≥ the pre-Phase-5
      baseline (20 MPI tests + the documented CPU suite size
      at the parent spec's Status snapshot timestamp).
      (Closes parent ac-025.)
    status: verified
    last_checked: 2026-05-19

  - id: ac-002
    summary: parent parallel-md-engine.md flipped to retired
    type: docs
    pass_when: |
      `.claude/specs/parallel-md-engine.md` frontmatter
      contains `status: retired` (was `in-progress`). The
      body retains a closure note dated 2026-05-19 or later
      that references the six chain sub-specs by full slug.
    status: verified
    last_checked: 2026-05-19

  - id: ac-003
    summary: All chain sub-specs flipped to retired in same commit
    type: docs
    pass_when: |
      Each of `parallel-md-finish-{01,02,03,04,05,06}-*.md`
      frontmatter has `status: retired` after the chain
      closure commit. The git log shows a single commit
      whose subject matches `(?i)retire.*chain` and whose
      diff touches all 7 files (parent + 6 chain).
    status: verified
    last_checked: 2026-05-19

  - id: ac-004
    summary: INDEX.md reflects retirement
    type: docs
    pass_when: |
      `.claude/specs/INDEX.md` contains the substring
      `[retired — see parallel-md-finish` (the chain
      pointer) AND no longer contains `parallel-md-engine`
      under any `in-progress` heading. The six chain
      entries are listed (either indented under the parent
      line or as their own section).
    status: verified
    last_checked: 2026-05-19
---

# Acceptance criteria

- **ac-001** is the load-bearing acceptance: it asserts that everything still works after the Phase 5 work. The test-count floor (≥ pre-Phase-5 baseline) guards against silently dropped tests, and the explicit `test_nose_hoover_t_ramp` + `test_md_no_alloc` mentions lock in the new infrastructure tests.
- **ac-002 / ac-003** retire the parent and chain together, atomically. Either both happen or neither — the chain shouldn't outlive a still-active parent.
- **ac-004** ensures the INDEX (the single source of truth read by humans + downstream specs) reflects reality.

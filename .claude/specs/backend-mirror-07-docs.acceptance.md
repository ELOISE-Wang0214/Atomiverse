---
spec: backend-mirror-07-docs
created: 2026-07-22
criteria:
  - id: ac-001
    summary: CLAUDE Architecture matches new layout vocabulary
    type: docs
    pass_when: |
      CLAUDE.md Architecture section describes frame as the data pillar, core contracts
      (no facade/), full-word backend dirs (dynamics, communication, longrange, …),
      Driver/Criterion (not CPUDriver/LoopHandler), and no runtime/loop/shared layers.
    status: pending
  - id: ac-002
    summary: AGENTS Architecture matches CLAUDE (no Handler dual story)
    type: docs
    pass_when: |
      AGENTS.md Architecture uses Criterion/add_criterion and Driver naming consistent
      with CLAUDE; does not document LoopHandler/add_handler or CPUDriver/CUDADriver
      as current API.
    status: pending
  - id: ac-003
    summary: architecture.md inventory refreshed
    type: docs
    pass_when: |
      .claude/notes/architecture.md inventory lists the post-rename directories and
      does not present facade/, top-level longrange/, or cuda/kspace/ as current.
    status: pending
  - id: ac-004
    summary: banned legacy terms absent from arch docs
    type: code
    pass_when: |
      In CLAUDE.md and AGENTS.md Architecture sections, rg finds no current-API claims
      for include paths atomiverse/facade, cuda/kspace as canonical, or CPUDriver as
      the user-facing loop type.
    status: pending
  - id: ac-005
    summary: suite still green after doc edits
    type: runtime
    pass_when: |
      Full check + ctest pass (docs-only change must not break build).
    status: pending
out_of_scope:
  - Code layout changes; torch redesign; per-example tutorial rewrites beyond path fixes already done
---

# Acceptance — backend-mirror-07-docs

Done means harness docs tell one story that matches the mirrored full-word layout and new Driver/Criterion vocabulary.

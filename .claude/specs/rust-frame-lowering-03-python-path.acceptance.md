---
slug: rust-frame-lowering-03-python-path
criteria:
  - id: ac-001
    summary: run() accepts a molrs.Frame and the Frame binding is removed
    type: code
    pass_when: |
      The user-facing driver run() entry accepts a molrs.Frame handle;
      python/src/atomiverse/_convert.py no longer exists; no CpuFrame type is
      exposed on the Atomiverse Python binding surface.
    status: pending
  - id: ac-002
    summary: a Python MD run starts from a molrs.Frame end to end
    type: runtime
    pass_when: |
      A Python integration test constructs a molrs.Frame, runs cpu::md::MD for
      several steps, and finds finite energy plus updated positions written back
      into that same molrs.Frame; the test passes.
    status: pending
  - id: ac-003
    summary: AIMD still conserves energy through the new ingress path
    type: runtime
    pass_when: |
      An H2 STO-3G NVE AIMD run (MD + SCFEvaluator on one materialized
      cpu::Frame) keeps the E_total peak-to-peak envelope within the aimd
      spec bound; the test passes under ctest.
    status: pending
  - id: ac-004
    summary: full test suite stays green
    type: runtime
    pass_when: |
      `cmake --build build` succeeds and `ctest --test-dir build` reports all
      tests passing.
    status: pending
---

# Acceptance criteria

- **ac-001** is the structural deliverable — `run()` is on `molrs.Frame`, and the
  old hand-rolled Frame binding (`_convert.py` + `CpuFrame`) is gone. Verified by
  inspecting the source tree, hence `code`.
- **ac-002** proves the end-to-end Python ingress path works on the real
  shared-handle mechanism.
- **ac-003** guards the AIMD shared-frame invariant — the decisive constraint
  behind one Frame per backend — through the rewritten ingress.
- **ac-004** is the regression gate.

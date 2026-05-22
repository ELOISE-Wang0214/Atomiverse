---
slug: rust-frame-lowering-01-molrs-bridge
criteria:
  - id: ac-001
    summary: molrs-cxxapi bridge exposes a zero-copy Frame read + create/write interface
    type: code
    pass_when: |
      molrs-cxxapi/src/bridge.rs declares Frame column-read functions returning
      rust::Slice per dtype (f64/i64/...), block/column introspection, a simbox
      accessor, plus frame_new and create-or-update frame_set_column_* writers;
      molrs-cxxapi/src/lib.rs implements them.
    status: pending
  - id: ac-002
    summary: molrs-python exposes the underlying Frame handle as a PyCapsule
    type: code
    pass_when: |
      molrs-python's Frame pyo3 class has a method returning a PyCapsule that
      carries the underlying *mut Frame, usable to borrow the same Frame from
      another FFI.
    status: pending
  - id: ac-003
    summary: bridge Frame read/write round-trips under cargo test
    type: runtime
    pass_when: |
      A Rust test in the molrs workspace builds a Frame via frame_new, writes
      f64 and i64 columns plus a simbox, reads them back through the bridge
      accessors, and asserts element-wise equality; `cargo test` passes.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002** are the deliverable surface — the CXX bridge gains a Frame
  ingress it never had, and molrs-python gains the shared handle `-03` needs.
  Both are structural, verified by reading the source, hence `code`.
- **ac-003** proves the bridge actually round-trips data, exercised by the molrs
  workspace's own `cargo test`.

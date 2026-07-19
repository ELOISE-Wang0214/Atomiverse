---
spec: python-on-demand-install-04-docs-regression
created: 2026-07-19
criteria:
  - id: ac-001
    summary: getting_started documents pure-pip then on-demand install
    type: docs
    pass_when: |
      Installation section states pure pip does not compile native code and
      documents second-phase python -m atomiverse install with copy-paste example.
    status: pending
  - id: ac-002
    summary: Document --config XOR --md/--global mutual exclusion
    type: docs
    pass_when: |
      Docs state mutual exclusion and show at least one illegal combination.
    status: pending
  - id: ac-003
    summary: Document version lock and --src/--url
    type: docs
    pass_when: |
      Docs mention version-lock rejection and document --src and --url.
    status: pending
  - id: ac-004
    summary: Document module list and --hf rejection
    type: docs
    pass_when: |
      Docs list installable modules consistent with 03 and reject --hf.
    status: pending
  - id: ac-005
    summary: Document single-env single-binary overwrite
    type: docs
    pass_when: |
      Docs state re-install overwrites the single native binary in the env.
    status: pending
  - id: ac-006
    summary: Root README points to Python install docs
    type: docs
    pass_when: |
      README.md links or points to python/docs/getting_started.md for Python install.
    status: pending
  - id: ac-007
    summary: Docs flags/fields match 03 unified config schema
    type: docs
    pass_when: |
      Every install flag/key in getting_started matches 03 schema surface.
    status: pending
  - id: ac-008
    summary: install_cli regression passes without full native compile
    type: runtime
    pass_when: |
      reg_install_cli.sh exits 0 in ≤60s asserting help, mutual exclusion,
      --hf reject, version mismatch; no full native compile required.
    status: pending
  - id: ac-009
    summary: Regression wired into regressions target
    type: code
    pass_when: |
      regressions/CMakeLists.txt defines run_reg_python_install_cli dependency
      of regressions target; README table lists python/install_cli.
    status: pending
  - id: ac-010
    summary: REFERENCE.md records canonical CLI contract
    type: docs
    pass_when: |
      REFERENCE.md states literature N/A and expected exit semantics for
      mutual exclusion, --hf, version mismatch.
    status: pending
out_of_scope:
  - Implementing CLI/CMake (01-03)
  - manylinux / SCF enablement
---

# Acceptance — python-on-demand-install-04-docs-regression

完成态：文档描述两阶段安装契约；`regressions/python/install_cli` 钉住公开 CLI 行为且无需完整原生编译。

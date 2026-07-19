---
spec: python-on-demand-install-02-pure-wheel
created: 2026-07-19
criteria:
  - id: ac-001
    summary: Unit tests cover pure import and install-message ImportError
    type: code
    pass_when: |
      python/tests/test_pure_import.py asserts soft import, __version__,
      ImportError containing "python -m atomiverse install", and extension-present path.
    status: pending
  - id: ac-002
    summary: _binding soft-imports missing extension without raising at import
    type: code
    pass_when: |
      Importing atomiverse._binding without extension does not raise;
      require_cpp / _cpp attr raises install-oriented ImportError.
    status: pending
  - id: ac-003
    summary: Package root soft-loads and exposes __version__
    type: code
    pass_when: |
      Without extension, import atomiverse succeeds with __version__; from
      atomiverse import MD raises install ImportError. With extension, surface
      and PW hasattr behavior preserved.
    status: pending
  - id: ac-004
    summary: python -m atomiverse stub is importable and install-oriented
    type: runtime
    pass_when: |
      python -m atomiverse runs without missing-so ImportError; output mentions
      install guidance; no uncaught traceback.
    status: pending
  - id: ac-005
    summary: Default wheel build is pure Python without CMake or nanobind
    type: runtime
    pass_when: |
      wheel.cmake false; build-system.requires lacks nanobind; backend remains
      scikit_build_core; pure wheel has no _atomiverse shared library.
    status: pending
  - id: ac-006
    summary: Regression pure_import passes in isolated pure-Python path
    type: runtime
    pass_when: |
      regressions/python/pure_import/reg_pure_import.sh exits 0 with isolated
      pure sources; verifies import, compute ImportError, and -m guidance.
    status: pending
  - id: ac-007
    summary: Existing extension-present Python tests still pass
    type: runtime
    pass_when: |
      With native extension available, existing python/tests surface suite passes.
    status: pending
  - id: ac-008
    summary: Full library check and test suite green
    type: runtime
    pass_when: |
      Project check/test path used by implementer is green including pure-import tests.
    status: pending
out_of_scope:
  - Full install CLI (03)
  - CMake modularization (01)
---

# Acceptance — python-on-demand-install-02-pure-wheel

完成态：pure wheel 默认可装；无 so 时包/CLI stub 可用；计算路径 fail-loud 指向 install。

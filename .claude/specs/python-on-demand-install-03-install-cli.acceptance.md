---
spec: python-on-demand-install-03-install-cli
created: 2026-07-19
criteria:
  - id: ac-001
    summary: Mutual exclusion of --config vs inline module flags fails loud
    type: runtime
    pass_when: |
      install with both --config and --md exits non-zero naming mutual exclusion.
    status: verified
    last_checked: 2026-07-19
  - id: ac-002
    summary: Empty install with no modules and no config fails loud
    type: runtime
    pass_when: |
      install with neither --config nor enabling flags exits non-zero.
    status: verified
    last_checked: 2026-07-19
  - id: ac-003
    summary: --hf / modules.hf rejected until SCF Python bindings exist
    type: runtime
    pass_when: |
      Any hf enable path exits non-zero stating SCF/HF bindings unavailable.
    status: verified
    last_checked: 2026-07-19
  - id: ac-004
    summary: Inline k:v parse and bool synonyms map via binding table
    type: runtime
    pass_when: |
      --md mpi:yes and true/on/1 map to ON via BINDING_TABLE; first-colon split tested.
    status: verified
    last_checked: 2026-07-19
  - id: ac-005
    summary: md auto-enable closure turns on potentials and neighbor
    type: runtime
    pass_when: |
      resolve_defines for md includes potentials+neighbor ATV_* defines (01 DAG).
    status: verified
    last_checked: 2026-07-19
  - id: ac-006
    summary: Bound keys cannot be passed as unbound -D
    type: runtime
    pass_when: |
      Binding-table key under cmake.unbound or -D exits non-zero.
    status: verified
    last_checked: 2026-07-19
  - id: ac-007
    summary: Version from importlib.metadata; default tarball URL uses it
    type: code
    pass_when: |
      source.py uses importlib.metadata.version; default URL
      https://github.com/Roy-Kid/Atomiverse/archive/refs/tags/v{version}.tar.gz;
      no install --version flag.
    status: verified
    last_checked: 2026-07-19
  - id: ac-008
    summary: Source pyproject version mismatch rejects install
    type: runtime
    pass_when: |
      Mismatched python/pyproject.toml version fails before pip is invoked.
    status: verified
    last_checked: 2026-07-19
  - id: ac-009
    summary: Install always sets ATV_BUILD_PYTHON=ON in cmake.define
    type: code
    pass_when: |
      install uses pip -C cmake.define with ATV_BUILD_PYTHON=ON; no hand cmake -B path.
    status: verified
    last_checked: 2026-07-19
  - id: ac-010
    summary: status prints enabled modules and define fingerprint
    type: runtime
    pass_when: |
      python -m atomiverse status exits 0 with module list and define fingerprint.
    status: verified
    last_checked: 2026-07-19
  - id: ac-011
    summary: CLI package loads without importing the C++ extension
    type: code
    pass_when: |
      cli modules import without atomiverse._binding / _atomiverse at import time.
    status: verified
    last_checked: 2026-07-19
  - id: ac-012
    summary: Regression dry-run exercises public CLI without C++ compile
    type: runtime
    pass_when: |
      regressions/python/install_cli_dryrun runs without C++ compile; checks help,
      mutual exclusion, empty/HF reject, dry-run defines for --md including closure
      and ATV_BUILD_PYTHON=ON per REFERENCE.md.
    status: verified
    last_checked: 2026-07-19
  - id: ac-013
    summary: Full Python install-cli unit suite green
    type: runtime
    pass_when: |
      pytest python/tests/test_install_cli.py passes.
    status: verified
    last_checked: 2026-07-19
out_of_scope:
  - manylinux prebuilt wheels
  - multi-profile cache
  - SCF/HF bindings
---

# Acceptance — python-on-demand-install-03-install-cli

完成态：install/status 在纯 Python 层正确可测；scikit-build 为唯一编译入口；配置统一且 fail-loud。

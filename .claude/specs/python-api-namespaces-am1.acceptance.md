---
slug: python-api-namespaces-am1
criteria:
  - id: ac-001
    summary: Namespaced package layout exists; flat md.py deleted
    type: code
    pass_when: |
      python/src/atomiverse/ contains md/{__init__,integrator,thermostat}.py,
      neighbor.py, potentials/{__init__,pair,bond,angle,dihedral}.py,
      recorder.py, driver.py, semiempirical/__init__.py; the flat file
      python/src/atomiverse/md.py no longer exists.
    status: verified
    last_checked: 2026-07-15
  - id: ac-002
    summary: Top-level names are identity re-exports of deep canonical paths
    type: runtime
    pass_when: |
      python/tests/test_namespaces.py passes: for every row of the spec's
      mapping table (16 entries), `atv.<TopName> is <deep.path.Name>` is True
      (e.g. atv.PairLJCut is atv.potentials.pair.LJCut,
      atv.VelocityVerlet is atv.md.integrator.VelocityVerlet).
    status: verified
    last_checked: 2026-07-15
  - id: ac-003
    summary: __all__ equals exactly the 16-name honest top-level surface
    type: runtime
    pass_when: |
      Updated python/tests/test_surface.py passes with set(atomiverse.__all__)
      == {MD, VelocityVerlet, InitVelocity, NoseHoover, PairLJCut,
      NeighborList, BondHarmonic, AngleHarmonic, DihedralMultiHarmonic,
      Thermo, pw_scf, PWResult, AM1Config, AM1Evaluator, AM1Result,
      single_point}; dead-symbol and fail-loud gates still pass.
    status: verified
    last_checked: 2026-07-15
  - id: ac-004
    summary: device= parameter — "cpu" default works, "cuda" raises clear error
    type: runtime
    pass_when: |
      atv.md.MD(dt=0.001) and atv.md.MD(dt=0.001, device="cpu") run the
      conftest LJ smoke; atv.md.MD(device="cuda") and
      atv.driver.single_point(frame, ev, device="cuda") raise ValueError whose
      message states CUDA bindings are not yet available.
    status: verified
    last_checked: 2026-07-15
  - id: ac-005
    summary: AM1 single_point plumbing — energy on scalars, charges written back
    type: runtime
    pass_when: |
      On the tetrahedral CH4 molrs.Frame, atv.single_point(frame, AM1Evaluator())
      leaves frame scalars energy finite and < 0 [eV], and the frame's atoms
      charge column plus evaluator.last_result().mulliken_charges both carry
      5 charges; last_result().converged is True.
    status: verified
    last_checked: 2026-07-15
  - id: ac-006
    summary: AM1 CH4 Mulliken charges match the pinned AM1 reference oracle
    type: scientific
    pass_when: |
      Via the Python single_point entry, tetrahedral CH4 (C-H 1.0897 A) gives
      Mulliken charges C = -0.262125 e and H = +0.065531 e within 2.5e-3 e
      (AM1: Dewar et al. 1985, DOI 10.1021/ja00299a024; oracle values already
      cited in tests/cpu/semiempirical/test_am1_scaffold.cpp — no new
      reference values introduced).
    status: pending
  - id: ac-007
    summary: Unsupported element fails loud across the binding boundary
    type: runtime
    pass_when: |
      single_point on a frame containing Z=99 raises RuntimeError in Python
      (no crash, no fallback), per the no-fallback rule.
    status: verified
    last_checked: 2026-07-15
  - id: ac-008
    summary: Existing Python MD suite green on the restructured package
    type: runtime
    pass_when: |
      pytest python/tests passes in the ~/work/.x86_64 venv with the rebuilt
      extension: conftest.py, test_md.py (deep-path rewrite + one top-level
      convenience case), test_namespaces.py, test_semiempirical.py,
      test_surface.py all green.
    status: verified
    last_checked: 2026-07-15
  - id: ac-009
    summary: molpy + atomiverse interop check in the x86_64 venv
    type: runtime
    evaluator_hint: acceptance-manual (venv ~/work/.x86_64; never ctest)
    pass_when: |
      In the ~/work/.x86_64 venv, python/tests/test_molpy_interop.py passes:
      molpy (source ~/work/molcrafts/molpy) and atomiverse import together,
      a molpy-built system is handed to atv via the shared molrs Frame, a
      10-step LJ MD runs, and the written-back energy is finite. The test
      skips (importorskip) where molpy is absent.
    status: verified
    last_checked: 2026-07-15
  - id: ac-010
    summary: C++ build and full unit suite unaffected; bind_md.cpp untouched
    type: code
    pass_when: |
      cmake build with ATV_BUILD_PYTHON=ON compiles bind_semiempirical.cpp;
      ctest --test-dir build is fully green with no new/modified tests under
      tests/; git diff for python/src/bindings/bind_md.cpp is empty (each
      add_* remains a discrete lambda site for chain link 04-embed).
    status: verified
    last_checked: 2026-07-15
  - id: ac-011
    summary: Mapping table and canonical paths documented
    type: docs
    pass_when: |
      atomiverse/__init__.py docstring contains the full deep-path <-> top-level
      mapping table; python/docs/getting_started.md shows both the canonical
      deep-path style and the top-level convenience style, plus the AM1
      single-point and device="cpu" examples; every new public symbol in
      bind_semiempirical.cpp carries @brief/@param (with units)/@return.
    status: pending
---

# Acceptance criteria

- **ac-001/002/003** 锁定重构本体：深层为真源、顶层仅同一性再导出、诚实面恰 16 名。
- **ac-004** 锁定 device= 用户裁决（默认 cpu 可用、cuda 明确报错、无回退）。
- **ac-005/006/007** 锁定 AM1 单点：管路（能量/电荷双通道 + write_back 契约）、科学 oracle（沿用已引证参考值）、fail-loud 边界。
- **ac-008/009** 是 venv 级运行时验收（含 molpy 互操作，绝不进 ctest）。
- **ac-010** 保护 C++ 三层测试铁律与 launcher 链 04-embed 的绑定位点前提。
- **ac-011** 文档合同：映射表、双写法示例、Doxygen 单位。

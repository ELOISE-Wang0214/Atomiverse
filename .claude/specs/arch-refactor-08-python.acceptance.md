---
spec: arch-refactor-08-python
created: 2026-07-10
criteria:
  - id: ac-001
    summary: "Dead SCF/optimizer facade files are gone; bind.cpp comment updated"
    type: code
    pass_when: "python/src/atomiverse/scf.py, python/src/atomiverse/optimizer.py, python/src/bindings/bind_scf.cpp and python/src/bindings/bind_optimizer.cpp do not exist in the tree, and python/src/bindings/bind.cpp no longer claims bind_scf/bind_optimizer are 'disabled pending migration' (comment states SCF/optimizer bindings are future work)."
    status: pending
  - id: ac-002
    summary: "No dead _cpp references or legacy/dead symbols anywhere under python/"
    type: code
    evaluator_hint: "rg over python/ excluding python/tests/test_surface.py"
    pass_when: "rg -n '_cpp\\.(LBFGS|LBFGSStep|MinimizerPhase|SCFConfig|SpinType|MeanFieldType|XCType|scf_solve|scf_gradient|build_integrals_sto3g)|av\\.(SCF|RHF|UHF|UKS|GeometryOptimizer|Frame)\\b|CpuFrame' python/ --glob '!tests/test_surface.py' returns no matches (exit code 1)."
    status: pending
  - id: ac-003
    summary: "Python MD class names match canonical C++ names, no legacy aliases"
    type: code
    pass_when: "python/src/bindings/bind_md.cpp binds the Python class names 'PairLJCut', 'BondHarmonic', 'AngleHarmonic', 'DihedralMultiHarmonic' (and none of 'LJ', 'HarmonicBond', 'HarmonicAngle', 'MultiHarmonicDihedral' as nb::class_ names); atomiverse.__all__ lists exactly the canonical MD surface (MD, VelocityVerlet, InitVelocity, NoseHoover, PairLJCut, NeighborList, BondHarmonic, AngleHarmonic, DihedralMultiHarmonic, Thermo) plus pw_scf and PWResult."
    status: pending
  - id: ac-004
    summary: "Default-config wheel builds from python/ via scikit-build-core"
    type: runtime
    evaluator_hint: "pip install --no-build-isolation -ve python/"
    pass_when: "pip install --no-build-isolation -ve python/ exits 0 with the default CMake defines (ATV_BUILD_RHF=ON, ATV_BUILD_PW=OFF) and `python -c \"import atomiverse\"` succeeds afterwards."
    status: pending
  - id: ac-005
    summary: "Full Python test suite green, including the new surface smoke tests"
    type: runtime
    evaluator_hint: "python -m pytest python/tests"
    pass_when: "pytest python/tests exits 0: all tests in test_md.py pass under the canonical names, test_surface.py's MD LJ 10-step smoke passes (positions change, energy and kinetic finite), and the pw_scf toy smoke either passes (PW built) or is skipped exactly by the hasattr(_cpp, 'pw_scf') gate (PW off) — never errors."
    status: pending
  - id: ac-006
    summary: "Dead surface fails loud — imports and attributes raise immediately"
    type: runtime
    evaluator_hint: "pytest node: python/tests/test_surface.py (dead-surface tests)"
    pass_when: "importlib.import_module('atomiverse.scf') and importlib.import_module('atomiverse.optimizer') raise ModuleNotFoundError; getattr(atomiverse, name) raises AttributeError for each of SCF, RHF, UHF, UKS, GeometryOptimizer, Frame, LJ, HarmonicBond, HarmonicAngle, MultiHarmonicDihedral; when the extension lacks pw_scf, calling atomiverse.pw_scf(...) raises RuntimeError whose message mentions ATV_BUILD_PW."
    status: pending
  - id: ac-007
    summary: "Python docs describe only the honest MD+PW surface"
    type: code
    pass_when: "python/docs/getting_started.md and the atomiverse package docstring build the quickstart on molrs.Frame with canonical component names, contain no occurrence of CpuFrame, av.Frame, RHF/SCF/GeometryOptimizer usage, and state the roadmap note that SCF and geometry-optimization bindings are future work."
    status: pending
out_of_scope:
  - "正经的 Minimize/SCF Python 绑定（后续独立 spec）"
  - "链内 C++ 侧改动（include/ src/ tests/ 不触碰）"
  - "CUDA/torch wheel 变体、PyPI 发布流程"
---

# Acceptance — arch-refactor-08-python

- **ac-001/ac-002** 拆除面完整性：文件级删除 + 全树 token 门禁（token 集合与
  test_surface.py 内置 grep 门禁一致；排除该测试文件自身以避免模式串自匹配）。
- **ac-003** clean-break 命名：Python 名一次性对齐 canonical C++ 名，无 alias。
- **ac-004/ac-005** 构建与运行契约：默认 wheel 构建 + pytest 全绿；PW smoke
  的 pass/skip 语义由 `hasattr(_cpp, "pw_scf")` 精确决定。
- **ac-006** fail-loud 契约：死面必须立即、明确地报错，杜绝假可用回潮。
- **ac-007** 文档诚实性：文档面与实际绑定面一致，并明示 roadmap。

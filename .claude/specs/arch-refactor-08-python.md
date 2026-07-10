---
title: arch-refactor-08-python — Python 绑定面收口：删除假可用 SCF/optimizer 面，canonical 命名对齐
status: approved
date: 2026-07-10
chain: arch-refactor (08/09)
depends_on: []
---

# arch-refactor-08-python — Python 绑定面收口

## Summary

Python wheel 当前暴露一个大面积"假可用"表面：`_atomiverse` 扩展只编译
bind.cpp / bind_md.cpp / bind_pw.cpp（根 CMakeLists.txt:718-722），但
`atomiverse.scf`（import 即在 :209-214 处因未绑定的 `_cpp.XCType` 炸掉）、
`atomiverse.optimizer`（`GeometryOptimizer.run()` 调用不存在的
`_cpp.LBFGS`/`_cpp.LBFGSStep`/`_cpp.MinimizerPhase`，并硬编码返回
`converged=True, final_energy=0.0`）以及 5 个引用 `av.RHF`/`av.UHF`/`av.UKS`/
`av.SCF`/`av.GeometryOptimizer`/`av.Frame` 的 examples 仍留在树上；未编译的
bind_scf.cpp / bind_optimizer.cpp 已与 C++ 现状脱节（把 `cpu::minimizer::LBFGS`
组件当 driver 绑、引用已不存在的 `LBFGSStep`、参数错位）。同时 bind_md.cpp 把
canonical C++ 名映射回 legacy Python 名（`PairLJCut`→"LJ" 等），违反
"LAMMPS-style names are canonical" 既定规则。本 spec 一次性 clean-break：
删除全部假可用面及其 examples，Python 类名对齐 canonical C++ 名（无 alias
过渡），`__init__.py`/文档收口为诚实的 MD+PW 面，并加 Python smoke 测试 +
grep 门禁防止假可用面回潮。逻辑上独立于链内 02-07，可并行实施。

## Domain basis

本 spec 不引入新物理，纯 API/绑定面重构。smoke 测试沿用已验证的物理参数与
单位约定（eV, Å, fs, amu, K；CLAUDE.md § Conventions）：

- MD smoke：Ar Lennard-Jones，ε = 0.0104 eV，σ = 3.4 Å（现有
  python/tests/conftest.py fixture 参数），只断言能量有限、KE > 0、坐标演化，
  不引入新的参考值。
- PW smoke：空 `upf_path` 走 toy 局域势（CLAUDE.md 明示 "runs, unphysical"），
  只断言可运行（`iterations >= 1`、能量为有限 float [Ry]），不做数值物理断言。

## Design

**删除的实体（假可用面）：**

- `atomiverse.scf` 模块（`SCF`/`RHF`/`UHF`/`UKS`/`Integrals`/`SCFResult`/
  `GradientResult` 全部死代码）与 `atomiverse.optimizer`
  （`GeometryOptimizer`/`OptResult`/`OptRecord`）。
- `python/src/bindings/bind_scf.cpp`、`bind_optimizer.cpp`：未列入编译且三重
  过期——bind_optimizer.cpp 把 `cpu::minimizer::LBFGS`（lbfgs.h:79，
  ComponentBase）绑成带 `.run()/.build()/.on()` 的 driver（违反"算法是组件"，
  真 driver 是 `cpu::minimizer::Minimize`，minimize.h:76）、绑定 C++ 中已不
  存在的 `LBFGSStep`、且 `LBFGSStep(memory, alpha0, ftol, max_ls)` 与
  optimizer.py 传 `(memory, ftol)` 参数错位（ftol 落到 alpha0）。将来重绑是
  独立 spec（见 Out of scope），届时以 `Minimize` driver + `LBFGS` 组件 +
  `SCFEvaluator` 的正确形状起草。
- 引用不存在符号的 examples，逐一定案：
  - `01_h2_rhf.py`（av.RHF）、`02_h2_uhf.py`（av.UHF/av.RHF，已核实）、
    `03_uks.py`（av.UKS/av.UHF，已核实）、`04_aspirin_opt.py`
    （av.SCF/av.GeometryOptimizer）、`single_point_bench.py`（av.SCF）→ **删**。
  - `04_aspirin_bench.py` → **删**（经 importlib 硬依赖 04_aspirin_opt.py）。
  - `04_aspirin_opt_pyscf.py` → **删**（PySCF-only 对照基线，其 Atomiverse
    对照物已删，单独保留具有误导性）。
  - `pw_qe_compare.py` → **改写**（PW 面真实可用；仅其 `build_si_diamond`
    使用已删除的 `av.Frame` CpuFrame API，改为 molrs.Frame 构造，
    `pseudopotentials/` 数据保留）。
  - `examples/utils/`（timing.py）→ **删**（已核实无任何 example 实际 import，
    孤儿助手）。

**命名 clean-break（Python 名 = canonical C++ 名）：**

| 旧 Python 名 | 新 Python 名（= C++ 名） |
|---|---|
| `LJ` | `PairLJCut` |
| `HarmonicBond` | `BondHarmonic` |
| `HarmonicAngle` | `AngleHarmonic` |
| `MultiHarmonicDihedral` | `DihedralMultiHarmonic` |

`MD`/`VelocityVerlet`/`InitVelocity`/`NoseHoover`/`NeighborList`/`Thermo`/
`Phase` 已与 C++ 名一致，保持不变。bind_md.cpp 中 PhaseBuilder 的内部注册
方法名同步（`add_lj`→`add_pair_lj_cut`、`add_harmonic_bond`→`add_bond_harmonic`、
`add_harmonic_angle`→`add_angle_harmonic`、
`add_multi_harmonic_dihedral`→`add_dihedral_multi_harmonic`），md.py 的
`_PHASE_DISPATCH`/`_DRIVER_DISPATCH` 随之更新。无 alias、无 deprecation 过渡。

**诚实面收口：**

- `__init__.py` 导出 = MD 面（driver + 组件，canonical 名）+ PW 面
  （`pw_scf`/`PWResult`），docstring 说明当前能力边界与 roadmap（SCF/Minimize
  绑定为后续独立 spec）。
- `pw.py`：默认 wheel（`ATV_BUILD_PW=OFF`）下 `_cpp.pw_scf` 不存在，现状是
  调用时抛裸 `AttributeError`——改为显式守卫，抛带指引的 `RuntimeError`
  （"rebuild with ATV_BUILD_PW=ON"），消除最后一块半假可用面。
- `python/docs/getting_started.md` 现仍教已删除的 `av.CpuFrame` API，重写为
  molrs.Frame 输入 + canonical 组件名 + roadmap 注记。
- 新增 `python/tests/test_surface.py`：诚实面 smoke（canonical 名 MD LJ 短跑、
  pw_scf toy 按 `hasattr(_cpp, "pw_scf")` 条件跳过）+ 死面 fail-loud 断言
  （`import atomiverse.scf` → `ModuleNotFoundError`，`av.SCF` 等 →
  `AttributeError`）+ 树内死引用 grep 门禁（pytest 内实现，扫描 python/ 下
  `_cpp.LBFGS`、`_cpp.SCFConfig` 等禁用 token；模式串在测试内拼接构造，
  避免自匹配）。

## Files to create or modify

- python/tests/test_surface.py (new)
- python/src/atomiverse/scf.py（删除）
- python/src/atomiverse/optimizer.py（删除）
- python/src/bindings/bind_scf.cpp（删除）
- python/src/bindings/bind_optimizer.cpp（删除）
- python/src/bindings/bind.cpp（更新 :21-22 过期注释：SCF/optimizer 绑定为
  future work，指向 Out of scope）
- python/examples/01_h2_rhf.py（删除）
- python/examples/02_h2_uhf.py（删除）
- python/examples/03_uks.py（删除）
- python/examples/04_aspirin_opt.py（删除）
- python/examples/04_aspirin_bench.py（删除）
- python/examples/04_aspirin_opt_pyscf.py（删除）
- python/examples/single_point_bench.py（删除）
- python/examples/utils/__init__.py（删除）
- python/examples/utils/timing.py（删除）
- python/examples/pw_qe_compare.py（改写 frame 构造为 molrs.Frame）
- python/src/bindings/bind_md.cpp（canonical Python 类名 + 注册方法名）
- python/src/atomiverse/md.py（dispatch 表 + re-export + docstring 对齐）
- python/src/atomiverse/__init__.py（诚实导出面 + docstring/roadmap）
- python/src/atomiverse/pw.py（`_cpp.pw_scf` 缺失时的显式守卫）
- python/tests/test_md.py（canonical 名迁移）
- python/tests/conftest.py（canonical 名迁移）
- python/benchmarks/bench_md.py（canonical 名迁移）
- python/docs/getting_started.md（重写为当前诚实面 + roadmap）

## Tasks

- [ ] Write failing tests for the honest Python surface (python/tests/test_surface.py): canonical names importable (PairLJCut, BondHarmonic, AngleHarmonic, DihedralMultiHarmonic), legacy names and SCF/RHF/UHF/UKS/GeometryOptimizer/Frame absent, `import atomiverse.scf` / `atomiverse.optimizer` raise ModuleNotFoundError, MD LJ 10-step smoke on an 8-atom Ar frame, pw_scf toy smoke gated by `hasattr(_cpp, "pw_scf")`, and an in-tree dead-reference grep gate over python/
- [ ] Delete the dead SCF/optimizer facade: python/src/atomiverse/scf.py, python/src/atomiverse/optimizer.py, python/src/bindings/bind_scf.cpp, python/src/bindings/bind_optimizer.cpp; update the stale :21-22 comment in python/src/bindings/bind.cpp
- [ ] Delete dead examples and orphaned helpers: python/examples/{01_h2_rhf,02_h2_uhf,03_uks,04_aspirin_opt,04_aspirin_bench,04_aspirin_opt_pyscf,single_point_bench}.py and python/examples/utils/ (keep pseudopotentials/)
- [ ] Rewrite python/examples/pw_qe_compare.py frame construction from the removed `av.Frame` API to molrs.Frame/molrs.Block (keep the QE comparison flow and pseudopotentials/ paths)
- [ ] Rename Python-exposed MD classes to canonical names in python/src/bindings/bind_md.cpp (LJ→PairLJCut, HarmonicBond→BondHarmonic, HarmonicAngle→AngleHarmonic, MultiHarmonicDihedral→DihedralMultiHarmonic) and rename the PhaseBuilder registration methods to match (add_pair_lj_cut, add_bond_harmonic, add_angle_harmonic, add_dihedral_multi_harmonic)
- [ ] Update python/src/atomiverse/md.py (dispatch tables, re-exports, docstrings), python/src/atomiverse/__init__.py (honest MD+PW exports + roadmap docstring), and python/src/atomiverse/pw.py (explicit RuntimeError guard with an ATV_BUILD_PW=ON hint when `_cpp.pw_scf` is absent)
- [ ] Update surviving consumers to canonical names: python/tests/test_md.py, python/tests/conftest.py, python/benchmarks/bench_md.py
- [ ] Rewrite python/docs/getting_started.md against the current surface (molrs.Frame input, canonical component names, honest capability boundary + roadmap note)
- [ ] Run full check + test suite (editable wheel build via scikit-build-core, then pytest python/tests; C++ ctest gate unaffected)

## Testing strategy

- **Happy path**：默认配置（`ATV_BUILD_RHF=ON`、`ATV_BUILD_PW=OFF`）
  `pip install --no-build-isolation -ve python/` 构建成功；
  `pytest python/tests` 全绿——现有 test_md.py 全套（canonical 名迁移后）+
  test_surface.py 的 MD LJ 短跑（8 原子 Ar，10 步，坐标演化、能量/KE 有限，
  单位 eV/Å/fs）。
- **PW 面**：`ATV_BUILD_PW=ON` 构建时 pw_scf toy smoke 实跑
  （空 upf_path，`iterations >= 1`、能量有限 [Ry]）；`OFF` 时该测试按
  `hasattr(_cpp, "pw_scf")` 干净跳过，且 `av.pw_scf(...)` 抛带
  ATV_BUILD_PW 指引的 RuntimeError（fail-loud，非裸 AttributeError）。
- **Edge cases（防假可用回潮）**：`import atomiverse.scf` /
  `import atomiverse.optimizer` → ModuleNotFoundError；`av.SCF`、`av.RHF`、
  `av.GeometryOptimizer`、`av.Frame`、`av.LJ`、`av.HarmonicBond` →
  AttributeError；`_cpp` 上无 LBFGS/LBFGSStep/MinimizerPhase/SCFConfig/
  XCType/scf_solve/scf_gradient/build_integrals_sto3g 属性。
- **grep 门禁**：test_surface.py 内的树扫描测试——python/ 下（排除
  test_surface.py 自身）不得出现 `_cpp.LBFGS`、`_cpp.LBFGSStep`、
  `_cpp.MinimizerPhase`、`_cpp.SCFConfig`、`_cpp.SpinType`、
  `_cpp.MeanFieldType`、`_cpp.XCType`、`_cpp.scf_solve`、`_cpp.scf_gradient`、
  `_cpp.build_integrals_sto3g`，以及 `CpuFrame`、`av.SCF`/`av.RHF`/`av.UHF`/
  `av.UKS`/`av.GeometryOptimizer`/`av.Frame` 等死符号 token。
- **域验证**：不引入新物理，无新参考值；沿用既有 LJ-Ar 参数与 toy-PW
  "可运行、不物理"契约（见 Domain basis）。

## Out of scope

- **正经的 Minimize/SCF Python 绑定**：以正确形状（`cpu::minimizer::Minimize`
  driver + `LBFGS`/`CG` Optimizer 组件 + `cpu::scf::SCFEvaluator`/
  `build_integrals`/`rhf_solve` 家族）重绑 bind_optimizer/bind_scf 是后续
  独立 spec；本篇只拆除坏死面，不重建。
- 09-harness 及链内 02-07 的任何 C++ 侧改动；本篇不触碰
  include/、src/、tests/（C++）。
- `python/src/atomiverse/_device.py` molpy 设备注册占位（fail-loud 存根、
  非自动导入，不属假可用面；其 docstring 修辞留待 SCF 重绑 spec 一并清理）。
- CUDA/torch wheel 变体、PyPI 发布流程、根目录 docs/ 中的叙事文档。
- 新的 MD/PW 教程 example 序列（examples 收口后仅存 pw_qe_compare.py +
  pseudopotentials/；教程重建随 SCF 重绑 spec 走）。

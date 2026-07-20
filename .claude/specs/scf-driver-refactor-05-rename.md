---
title: SCF 驱动化重构 · 05 编译开关改名（ATV_BUILD_RHF→ATV_BUILD_HF，ATV_HAS_RHF→ATV_HAS_HF）
slug: scf-driver-refactor-05-rename
status: approved
chain: scf-driver-refactor
link: 05/08
created: 2026-07-20
---

# SCF 驱动化重构 · 05 编译开关改名

## Summary
纯机械改名：CMake 选项 `ATV_BUILD_RHF`→`ATV_BUILD_HF`、C++ 宏 `ATV_HAS_RHF`→`ATV_HAS_HF`，跨两后端全树同步（~67 个非本链-spec 文件：根 CMake 选项/门/define/GPW 依赖检查、CPU+CUDA 所有 `#ifdef`/`#define`、tests/examples CMake、tools/ci + slurm 脚本、python/pyproject、CLAUDE.md/AGENTS.md/README/docs 的选项表与构建配方）。不设按方法编译开关。无任何行为变化——同样的源、同样的测试，只是开关名从"RHF"（其实早已编译全部 RHF/UHF/RKS/UKS）正名为"HF"。这是决策⑤"单一 `ATV_BUILD_HF`、按文件拆分"的最后一块（拆分在 03a 已落）。

## Domain basis
纯机械改名，无物理、无数值路径变动；领域锚点是**行为保持**——改名后在新开关下构建+全测试与改名前逐位一致。开关旧名 `RHF` 是历史遗留（它一直门控着 RHF/UHF/RKS/UKS 全部高斯均场），`HF` 正名反映其真实语义：高斯基均场（HF + 杂化 KS）。沿用决策日志⑤。

## Design
两个 token 的一对一替换，分四组落地（每组内纯文本替换，不改逻辑）：

**1. 构建系统（权威定义处）：** 根 `CMakeLists.txt`——选项声明（约 65 行）、`ATV_HAS_RHF` 的 `target_compile_definitions`（约 380/427/431 行）、GPW 依赖检查（`ATV_BUILD_GPW` 需 `ATV_BUILD_RHF`→`ATV_BUILD_HF`）；`tests/CMakeLists.txt`、`examples/CMakeLists.txt`；`tests/architecture/module_dag_auto_enable.cmake`、`module_off_excludes_sources.cmake`（模块 DAG/排除规则里的开关名）。

**2. 源码 `#ifdef ATV_HAS_RHF`（两后端）：** CPU `include/atomiverse/cpu/observer/density.h`；CUDA `cuda/scf/{scf_solver,scf_evaluator,internal/rhf_density,internal/rhf_fock}.cuh` + `src/cuda/scf/{scf_solver,rhf_fock,scf_evaluator,rhf_density}.cu` + `cuda/observer/mulliken.cuh` + `cuda/recorder/density_recorder.cuh`。**只改宏名，不动 CUDA 逻辑**（CUDA SCF 保持 inner-solver）。

**3. 脚本 + 打包：** `tools/ci_cpu_unit.sh`（约 36 行 `-DATV_BUILD_RHF=ON`）、`tools/gpu_job.sh`（约 18 行）、`slurm/{job_cpu_full,job_cpu_gpw_test,job_cpu_xc_smoke}.sh`、`python/pyproject.toml`。

**4. 文档 + 示例/测试 prose（flagged docs）：** `CLAUDE.md`/`AGENTS.md` 选项表 + 构建配方、`README.md`、`.claude/architecture.md`、`docs/{user-guide/methods,user-guide/getting-started,tutorial/ch02-rhf,ch03-open-shell,ch05-kohn-sham-dft,ch09-gpw}.md`、`examples/README.md`、`examples/results/README.md`；示例/测试里出现的 `-DATV_BUILD_RHF` 构建注释与 `#ifdef`（`examples/{02_rhf.cpp,02_rhf.cu,05_geom_opt,06_aimd,09_gpw}.cpp`、`tests/{test_xc_reference,cpu/test_evaluate_driver,cpu/scf/test_uks_solver,cpu/scf/test_koopmans}.cpp` + CUDA 测试 `test_async_rhf_headers/test_device_integrals/test_rhf/test_driver_binding_kind`）。

**不改：** `.claude/specs/*`（本链自己的 spec/acceptance，`done` 时删除；grep 验收排除之）与 git 历史。

**Reuse decision：** 纯改名，无新符号、无 `make_*`、无新 `.sh` 门；不新增/删除任何开关（不设按方法开关）。所有改动都是同一对 token 的文本替换。

## Files to create or modify
- `CMakeLists.txt` — 选项/define/GPW 依赖检查改名（约 65/380/427/431 行）。
- `tests/CMakeLists.txt`、`examples/CMakeLists.txt`、`tests/architecture/module_dag_auto_enable.cmake`、`tests/architecture/module_off_excludes_sources.cmake` — 开关名。
- `include/atomiverse/cpu/observer/density.h` — `#ifdef ATV_HAS_RHF`→`ATV_HAS_HF`。
- `include/atomiverse/cuda/scf/{scf_solver,scf_evaluator}.cuh`、`cuda/scf/internal/{rhf_density,rhf_fock}.cuh`、`cuda/observer/mulliken.cuh`、`cuda/recorder/density_recorder.cuh` + `src/cuda/scf/{scf_solver,rhf_fock,scf_evaluator,rhf_density}.cu` — 宏名。
- `tools/ci_cpu_unit.sh`、`tools/gpu_job.sh`、`slurm/{job_cpu_full,job_cpu_gpw_test,job_cpu_xc_smoke}.sh`、`python/pyproject.toml` — 构建标志。
- `CLAUDE.md`、`AGENTS.md`、`README.md`、`.claude/architecture.md`、`docs/user-guide/{methods,getting-started}.md`、`docs/tutorial/{ch02-rhf,ch03-open-shell,ch05-kohn-sham-dft,ch09-gpw}.md`、`examples/README.md`、`examples/results/README.md` — 选项表/配方（docs）。
- `examples/{02_rhf.cpp,02_rhf.cu,05_geom_opt,06_aimd,09_gpw}.cpp` + `tests/{test_xc_reference.cpp,cpu/test_evaluate_driver.cpp,cpu/scf/test_uks_solver.cpp,cpu/scf/test_koopmans.cpp}` + CUDA 测试 `tests/cuda/{async/test_async_rhf_headers,scf/test_device_integrals,scf/test_rhf,test_driver_binding_kind}.cu` — 构建注释/`#ifdef`。

## Tasks
- [ ] Rename the CMake option, ATV_HAS_HF define, and GPW dependency checks in CMakeLists.txt, tests/CMakeLists.txt, examples/CMakeLists.txt, and tests/architecture/*.cmake
- [ ] Rename ATV_HAS_RHF → ATV_HAS_HF in all CPU and CUDA source #ifdef/#define sites (include/ + src/, both backends), touching macro names only
- [ ] Rename the -DATV_BUILD_RHF flag in tools/ci_cpu_unit.sh, tools/gpu_job.sh, slurm/{job_cpu_full,job_cpu_gpw_test,job_cpu_xc_smoke}.sh, and python/pyproject.toml
- [ ] Rename the flag in the doc/example/test option tables + build recipes (CLAUDE.md/AGENTS.md/README.md/.claude/architecture.md/docs/*/examples README + example/test build comments) — docs
- [ ] Run full check + test suite under -DATV_BUILD_HF=ON

## Testing strategy
- **改名后构建+全测试**：在 `-DATV_BUILD_HF=ON` 下 `cmake` + `build` + `ctest` 与改名前逐位一致（同样的目标编译、同样的测试通过）。
- **零残留**：evaluator grep `ATV_BUILD_RHF|ATV_HAS_RHF` 在 `include/ src/ tests/ examples/ tools/ slurm/ docs/ cmake/ python/` + `CMakeLists.txt`/`CLAUDE.md`/`AGENTS.md`/`README.md`/`.claude/architecture.md` 下**计数为零**（排除 git 历史与 `.claude/specs/*`）。
- **GPW 依赖**：`ATV_BUILD_GPW` 仍正确要求 `ATV_BUILD_HF`（+`ATV_BUILD_PW`）；开错组合仍报错。
- **CUDA**：宏改名不动 CUDA 逻辑；CUDA 构建在 GPU 节点验证（本环 CPU 登录节点只验证 CPU 全绿 + grep 清零）。
- **Domain validation**：无数值改动，行为保持即验证。
- **Regression example**：本环**不**含 `regressions/` 案例——纯改名；本链回归在 03b。runtime 关口是新开关下的完整 ctest 快套件。

## Out of scope
- Handler→Criterion 迁移 + 删 `LoopHandler` —— 06。
- 按方法编译开关（不设）。
- 本链 spec/acceptance 文件里的旧 token（`.claude/specs/*`，`done` 时删）。

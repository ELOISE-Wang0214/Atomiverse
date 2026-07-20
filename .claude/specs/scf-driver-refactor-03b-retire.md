---
title: SCF 驱动化重构 · 03b 退役旧机器 + 首个高斯 SCF 回归（删 SCFSolver/SCFConfig，梯度解耦，文献值 e2e）
slug: scf-driver-refactor-03b-retire
status: approved
chain: scf-driver-refactor
link: 03b/06
created: 2026-07-20
---

# SCF 驱动化重构 · 03b 退役旧机器 + 首个高斯 SCF 回归

## Summary
03a 已把 `RHF`/`UHF`/`RKS`/`UKS` 改成 `CPUDriver<scf::Phase>` 驱动、并对旧 `SCFSolver` 参考逐位对拍通过；本环删除这套已被证冗余的旧机器——`SCFSolver`（`internal/scf_solver.{h,cpp}`）、`ScfLoop`/`make_config`/`config_sp_` 残余、`SCFConfig`（`internal/scf_types.h`，保留 `SpinResolution`）、`HamiltonianBuilder` 的 SCFConfig 重载——并把最后几处 `SCFConfig` 消费方（解析梯度 `scf_gradient`/`scf_gradient_fd` + `integral_derivatives`、`SCFResultHandle.config`、`vib/hessian.h`+`evaluate.h` 文档示例、`delta_scf_reference.h`）迁到离散参数/成员 API。退役 `test_scf_solver_internals`/`_internal_isolation`（有意义的 helper 断言并入 02 的 stage 测试）、删除 `tests/cuda/scf_solver_isolation.sh` bash 门（存活的隔离关注改为验收 evaluator grep，**不设替代脚本**）。最后交付本链**唯一**的高斯 SCF 端到端回归 `regressions/scf/h2o_rhf_sto3g/`——库外消费方、对**文献值**核对、≤60s、`ATV_BUILD_REGRESSIONS` 门控、只经 `regressions` target 跑。整套快套件秒级全绿。

## Domain basis
本环两个领域锚点：
1. **行为保持**——删除/迁移都是机械改写，不改任何数值路径：解析梯度改为吃离散成员而非 `SCFConfig`（读的字段不变），故既有 `test_gradient*`/`test_gradient_vs_libcint` 力值不变；四方法能量在 03a 已逐位一致。参考实现在此才删，正因为 03a 已证毕（参考先于替换物存活）。
2. **文献值 e2e**——本链首个高斯 SCF 回归：H₂O/STO-3G 限制性 Hartree–Fock 总能量。引用值（实现期已按规格要求核实并修正：原 −74.9420 Ha/Szabo–Ostlund 配对系误引，实为 Crawford Project #3 拉伸几何值）：LibreTexts *Live Textbook of Physical Chemistry*（Peverati）§28.2；几何 R(OH)=0.95 Å、∠HOH=104.5°；**E_RHF/STO-3G(H₂O) = −74.961754 Hartree**，容差 **1×10⁻³ Ha**（端到端理智核对量级，非高精度基准）。精确数字 + 几何 + 引用逐字记入该案例 `REFERENCE.md`（遵 `regressions/README.md` 约定）；实现时对照所引版次核实该数位。梯度导数参考：Pulay, *Mol. Phys.* **17**, 197 (1969)；Szabo & Ostlund, App. C。

沿用整链决策日志（01/02/03a 已固化）。

## Design

**删除清单（不再被任何方法/消费方引用后）：** `SCFSolver`（`internal/scf_solver.{h,cpp}`）、`methods` 里的 `ScfLoop`/`make_config`/`publish_config_`/`config_sp_`（03a 仅为旧 handle 快照留存者）、`SCFConfig` 结构体（`internal/scf_types.h`）、`HamiltonianBuilder::prepare/step` 的 SCFConfig 重载（02 的离散入口存活）。

**`internal/scf_types.h` 保留什么：** 仅 `SpinResolution`（(n_electrons, multiplicity)→(n_alpha,n_beta) 解析器）+ 其必要 include；`SCFConfig` 结构体删除。文件从"SCFConfig + SpinResolution"缩为"SpinResolution"。`hamiltonian_builder.h`/`gradient.h`/`integral_derivatives.h` 改为直接 include `scf_result.h`（取 `MeanFieldType`/`SpinType`）+ `xc_functional.h`（取 `XCType`），不再经 `scf_types.h` 取 `SCFConfig`。

**根部 tombstone 决策：** `cpu/scf/scf_solver.h` 与 `cpu/scf/scf_types.h` 是 `#error` 墓碑（拦截旧公开 include）。**决定：两者仍保留为 `#error` 墓碑**（把误 include 旧路径重定向到 `methods.h` 比"file not found"友好），仅更新文案：`scf_solver.h` 去掉已删的 "Internal: cpu/scf/internal/scf_solver.h" 指引、改述"SCFSolver 已删除，SCF 现为 RHF/UHF/RKS/UKS 内的 `CPUDriver<scf::Phase>` 循环"；`scf_types.h` 改述"SCFConfig 已删除，参数是方法成员；SpinResolution 仍在 `cpu/scf/internal/scf_types.h`"。

**逐消费方：从 `config` 读什么 → 替换：**

| 消费方 | 从 `config` 读 | 替换 |
|---|---|---|
| `density_cube_writer.cpp:58` | `cfg.spin`（限制/非限制） | 新增 `SCFResult::spin`（方法后置填），writer 读 `result.spin`；去掉 null-config 检查 |
| `orbital_cube_writer` / `orbital_info_observer` | 无（仅 `result`/`integrals`/`ready`） | 不变 |
| `scf_gradient`（`gradient.cpp:481`，自由函数） | `method`/`spin`/`xc_type`/`n_alpha`/`n_beta`/`grid_level`/`n_radial`/`n_angular`/`hybrid_alpha_override`/`force_native_two_electron_deriv`/`xc_grad_project_translation`/`fixed_grid` | **类化为内部 `internal::ScfGradient` 类**（用户裁定无自由 helper 函数）：构造器吃这些离散方法成员，成员函数求梯度；方法私有持有，主体数值代码原样搬入不改 |
| `scf_gradient_fd`（`gradient.cpp:609`，FD 校验 oracle，自由函数） | 同上 + `charge`/`basis` | 并入 `ScfGradient` 的 FD 成员函数，或下沉进其消费测试（禁工厂、禁自由函数；由既有 FD 梯度测试保绿） |
| `integral_derivatives`（`.cpp:669`） | `force_native_two_electron_deriv` 等 | 其 `SCFConfig` 形参改离散参数（既有类/入口内改签名，不新增自由函数） |
| `SCFResultHandle.config`（`mean_field_method.h`） | 仅为 `density_cube_writer` 发布 | **删除该字段**；`config_sp_` + `make_config`/`publish_config_` 一并移除 |

**`SCFResultHandle` 收窄边界（与 04 分工，03a 已裁定）：** 本环删除 `config` 字段 + 迁移其唯一消费方；`integrals` 字段留到 **04** 收窄为方法自有具名 handle。

**分阶段 run/launch 债务不在本环结清：** `MeanFieldMethod` 仍是 `HF`/Component、仍有 `build`/`launch`——那由 **04** 删除。本环只删旧 SCF *机器*，不动方法的驱动/组件双身份。

**Reuse decision：**
- `reuse` 02 的 stage 组件 + 03a 的方法驱动 —— 回归与迁移后的梯度全走驱动化方法；不新增数值。
- `reuse` 现有梯度实现（`gradient.cpp`/`integral_derivatives.cpp`）—— 仅把 `SCFConfig` 形参换成离散成员，主体不动。
- `reuse` 现有回归 harness（`regressions/CMakeLists.txt`、`reg_*.cpp` 库外消费方模式、`REFERENCE.md` 约定）—— 新 `regressions/scf/h2o_rhf_sto3g` 仿 `reg_am1_water_force.cpp` 形状。
- `reuse` 02 `test_scf_stages` 覆盖 —— 吸收 `test_scf_solver_internals` 的 helper 断言。
- `new` —— 仅回归案例 + `SCFResult::spin` 字段。**无新自由函数、无 `make_*`、无新 `.sh` 门。**
- `delete`（非复用）：`SCFSolver`/`ScfLoop`/`SCFConfig`/`make_config`/`config_sp_` + HamiltonianBuilder SCFConfig 重载。

## Files to create or modify
- `include/atomiverse/cpu/scf/internal/gradient.h` + `src/cpu/scf/gradient.cpp` — `scf_gradient`/`scf_gradient_fd` 的 `SCFConfig` 形参改离散参数。
- `include/atomiverse/cpu/scf/internal/integral_derivatives.h` + `src/cpu/scf/integral_derivatives.cpp` — 同上。
- `include/atomiverse/cpu/scf/mean_field_method.h` (03a) — 删 `SCFResultHandle::config` + `config_sp_`。
- `src/cpu/scf/mean_field_method.cpp` (03a) — 删 `make_config`/`publish_config_`/`ScfLoop` 残余；后置填 `result.spin`；向 `scf_gradient` 传离散成员。
- `include/atomiverse/cpu/scf/scf_result.h` — `SCFResult` 增 `SpinType spin`。
- `src/cpu/recorder/density_cube_writer.cpp` — 读 `result.spin` 取代 `handle_.config`。
- `include/atomiverse/cpu/scf/internal/scf_solver.h` (delete) + `src/cpu/scf/scf_solver.cpp` (delete) — SCFSolver。
- `include/atomiverse/cpu/scf/internal/scf_types.h` — 删 `SCFConfig` 结构体（留 `SpinResolution` + 收窄 include）。
- `include/atomiverse/cpu/scf/internal/hamiltonian_builder.h` + `src/cpu/scf/internal/hamiltonian_builder.cpp` — 删 SCFConfig 重载（留 02 离散入口）。
- `include/atomiverse/cpu/scf/scf_solver.h` + `include/atomiverse/cpu/scf/scf_types.h` — 更新根部 `#error` 墓碑文案。
- `include/atomiverse/cpu/vib/hessian.h` + `include/atomiverse/cpu/driver/evaluate.h` — `@code` 文档示例迁到成员 API。
- `tests/cpu/scf/delta_scf_reference.h` — `SCFConfig`+`SCFSolver` → `RHF`/`UHF` 成员 API。
- `tests/cpu/scf/stages/test_scf_stages.cpp` (02) — 吸收 `test_scf_solver_internals` 的有意义 helper 断言。
- `tests/cpu/scf/test_scf_solver_internals.cpp` (delete) + `tests/cpu/scf/test_scf_solver_internal_isolation.cpp` (delete)。
- `tests/cuda/scf_solver_isolation.sh` (delete) + `tests/cuda/observer/rdf_molrs_migration.sh` — 扫掉第 9 行 "Mirrors ... scf_solver_isolation.sh" 悬挂注释。
- `CMakeLists.txt` (root) — 去掉 `src/cpu/scf/scf_solver.cpp`。
- `tests/CMakeLists.txt` — 去掉 `test_cpu_scf_solver_internals`、`test_cpu_scf_solver_internal_isolation`（:482-484）、`test_cuda_scf_solver_isolation`（:671-672）注册。
- `regressions/scf/h2o_rhf_sto3g/reg_h2o_rhf_sto3g.cpp` (new) + `regressions/scf/h2o_rhf_sto3g/REFERENCE.md` (new) + `regressions/CMakeLists.txt` (modify) — 注册该案例。

## Tasks
- [x] Decouple the gradient path from SCFConfig by classifying it as internal::ScfGradient (ctor takes the discrete method parameters; FD oracle becomes a member or relocates into its consuming test; no free helper functions) with integral_derivatives re-signatured to discrete parameters (existing gradient tests must keep identical forces)
- [x] Remove SCFResultHandle::config + config_sp_ + make_config/publish_config_ remnants; add SCFResult::spin (set in the method epilogue); migrate density_cube_writer to result.spin
- [x] Migrate the vib/hessian.h + evaluate.h @code doc examples and tests/cpu/scf/delta_scf_reference.h to the member API (no SCFConfig token)
- [x] Delete SCFSolver (internal/scf_solver.{h,cpp}), the SCFConfig struct (internal/scf_types.h, keeping SpinResolution), and HamiltonianBuilder's SCFConfig overloads; update the root scf_solver.h/scf_types.h tombstone text
- [x] Retire test_scf_solver_internals + test_scf_solver_internal_isolation, moving their meaningful helper assertions onto the 02 stage test; delete tests/cuda/scf_solver_isolation.sh + its ctest registration and sweep the stale rdf_molrs_migration.sh comment
- [x] Update the root CMakeLists (drop scf_solver.cpp) and tests/CMakeLists (drop the three retired registrations)
- [x] Add the regression regressions/scf/h2o_rhf_sto3g (reg_h2o_rhf_sto3g.cpp + REFERENCE.md) and register it in regressions/CMakeLists.txt
- [ ] Run full check + test suite

## Testing strategy
- **Regression（端到端，文献值，ATV_BUILD_REGRESSIONS，只经 `regressions` target）**：`reg_h2o_rhf_sto3g.cpp` 作**库外消费方**（仅 include 公开 `cpu/scf/methods.h`+`frame.h`，链 atomiverse）：构 `cpu::scf::RHF{basis=STO3G}`、`run(frame)`、断言 `last_result().energy` 与 `REFERENCE.md` 引用值 −74.961754 Ha 差 ≤1e-3 Ha，≤60s；绝不进默认 ctest。
- **梯度解耦不变性**：既有 `test_gradient*` / `test_gradient_vs_libcint` 以迁移后离散签名跑出**逐位相同**力值（证 `SCFConfig`→离散只是形参改写）。`scf_gradient_fd`（FD 校验 oracle）以离散参数或 SCF-runner 回调重写，由既有 FD 梯度测试保绿。
- **Handle 收窄**：`test_density_cube` / `test_orbital_info_observer` / example 16 编译+通过；`density_cube_writer` 经 `result.spin` 得到与旧 `cfg.spin` 相同的限制/非限制分支。
- **迁移测试**：`test_delta_scf`（经迁移后的 `delta_scf_reference.h`）、`test_koopmans`/`test_level_shift`/`test_accelerators`（本已成员 API）全绿、值不变。
- **删除后编译**：全树无 `SCFSolver`/`ScfLoop`/`SCFConfig`/`make_config` 活引用（除根墓碑 `#error` 串、CUDA 自有 `cuda::scf::SCFSolver` 及交叉引用注释）。
- **隔离关注存续**：退役的 `_internal_isolation` 关注改为验收 grep（见 acceptance ac-008），不落地脚本。
- **Domain validation**：回归的文献值核对即本环领域验证。
- **RED/GREEN**：本环为删除+迁移相，无新增 impl 符号需配对 RED 单测（梯度迁移由既有梯度测试覆盖）；唯一新增制品是端到端回归（交付时 runtime 验证）。
- **Regression example**：本环即整链**唯一** `regressions/scf/` 案例的交付点。

## Out of scope
- 删 `build`/`launch`（结束 run/launch 共存）、`HF` 收窄纯 tag、原生嵌套、删 `sub_driver.h`、`SCFResultHandle.integrals` 收窄为方法具名 handle + 其消费方迁移 —— 属 04。
- `ATV_BUILD_RHF→ATV_BUILD_HF` 改名 —— 属 05。
- 既有 `LoopHandler` 全树迁移到 `Criterion` 与删除 —— 属 06。
- CUDA `cuda::scf::SCFSolver`/`SCFConfig`-注释、其他文件里零散的 `@see cpu::scf::SCFSolver` 注释残迹 —— 随各自文件下次触碰时清扫；本环只迁 `vib/hessian.h`+`evaluate.h` 两处会误导用户的 `@code` 示例。
- 高精度基准（bench 仓，`scientific`）—— 本环只交付库内 e2e 理智回归（`runtime`）。

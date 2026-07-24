---
slug: scf-review-04-meanfield-dedup
title: "均场方法 FD Hessian 与梯度结果转换去重（scf-review-04）"
status: done
created: 2026-07-24
domain: scf / mean-field
queue: scf-review
depends_on: []
grilled: true
---

# 均场方法 FD Hessian 与梯度结果转换去重（scf-review-04）

## Summary

将 `cpu/scf/` 下三处于 2026-07-23/23/24 三次独立复制粘贴落地的重复代码收敛为共享抽象，行为完全不变（byte-identical replay）。其一，`rks.cpp`(43-167)、`uks.cpp`(46-183)、`rohf.cpp`(46-162) 中逐字节相同的"中心有限差分核 Hessian 驱动"上提为 `MeanFieldMethod<Derived>` 上的一个受保护 CRTP 模板方法，由各方法传入 `make_scratch` 仿函数与可选的位移前后状态钩子（UKS 的自旋密度初猜捕获）参数化，三个方法各自只保留约 5-10 行特化。其二，`rhf.cpp`(129-156) 与 `uhf.cpp`(45-72) 中仅类名字符串不同的 Hessian 前置守卫块（球谐基/ECP/未收敛抛出）下沉为同一 CRTP 基类上的私有辅助函数。其三，`gradient.cpp`、`internal/mp2_gradient.cpp`、`internal/ccsd_gradient.cpp` 三份以不同名字（`make_result`/`make_force_result`，均违反工厂函数禁令）实现的 Ha/Bohr→eV/Å 单位转换+符号翻转胶水，统一为 `GradientResult` 上的一个带标签消歧的转换构造函数。用户可见的公开 API 与所有已 pin 的 PySCF 数值一律不变。

## Domain basis

核 Hessian 采用解析力的中心有限差分（与离线 PySCF 2.14.0 pin 同一算法类）：

- 位移步长 **h = 1e-4 Bohr**（`h_ang = h · kBohrInAngstrom`），与 PySCF pin 严格一致。
- 梯度关系 `g = −F`；Hessian 元 `H_{row,col} = −∂F_row/∂R_col`，由 `(F(+h) − F(−h)) / (2h)` 中心差分得到，随后对称化 `½(H + Hᵀ)` 以吸收 FD 噪声。
- 单位常量（Atomiverse CODATA）：`kHartreeInEv = 27.211386245988`，`kBohrInAngstrom = 0.529177210903`；`f_to_ha_bohr = kBohrInAngstrom / kHartreeInEv`（eV/Å → Ha/Bohr），`hess_to_ev_ang2 = kHartreeInEv / kBohrInAngstrom²`（Ha/Bohr² → eV/Å²）。
- 梯度结果转换：力 `F[eV/Å] = −convert(g[Ha/Bohr], HartreePerBohr, eVPerAngstrom)`，能量 `E[eV] = convert(E[Ha], Hartree, eV)`；`units::convert` 为唯一原语。

**科学不变量（编码进验收）**：本次上提"逐字节等价"当且仅当以下四项被冻结——(1) 步长 h = 1e-4 Bohr；(2) 每个位移使用全新求解器语义（`make_scratch` 每次 ±h 构造新的 `Derived`，无 DIIS/MO 暖启动串扰）；(3) 位移遍历顺序（原子为外层、轴为内层；`col = 3A+d`，`row = 3B+e`；先 +h 后 −h）；(4) BLAS 线程配置不变、位移循环保持串行（不得引入会改变浮点归约次序的新并行）。规范正文冻结全部四项，验收门 = 现有 FD-Hessian 与梯度回归/单测复现完全相同的数值（byte-identical replay）。

参考：B. G. Johnson, J. A. Pople, P. M. W. Gill, M. J. Frisch, *Chem. Phys. Lett.* **220**, 377 (1994)（解析 DFT 二阶导，目标闭式）；J. A. Pople 等, *Int. J. Quantum Chem. Symp.* **13**, 225 (1979)（MP2 梯度）；J. Gauss & J. F. Stanton, *J. Chem. Phys.* **104**, 2574 (1996)（CCSD 梯度）；PySCF 2.14.0（`hessian.rks` / `grad.mp2` / `grad.ccsd` 离线 oracle）。

## Design

**实体与去重决策（三项，全部行为保持）。**

1. **FD Hessian 模板（Item 3）** —— 在 `MeanFieldMethod<Derived>`（`mean_field_method.h:99`，已有受保护辅助群 `validate_components`/`resolve_occupation`/`set_predictor_`/`integral_workspace_` 于 :604+）新增一个受保护成员模板方法 `central_difference_nuclear_hessian(Frame& frame, const std::function<GradientResult(Frame&)>& forces_at)`。模板独占"共享骨架"：`r0` 几何快照、`coord()` span 存取 lambda、原子×轴 ±h 位移双重循环、`(Fp−Fm)/(2h)` 差分与 `f_to_ha_bohr`/`hess_to_ev_ang2` 单位换算、几何还原、`½(H+Hᵀ)` 对称化、组装 `vib::HessianResult`。方法专属部分（`make_scratch` 的字段拷贝、以及 UKS 的自旋密度初猜捕获）由调用方以 `forces_at` 仿函数携带——UKS 在仿函数内以 `&ref_guess` 调用 `forces`，RKS/ROHF 以默认 `nullptr`。`RKS::hessian`/`UKS::hessian`/`ROHF::hessian` 各自收缩为约 5-10 行：RI 前置守卫（保留内联，见下）、参考 `run()`+收敛检查、（UKS）捕获参考自旋密度、构造 `forces_at`、调用模板、还原后再 `run()` 重建参考波函数。RHF/UHF 使用解析 CPHF Hessian（`rhf.h:372`、`uhf.h:61`），不走 FD 模板。该模板体定义于 `mean_field_method.cpp`，由该文件已有的 `RHF/UHF/RKS/UKS/ROHF` 显式实例化覆盖。
   - *行为说明*：位移求解失败（不收敛）的异常文本采用方法名标记的统一措辞——该属于已抛出的错误路径，不被任何回归/单测断言，成功路径数值 byte-identical 不受影响。若需保留逐字节相同的"+h/−h atom A axis d"文本，实现方可将 `(A, d, sign)` 透传进 `forces_at`；本规范不强制。

2. **Hessian 前置守卫辅助（Item 3b）** —— 在同一 CRTP 基类新增私有 `check_hessian_preconditions(std::string_view method_name, const IntegralPack& pack) const`，集中 `rhf.cpp:129-156` / `uhf.cpp:45-72` 中仅类名不同的三段抛出：未收敛（`converged()`/`stop_reason()`）、球谐基 pack（`pack.spherical()`）、ECP pack（`pack.has_ecp`），方法名内插进消息以保持文本一致。`RHF::hessian`/`UHF::hessian` 在参考 `run()` 与取得 `integral_workspace_()` 之后调用之。RI 守卫（`ri != RI::None`）保留为各方法 `run()` 之前的内联单行，以严格保留"建 SCF 前快速失败"的既有顺序与逐字节消息（RI 消息本就按类名不同，且属唯一 3 行差异）——不并入辅助函数，避免任何执行顺序变更。

3. **梯度结果转换构造函数（Item 7）** —— `GradientResult`（`gradient_result.h:15`，当前为纯聚合 `{xt::xtensor<double,2> forces [eV/Å]; double energy [eV];}`）新增：`struct FromHartreeBohr {};` 标签（沿用 `SlaterDeterminant::Generalized{}` 于 `slater_determinant.h:59` 的先例），与转换构造 `GradientResult(FromHartreeBohr, const xt::xtensor<double,2>& grad_ha_bohr, double energy_ha)`——内部以 `units::convert` 为原语做符号翻转+单位换算，`natm` 由 `grad_ha_bohr.shape()` 推得。因用户声明构造函数会移除聚合初始化，需同时补回 `GradientResult() = default;`（现有 `GradientResult out;` 默认构造点：`gradient.cpp`、`mp2_gradient.cpp`、`ccsd_gradient.cpp`、`mean_field_method.cpp:1099`、`test_scf_type_homes.cpp:44` 均为默认构造+成员赋值，需要默认构造函数）以及一个显式的逐成员 eV 构造 `GradientResult(xt::xtensor<double,2> forces_ev, double energy_ev)`，以保住文档 `docs/tutorial/methods/ccsd-grad-h2o-sto3g.md:17` 所述 `GradientResult{forces, energy}` 花括号语法。随后删除三处文件局部辅助 `make_result`（`gradient.cpp:84-98`，`static`，兼违工厂禁令）与 `make_force_result`（`mp2_gradient.cpp:51-62`、`ccsd_gradient.cpp:35-46`），调用点改为 `GradientResult{FromHartreeBohr{}, grad, energy_ha}`。

**Reuse decision（消解 librarian 建议）：**
- `reuse units::convert` —— 转换构造函数以其为唯一换算原语（房内先例，直接复用）。
- `generalize GradientResult`（`gradient_result.h`，标签消歧构造）—— 采纳 librarian 决策，将三份 `make_*` 胶水提升为该类型上的一个标签构造函数以同时服务三个调用方（对应"Generalize `make_*` glue"任务，非并列的"Implement"任务）。
- `reuse SlaterDeterminant::Generalized{}` 标签先例 —— 新标签 `FromHartreeBohr{}` 的命名/构造/风格沿用之，使新代码与既有代码同读。
- `new — central_difference_nuclear_hessian`：FD 骨架此前无任何可复用宿主（三份复制粘贴），新受保护模板方法是最贴近的"既有受保护辅助群"模式（`mean_field_method.h:604+`）；命名/风格沿用之。
- `new — check_hessian_preconditions`：同上，无既有宿主，采纳 librarian 的"同 CRTP 基类私有辅助"决策。

**Out-of-scope 深层替代（显式推迟）：** 将三个 FD Hessian 改为委托 `cpu/vibration` 的 `vib::Hessian` FD 引擎——见 Out of scope 三条摩擦点。

## Files to create or modify

- `include/atomiverse/cpu/scf/mean_field_method.h` — 声明 `central_difference_nuclear_hessian`（受保护）与 `check_hessian_preconditions`（私有）
- `src/cpu/scf/mean_field_method.cpp` — 定义二者（由已有显式实例化覆盖）
- `src/cpu/scf/rks.cpp` — `hessian()` 委托模板
- `src/cpu/scf/uks.cpp` — `hessian()` 委托模板，保留自旋密度初猜捕获于调用方
- `src/cpu/scf/rohf.cpp` — `hessian()` 委托模板
- `src/cpu/scf/rhf.cpp` — 调用 `check_hessian_preconditions`
- `src/cpu/scf/uhf.cpp` — 调用 `check_hessian_preconditions`
- `include/atomiverse/cpu/scf/gradient_result.h` — `FromHartreeBohr` 标签 + 转换/默认/逐成员构造
- `src/cpu/scf/gradient.cpp` — 删除 `make_result`，改用标签构造
- `src/cpu/scf/internal/mp2_gradient.cpp` — 删除 `make_force_result`，改用标签构造
- `src/cpu/scf/internal/ccsd_gradient.cpp` — 删除 `make_force_result`，改用标签构造
- `tests/cpu/scf/test_fd_hessian_template.cpp` (new) — FD 模板单测（解析二次型 mock 力场）
- `tests/cpu/scf/test_gradient_result.cpp` (new) — 转换构造函数单测
- `tests/architecture/scf_dedup_scrub.sh` (new) — grep 门
- `tests/CMakeLists.txt` — 注册两个新单测与 scrub（`atv_add_cpp_test` / `add_test` 于 HF 门与 scrub 段）

## Tasks

- [x] Write failing unit tests for the central-difference nuclear-Hessian template in tests/cpu/scf/test_fd_hessian_template.cpp (analytic quadratic mock Derived, no SCF; register via atv_add_cpp_test in tests/CMakeLists.txt)
- [x] Implement central_difference_nuclear_hessian as a protected member template on MeanFieldMethod<Derived> (include/atomiverse/cpu/scf/mean_field_method.h, src/cpu/scf/mean_field_method.cpp), freezing h=1e-4 Bohr, fresh-scratch-per-displacement, atom-major/axis-minor ordering, serial loop
- [x] Refactor RKS/UKS/ROHF hessian() to delegate to the template (src/cpu/scf/rks.cpp, uks.cpp, rohf.cpp), keeping the UKS spin-density SCFInitialGuess capture in its caller functor
- [x] Implement check_hessian_preconditions on MeanFieldMethod<Derived> and call it from RHF/UHF hessian() (include/atomiverse/cpu/scf/mean_field_method.h, src/cpu/scf/mean_field_method.cpp, rhf.cpp, uhf.cpp); keep the RI guard inline pre-run
- [x] Write failing unit test for the GradientResult(FromHartreeBohr{}, …) converting ctor in tests/cpu/scf/test_gradient_result.cpp (register via atv_add_cpp_test in tests/CMakeLists.txt)
- [x] Generalize the Ha/Bohr→eV/Å gradient glue into GradientResult: add FromHartreeBohr tag + converting/default/member-wise ctors (include/atomiverse/cpu/scf/gradient_result.h) and delete make_result/make_force_result, rewiring gradient.cpp, internal/mp2_gradient.cpp, internal/ccsd_gradient.cpp
- [x] Add grep-gate scrub tests/architecture/scf_dedup_scrub.sh (no FD displacement loop in rks/uks/rohf.cpp; no make_result/make_force_result anywhere under src/cpu/scf) and register it in tests/CMakeLists.txt
- [x] Add docstrings per doxygen with units for central_difference_nuclear_hessian, check_hessian_preconditions, and the new GradientResult ctors
- [x] Verify the existing FD-Hessian and MP2/CCSD-gradient regressions reproduce byte-identical values (regressions/scf/rks-hessian-h2o-sto3g, uks-hessian-o2-sto3g, rohf-hessian-o2, uhf-hessian-o2, mp2-grad-h2o-sto3g, ump2-grad-o2-sto3g, ccsd-grad-h2o-sto3g)
- [x] Run full check + test suite

## Testing strategy

单测（仅 `tests/`，镜像源布局，单函数/方法，秒级；无完整 SCF）：

- `tests/cpu/scf/test_fd_hessian_template.cpp → TestFdHessianTemplate`（目标：`MeanFieldMethod<Derived>::central_difference_nuclear_hessian`）
  - happy path：以一个解析二次型 mock 力场（`F(R) = −K·(R − R0)`，`K` 为硬编码常量刚度矩阵）作为 `forces_at` 仿函数，验证 FD 中心差分复现该常量 Hessian `K`（换算至 eV/Å²），硬编码期望值，容差取 FD 截断+浮点级 1e-8。
  - edge cases：单原子（n3=3）；已对称输入验证 `½(H+Hᵀ)` 幂等；验证遍历顺序为 `col=3A+d` 外层原子/内层轴，且退出时几何还原到 `r0`（span 逐元素相等）。
  - 冻结不变量：断言步长常量为 1e-4 Bohr 对应的位移（对二次型，H 与 h 无关，用两个不同 h 得到相同 K 佐证中心差分正确性）。
- `tests/cpu/scf/test_gradient_result.cpp → TestGradientResult`（目标：`GradientResult(FromHartreeBohr, …)`）
  - happy path：输入 `grad_ha_bohr`、`energy_ha` 硬编码单元素，断言 `forces = −convert(...)`、`energy = convert(...)` 命中硬编码 eV/Å 与 eV 字面量（用 `kHartreeInEv=27.211386245988`、`kBohrInAngstrom=0.529177210903` 手算）；符号翻转到位。
  - edge cases：默认构造 `GradientResult{}`（`energy==0`，`forces` 空）；逐成员构造 `GradientResult{forces_ev, energy_ev}` 花括号语法保持可编译（保住文档契约）。

回归例（本 spec 不新增 `regressions/<slug>`）：本次为行为保持的去重，涉及的每个特性均已各自 ship 独立回归；按 Atomiverse"一特性一回归/不重复覆盖"，绑定的 runtime 门是现有回归以 **byte-identical replay** 复现其已 pin 的 PySCF 2.14.0 字面值：`regressions/scf/rks-hessian-h2o-sto3g`（如 `H_xx(0,0)=−9.632120712096 eV/Å²`，门 2e-2 eV/Å²）、`uks-hessian-o2-sto3g`、`rohf-hessian-o2`、`uhf-hessian-o2`（解析 Hessian，走前置守卫）、`mp2-grad-h2o-sto3g`（O `∂E/∂z=−9.528798092138e-02 Ha/Bohr`，门 1e-5）、`ump2-grad-o2-sto3g`、`ccsd-grad-h2o-sto3g`（O `∂E/∂z=−1.106969560091e-01 Ha/Bohr`，门 1e-5）。全部为库外消费者、字面值硬编码、无运行时第三方、≤60 s。

grep 门（`tests/architecture/scf_dedup_scrub.sh`，纯 grep，秒级，全配置运行）：`rks.cpp`/`uks.cpp`/`rohf.cpp` 不再含 FD 位移循环哨兵标记（`hess_to_ev_ang2`/`f_to_ha_bohr`/`(2.0 * h_bohr)`）；`src/cpu/scf/` 任何位置不再含 `make_result`/`make_force_result`。

## Out of scope

- 将三个 FD Hessian 委托给 `cpu/vibration` 的 `vib::Hessian` FD 引擎（更深层的统一）——显式推迟，三点摩擦：(1) `vib::Hessian` 固定 `kHessianStepAngstrom = 1e-3 Å`，与本路径 pin 的 1e-4 Bohr PySCF 平价步长冲突；(2) `vib::Hessian` 走组件 re-launch，会引入暖启动状态串扰，而本路径要求每位移全新方法实例；(3) 单位链为 eV/Å²，而本路径内部以 Ha/Bohr 差分再换算。三者任一未解前不做该合并。
- RHF/UHF 解析 CPHF Hessian 本体（不走 FD 模板，本 spec 不触及其数值）。
- 全解析 CPKS/UCPKS（RKS/UKS 闭式二阶导）——既有 follow-up，不在本次。
- RI/球谐基/ECP 的 Hessian 支持（当前仍按既有消息抛出，本 spec 仅去重守卫，不新增能力）。
- 位移失败路径异常文本的逐字节保真（错误路径，无门断言；如需保真由实现方透传位移标识）。

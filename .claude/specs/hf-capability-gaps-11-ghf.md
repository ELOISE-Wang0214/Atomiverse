---
title: GHF — 广义 Hartree-Fock(实二分量旋量)
slug: hf-capability-gaps-11-ghf
status: approved
chain: hf-capability-gaps
link: 11/21
external-depends: scf-slater-determinant-03-unbag
created: 2026-07-20
---

# GHF — 广义 Hartree-Fock(实二分量旋量)

## Summary

在 `cpu::scf` 中新增具体方法 `GHF`(CRTP 于 `MeanFieldMethod<GHF>`),实现实代数广义 Hartree-Fock:旋量分子轨道 φ_i = (φ_i^α; φ_i^β)(每个分量为完整 AO 展开),密度与 Fock 矩阵为 2n×2n(2×2 自旋分块)。这是本链中唯一强制推广行列式/密度存储的环节:两条独立自旋通道无法表示自旋非对角块 P^αβ,故必须先推广 `SlaterDeterminant` 为旋量形式——**该推广是本环节的架构决策,须与并行的 scf-slater-determinant 链负责人共同定案(倾向扩展现有类型而非分叉平行类型),并以 scf-slater-determinant-03-unbag 落地为前置**。用户可见成果:无自旋阻挫体系上 GHF 塌缩回 UHF 能量(正确性哨兵),三角 H3 等自旋阻挫体系上 GHF 通过非共线自旋降到 UHF 之下,并与外部钉定的 PySCF `scf.GHF` 一致。本环节仅实代数:不含复系数、不含自旋轨道耦合。依赖链内 10-rohf(教程约束阶梯叙事以 ROHF 先行)。

## Domain basis

- 旋量 MO:C 为 2n×N_mo 实矩阵;P = C_occ C_occᵀ(2n×2n),自旋分块 P = [[P^αα, P^αβ], [P^βα, P^ββ]]。
- Fock(实 HF):F^σσ′ = δ_σσ′ (h + J[P^αα+P^ββ]) − K[P^σσ′];Coulomb 只见电荷密度,交换作用于全部四块。E = ½Tr[P(H+F)] + E_nn,H = h⊗1₂;DIIS 误差用 S_spinor = S⊗1₂ 的 FPS−SPF。内部 Hartree,写回 Frame 为 eV(`core/units.h`)。
- 塌缩定理(作为测试编码):无自旋轨道项时所有 UHF 解都是 GHF 驻点;自旋无阻挫体系 GHF 极小值共线 ⇒ E_GHF = E_UHF。自旋阻挫(正三角 H3)时非共线解 E_GHF < E_UHF。初猜必须显式破坏共线性(α/β 块间小角度旋转混合),否则 DIIS 收敛到共线鞍点——此陷阱须由测试固定。
- 文献:Fukutome, Int. J. Quantum Chem. 20, 955 (1981), DOI 10.1002/qua.560200502(八类破缺对称分类,书目性引用);Helgaker–Jørgensen–Olsen, *Molecular Electronic-Structure Theory*, §10(一般自旋轨道);Seeger & Pople, J. Chem. Phys. 66, 3045 (1977), DOI 10.1063/1.434318(不稳定性关联,链 08)。
- Oracle:三角 H3/STO-3G,E_GHF ≤ E_UHF 且与外部运行 PySCF `scf.GHF` 钉定值符合到 1e-8 Ha;塌缩例 |E_GHF − E_UHF| ≤ 1e-8 Ha(PySCF 为外部 oracle,仅钉数值,无构建/测试期依赖)。

## Design

- `SlaterDeterminant`(`slater_determinant.h`)增加**广义(旋量)构造形态**:单通道 {C(2n×N_mo), ε, f, N},`density()` 返回 2n×2n,新增自旋分块访问器(P^αα/P^αβ/P^βα/P^ββ)与 `generalized()` 判别;`validate(S)` 用 S⊗1₂。既有受限/非受限语义与构造器不变(与 03-unbag 负责人协调签名,不分叉平行类型)。
- `SpinType` 增加 `Generalized`(`scf_result.h`);`SCFResult` 复用 α 通道字段承载旋量量(维度 2n),β 字段空。
- `internal::DensityBuilder` 增加旋量路径(2n×2n P + 分块抽取);`DensityStage` / `FockStage` 增加 generalized 分支(J 只由 P^αα+P^ββ 收缩,K 逐块收缩,AO 积分维度仍为 n)。`MeanFieldMethod` 既有成员矩阵(F_alpha_/C_alpha_/P_alpha_ 等)以 2n 维复用 α 通道,不新增状态袋。
- `GHF`:成员 `int n_electrons`(run 前必须设定)、`double guess_spin_rotation`(弧度,默认非零小角,破坏共线初猜;参数皆为成员,无 Config 袋、无工厂函数);`resolve_occupation` 置旋量占据。注册进 `methods.h` 伞头,显式实例化 `template class MeanFieldMethod<GHF>;`。本环节不做解析梯度(`compute_force` 路径对 GHF 抛出)。
- 生命周期/所有权同现有四方法:积分与 AO 工作区为方法私有,Frame 为唯一 I/O 通货。
- 最近模式参照:`uhf.h`/`uhf.cpp`(新公共方法接入 CRTP 底座的既有范例)。

## Files to create or modify

- include/atomiverse/cpu/scf/ghf.h (new)
- src/cpu/scf/ghf.cpp (new)
- include/atomiverse/cpu/scf/slater_determinant.h
- include/atomiverse/cpu/scf/scf_result.h
- include/atomiverse/cpu/scf/internal/density_builder.h
- include/atomiverse/cpu/scf/stages/density_stage.h
- include/atomiverse/cpu/scf/stages/fock_stage.h
- include/atomiverse/cpu/scf/methods.h
- src/cpu/scf/mean_field_method.cpp
- CMakeLists.txt
- tests/cpu/scf/test_ghf.cpp (new)
- tests/CMakeLists.txt
- regressions/scf/h3_ghf_sto3g/reg_h3_ghf_sto3g.cpp (new)
- regressions/scf/h3_ghf_sto3g/REFERENCE.md (new)
- regressions/CMakeLists.txt
- docs/tutorial/ch03-open-shell.md

## Tasks

- [ ] Write failing tests for GHF(tests/cpu/scf/test_ghf.cpp (new),注册于 tests/CMakeLists.txt):旋量 SlaterDeterminant 不变量(Tr(P·S⊗1₂)=N、幂等、分块访问器);塌缩测试 H2/STO-3G |E_GHF−E_UHF| ≤ 1e-8 Ha;非共线测试正三角 H3/STO-3G E_GHF ≤ E_UHF − 1e-6 Ha 且与钉定 PySCF 值符合 1e-8 Ha;共线初猜(guess_spin_rotation=0)停留 UHF 能量的陷阱测试。全部秒级。
- [ ] Implement 旋量(generalized)形态 in include/atomiverse/cpu/scf/slater_determinant.h 与 `SpinType::Generalized` in include/atomiverse/cpu/scf/scf_result.h(与 scf-slater-determinant-03-unbag 负责人协调定案后落地,扩展而非平行类型)
- [ ] Implement 旋量密度路径 in include/atomiverse/cpu/scf/internal/density_builder.h 与 include/atomiverse/cpu/scf/stages/density_stage.h(2n×2n P + 自旋分块抽取)
- [ ] Implement 旋量 Fock 组装 in include/atomiverse/cpu/scf/stages/fock_stage.h(J 由 P^αα+P^ββ、K 逐块;S⊗1₂ DIIS 误差)
- [ ] Implement `GHF` in include/atomiverse/cpu/scf/ghf.h + src/cpu/scf/ghf.cpp(成员 n_electrons、guess_spin_rotation;梯度路径抛出),注册 include/atomiverse/cpu/scf/methods.h,显式实例化于 src/cpu/scf/mean_field_method.cpp,源文件加入 CMakeLists.txt
- [ ] Add Doxygen(@brief/@param 含单位/@return)覆盖全部新公共符号(ghf.h、slater_determinant.h 新访问器、scf_result.h 枚举值)
- [ ] Implement 回归用例 regressions/scf/h3_ghf_sto3g/reg_h3_ghf_sto3g.cpp (new)(库外消费者、≤60 s、注册于 regressions/CMakeLists.txt),外部运行 PySCF 钉定 E_GHF/E_UHF 数值并连同脚本与引用写入 regressions/scf/h3_ghf_sto3g/REFERENCE.md (new)
- [ ] Write 教程节 "GHF 与自旋对称破缺" in docs/tutorial/ch03-open-shell.md(RHF→UHF→ROHF→GHF 约束移除阶梯、如实提及 Fukutome 八类、非共线自旋图像、三角 H3 上手示例 + molvis 自旋矢量渲染、章末 PySCF 对照表标量行)
- [ ] Run full check + test suite(含 tests/architecture/no_factory_scrub.sh 与 test_stage_typing_scrub)

## Testing strategy

- 快乐路径:H2/STO-3G GHF 收敛且能量 = UHF(塌缩);H3 非共线解低于 UHF 并匹配 PySCF。
- 边界:n_electrons 未设定 → 抛出;guess_spin_rotation=0 的共线陷阱行为固定;旋量 SlaterDeterminant 维度错配/占据越界 → `std::invalid_argument`;`forces()` 对 GHF 抛出(本环节无梯度)。
- 领域验证(science.required):塌缩定理 1e-8 Ha;H3 oracle 1e-8 Ha 对 PySCF `scf.GHF`;回归用例 ≤60 s、引用钉定值。
- tests/ 保持单元级秒量级;端到端仅入 regressions/。

## Out of scope

- 复系数 GHF、自旋轨道耦合、GHF 解析梯度/性质(后续链环)。
- Seeger–Pople 稳定性分析本体(链 08);ROHF 实现本体(链 10)。
- UHF→GHF 自动稳定性跟随(本环节初猜靠成员参数破坏共线性)。
- CUDA/torch 后端 GHF。

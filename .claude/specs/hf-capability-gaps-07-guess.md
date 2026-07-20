---
title: SCF 收敛层 I — Hückel/SAP 初猜与 EDIIS/ADIIS 远域混合器
slug: hf-capability-gaps-07-guess
status: approved
chain: hf-capability-gaps
link: 07/21
external-depends: scf-slater-determinant-02-rehome
created: 2026-07-20
---

# SCF 收敛层 I — Hückel/SAP 初猜与 EDIIS/ADIIS 远域混合器

## Summary
为 Gaussian-basis SCF 补齐"远离收敛区"的两条腿:初猜质量与远域混合器。新增扩展 Hückel 初猜与 SAP(原子势叠加)初猜,作为 `InitGuess` 枚举的新取值(选择仍是方法成员,非 Config 包);新增 EDIIS 与 ADIIS 混合器(`internal/diis_mixer.h` 的同级类),并在 `MixStage` 中接入 ADIIS/EDIIS→CDIIS 升级切换。交付后,H2O/def2-SVP 从 Hückel/SAP 出发的迭代数少于 Hcore 初猜,且一个 CDIIS 单独失速的困难体系经 ADIIS→CDIIS 收敛到 PySCF 钉定能量。同时填充 docs/tutorial/ch18 的"初猜与 DIIS 家族"半章(Jensen 3.8 风格叙事)。

## Domain basis
- **扩展 Hückel 初猜**:F_μν = ½K(H_μμ+H_νν)S_μν(μ≠ν),K=1.75;对角元 H_μμ 取 VSIP(价态电离势,静态表,Hartree)。对 F 一次对角化得初始密度。Lehtola, *J. Chem. Theory Comput.* **15**, 1593 (2019), DOI 10.1021/acs.jctc.8b01089。
- **SAP 初猜**:有效单电子势 V_SAP(r) = Σ_A V_A(|r−R_A|),V_A 来自 Lehtola 表列原子径向有效势(exchange-only LDA,静态数据头);F = T + V_SAP 在现有 `MolecularGrid` 上求积装配,一次对角化,无需密度。同上 DOI 10.1021/acs.jctc.8b01089。
- **EDIIS**:min_c Σc_iE_i − ½Σ_ij c_ic_j⟨D_i−D_j|F_i−F_j⟩,s.t. Σc=1, c≥0。Kudin, Scuseria, Cancès, *J. Chem. Phys.* **116**, 8255 (2002), DOI 10.1063/1.1470195。
- **ADIIS**:增广 Roothaan-Hall 能量 min_c 2Σ_i c_i⟨D_i−D_n|F_n⟩ + Σ_ij c_ic_j⟨D_i−D_n|F_j−F_n⟩,同一单纯形约束。Hu & Yang, *J. Chem. Phys.* **132**, 054109 (2010), DOI 10.1063/1.3304922。
- **适用域陷阱(设计约束 + 教程内容)**:EDIIS/ADIIS 只用于远离收敛区;近收敛区用其到底会失速。标准配方是 ADIIS/EDIIS→CDIIS 切换(按 DIIS 误差范数阈值一次性切换)。
- **开壳陷阱**:SAD 型初猜的原子占据平均问题——Hückel/SAP 产生自旋对称密度,开壳层沿用现有 MINAO 路径的 P/2 分裂 + HOMO-LUMO 混合规则,并在头文件与教程中记录该平均化局限。
- 单位:SCF 内部 Hartree;能量写回 Frame 为 eV(现状不变)。

## Design
- **前置**:`scf-slater-determinant-02-rehome` 先落地——`InitGuess` 已迁至 `scf_enums.h`(当前仍在 `scf_result.h:25`;若前置未落地则本环 blocked)。
- `internal::HueckelGuess`(构造函数取 shells/S/元素列表;成员 `RealMatrix F`)与 `internal::SapGuess`(构造函数取 shells/T/几何/网格;成员 `RealMatrix F`)——类+构造函数,禁止 `build_*` 自由函数(镜像 `internal/diis_mixer.h` + `stages/mix_stage.h` 的配对模式)。VSIP 表内联于 `hueckel_guess.h`;SAP 原子势系数在独立静态数据头 `sap_data.h`。不支持的 Z 抛异常,无回退值。
- `InitGuess` 增加 `Hueckel`、`SAP`;分派点是 `mean_field_method.cpp` 初猜块(:352-402),复用现有 `diagonalize` lambda 与 `DensityBuilder`。
- `internal::EDIISMixer` / `internal::ADIISMixer`:与 `DIISMixer` 同级,各自持 (F, D[, E]) 历史(深度 = `diis_max_vectors`);单纯形约束经 c_i = t_i²/Σt_j² 参数化后用小维数值最小化(确定性,维度 ≤ 历史深度);退化历史回退为返回最新 F。
- `MixStage` 升级接线:新枚举 `MixerKind { CDIIS, EDIIS_CDIIS, ADIIS_CDIIS }`(入 `scf_enums.h`);方法新成员 `mixer`(默认 `CDIIS`,行为零变化)与 `mixer_switch_tol`(max|DIIS err| 低于其值即永久切至 CDIIS)。EDIIS 所需 E_i 经 `MixStage` 构造函数新增 `const double&` 电子能引用传入(沿用命名引用风格,无 Context 切片)。
- 单元测试只验证单次调用的矩阵/系数数学;任何跑到收敛的 SCF 属 regressions(tests/ 禁 E2E)。

## Files to create or modify
- `include/atomiverse/cpu/scf/internal/hueckel_guess.h` (new) — `HueckelGuess` + VSIP 静态表
- `include/atomiverse/cpu/scf/internal/sap_guess.h` (new)、`include/atomiverse/cpu/scf/internal/sap_data.h` (new)
- `src/cpu/scf/internal/hueckel_guess.cpp` (new)、`src/cpu/scf/internal/sap_guess.cpp` (new)
- `include/atomiverse/cpu/scf/internal/ediis_adiis.h` (new)、`src/cpu/scf/internal/ediis_adiis.cpp` (new)
- `include/atomiverse/cpu/scf/scf_enums.h`(02-rehome 产物)— `InitGuess` 加值 + `MixerKind`
- `include/atomiverse/cpu/scf/mean_field_method.h`、`src/cpu/scf/mean_field_method.cpp` — 新成员 + 初猜分派 + MixStage 接线
- `include/atomiverse/cpu/scf/stages/mix_stage.h`、`src/cpu/scf/stages/mix_stage.cpp` — 升级切换
- `CMakeLists.txt`(根,HF 源列表 ~:402)— 登记三个新 .cpp
- `tests/cpu/scf/test_init_guess.cpp` (new)、`tests/cpu/scf/test_ediis_adiis.cpp` (new)、`tests/CMakeLists.txt`
- `regressions/scf/h2o_scf_convergence/reg_h2o_scf_convergence.cpp` (new)、`regressions/scf/h2o_scf_convergence/REFERENCE.md` (new)、regressions 构建登记
- `docs/tutorial/ch18-scf-techniques.md` — 填充"初猜与 DIIS 家族"半章(骨架来自链环 01;若缺则先建骨架)

## Tasks
- [ ] Write failing tests for HueckelGuess/SapGuess (tests/cpu/scf/test_init_guess.cpp:F_μν 公式与 VSIP 对角、SAP 矩阵对称性与 H 原子手算参考元、不支持 Z 抛异常;登记 tests/CMakeLists.txt)
- [ ] Write failing tests for EDIISMixer/ADIISMixer (tests/cpu/scf/test_ediis_adiis.cpp:两点历史闭式最优系数、Σc=1 且 c≥0、ADIIS 目标函数值、退化历史无 NaN;登记 tests/CMakeLists.txt)
- [ ] Implement HueckelGuess and SapGuess in include/atomiverse/cpu/scf/internal/{hueckel_guess.h,sap_guess.h,sap_data.h} + src/cpu/scf/internal/{hueckel_guess.cpp,sap_guess.cpp},并在 scf_enums.h 给 InitGuess 加 Hueckel/SAP;登记根 CMakeLists.txt
- [ ] Wire guess dispatch in src/cpu/scf/mean_field_method.cpp 初猜块(开壳沿用 P/2 分裂 + HOMO-LUMO 混合;头文件记录占据平均化陷阱)
- [ ] Implement EDIISMixer/ADIISMixer in include/atomiverse/cpu/scf/internal/ediis_adiis.h + src/cpu/scf/internal/ediis_adiis.cpp(t² 单纯形参数化求解;DOI 注释)
- [ ] Wire MixerKind escalation in stages/mix_stage.{h,cpp} + mean_field_method.{h,cpp}(成员 mixer/mixer_switch_tol,默认 CDIIS 行为零变化;能量经 const double& 引用入 MixStage)
- [ ] Implement regression regressions/scf/h2o_scf_convergence/(H2O/def2-SVP 各初猜迭代数表 + CDIIS 失速/ADIIS→CDIIS 收敛困难例;REFERENCE.md 钉定 PySCF 能量与迭代数;≤60 s;登记 regressions 构建)
- [ ] Write "初猜与 DIIS 家族" half of docs/tutorial/ch18-scf-techniques.md(Jensen 3.8 叙事:发散成因→初猜质量→damping/level-shift 回顾→DIIS→EDIIS/ADIIS 家族与适用域陷阱;molplot 收敛曲线图 + 章末 PySCF 迭代数对照表)
- [ ] Run full check + test suite

## Testing strategy
- **快乐路径(单元,秒级)**:Hückel F 的对角=VSIP、非对角=½K(H_μμ+H_νν)S_μν 逐元素校验;SAP V 矩阵对称、初猜密度 Tr(PS)=N_elec;EDIIS/ADIIS 两点闭式解匹配 ≤1e-10。
- **边界**:单原子/单 AO 历史深度 1;退化历史(重复 F)回退最新 F 无 NaN;不支持元素抛异常;`mixer=CDIIS` 默认下现有 test_scf_solver/test_accelerators 位级不变。
- **领域验证(regressions,非 tests/)**:迭代数 Hueckel<Hcore、SAP<Hcore(H2O/def2-SVP,对齐 PySCF init_guess='huckel'/'atom' 行为量级);困难例 CDIIS 单独在预算内不收敛而 ADIIS→CDIIS 收敛到 REFERENCE.md 钉定的 PySCF 能量(1e-6 Ha)。PySCF 只作外部 oracle 钉数,不作依赖。

## Out of scope
- SAD/原子密度初猜本体(MINAO 已覆盖其角色)、开壳分数占据平均化的修复(仅记录陷阱)。
- DIIS 家族其余成员(LIST、KAIN)、二阶 SCF/信赖域(后续链环)。
- CUDA 后端初猜/混合器移植;ch18 的另半章(二阶方法归 09 环)。

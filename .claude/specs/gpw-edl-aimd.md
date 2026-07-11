---
title: GPW production bridge — self-contained data + ZnO electric-double-layer AIMD
slug: gpw-edl-aimd
status: approved
science_required: true
created: 2026-06-20
revised: 2026-07-11
acceptance: gpw-edl-aimd.acceptance.md
---

# GPW 生产桥接 — 树内数据 + ZnO 电双层 AIMD

## Summary

目标不变：把 GPW 路径推到生产可用，最终在 Atomiverse 内启动 1214 原子 ZnO 电双层
周期性 NVT AIMD（PBE、GTH-PBE 赝势、DZVP-MOLOPT-SR-GTH、dt 0.5 fs、T 330 K、CSVR、
盒子 13.128×10.582×47.054 Å，PBC XYZ）。

**2026-07-11 重构。** 原 spec 无法完成，根因不是"物理没做完"，而是一处结构缺陷：

> **GPW 的每一个真实物理测试都从 `$ATV_CP2K_DATA_DIR/{GTH_POTENTIALS,BASIS_MOLOPT}`
> 读数据**（forces / evaluator / aimd / geoopt / diis / heavy_elements 及全部 parity
> 测试，共 14 个，`tests/CMakeLists.txt` 的 `ATV_CP2K_DATA_DIR` 门控内）。不指向一个
> **外部 CP2K 安装**，GPW 就只剩 12 个纯数学测试（积分、谐函数、collocate、FFT、
> Poisson、XC 格点、玩具 H 求解器）—— **一个真实分子都跑不了，零力校验**。

这与项目绑定约束直接冲突：**代码库不得出现、也不得调用任何第三方 package**，且
`ATV_CP2K_DATA_DIR` 在 CI 中必须保持未设置。因此原 spec 的 4 条 CP2K 判据
（ac-004/006/010/013）不是"暂时做不到"，而是**永远不该那样写**。

重构后的路线：**同样的物理数值，树内的出处**。GTH-PBE 赝势参数与 MOLOPT 基组都是
**已发表的表格**（CP2K 的数据文件本身就是这些论文的转录）；照 `cpu/scf/basis_data/`
的既有模式（`pople_631g.h` 等 6 张 constexpr 表，数值一次跑出后钉死 + 出处注释）把
它们做成树内数据。已跑出的 CP2K 能量/受力保留为 **constexpr 锚点常量**（钉死数值，
不调用 CP2K）—— 这正是"跑完然后 hardcode"规则，与 `pople_631g.h` 钉 PySCF 数值同构。
`dftd3_params.h` 已是树内参数表，先例俱在。

**"对齐 CP2K"的绑定依然成立 —— 对齐的是数值，不是依赖。**

## Domain basis

科学引用（`science_required: true`）。公开单位 eV/Å/fs；内部 Hartree/Bohr。

**GTH 赝势数据（M0 的核心）。** 局域部分
`V_loc(r) = -(Z_ion/r)·erf(r̄/√2) + exp(-r̄²/2)·Σ_k C_k r̄^{2(k-1)}`，`r̄ = r/r_loc`；
投影子 `p_i^l(r) = N_i^l · r^{l+2(i-1)} · exp(-r²/2r_l²)`（HGH 归一化）；
非局域 `V_nl = Σ_l Σ_ij h^l_ij |p_i^l⟩⟨p_j^l|`。PBE 用 Krack 优化参数集（非 HGH-LDA）。
EDL 五元素所需通道：H-q1、O-q6、Zn-q12（l=0,1,2 含 d）、Na-q9（半芯）、Cl-q7。
参考：Goedecker-Teter-Hutter *PRB* **54**, 1703 (1996)；Hartwigsen-Goedecker-Hutter
*PRB* **58**, 3641 (1998)；**Krack, *TCA* **114**, 145 (2005)**（PBE 参数表出处）。

**MOLOPT 基组数据（M0）。** DZVP-MOLOPT-SR-GTH 收缩系数与指数。
参考：**VandeVondele & Hutter, *JCP* **127**, 114105 (2007)**, DOI 10.1063/1.2770708。

**解析 GPW 力。** `F_A = F_HF(loc) + F_Pulay + F_nl + F_grid(H/XC) + F_Ewald`，
`F_Pulay = -Tr(P·dH_core/dR_A) + Tr(W·dS/dR_A)`，`W = Σ_i f_i ε_i c_i c_i†`。
参考：VandeVondele 等, *CPC* **167**, 103 (2005)。
**唯一的正确性判据是中心有限差分 oracle（步长 1e-4 Bohr）逐分量比对 —— 无参考、
自足、不需要任何外部引擎。** 这一直是项目内真正的力校验；CP2K 力参考从未增加过
信息量，只增加了一个被禁的依赖。

**CSVR / Bussi 恒温器。** 见 Bussi-Donadio-Parrinello *JCP* **126**, 014101 (2007)。
（已实现并验证，M4 完成。）

**DFT-D3(BJ)。** `E_disp = -Σ_{A<B} Σ_{n=6,8} s_n·C_n^AB / (R^n + (a1·R0+a2)^n)`；
解析力必须含 CN 导数第二遍。PBE-D3(BJ) 参数 s6=1.0, a1=0.4289, s8=0.7875, a2=4.4407。
参考：Grimme 等 *JCP* **132**, 154104 (2010)；Grimme-Ehrlich-Goerigk *JCC* **32**, 1456 (2011)。
参数与 C6 表已在树内（`gpw/internal/dftd3_params.h`）。

## Design

**M0 — 树内 GTH + MOLOPT 数据（阻塞性前置，本次重构新增）。**
新增 `include/atomiverse/cpu/gpw/basis_data/`，镜像 `cpu/scf/basis_data/` 的形态：
- `gth_pbe.h` — 五元素 GTH-PBE 赝势 constexpr 表（Z_ion、r_loc、C_k、每 l 的 r_l 与
  h^l_ij），出处注释指向 Krack 2005。
- `molopt_dzvp.h` — DZVP-MOLOPT-SR-GTH 收缩 constexpr 表，出处注释指向
  VandeVondele-Hutter 2007。
- `gpw_data_table.h` — 按元素+名字查表的入口，**不支持的元素/通道抛异常（无回退）**。

`cp2k_io.h` / `cp2k_basis.h` 的**解析器保留**（它们元素无关且已验证），但测试与生产
路径改走树内表；解析器降级为"可选摄入外部文件"的工具，不再是任何测试的前置。
删除 `tests/CMakeLists.txt` 的 `ATV_CP2K_DATA_DIR` 门控与 CMake 选项；14 个被门控的
测试按下述规则各归其位。

**M0b — 测试归位（服从仓库铁律）。** `tests/` 只放单元测试（单函数/单方法，ms 级）；
一切 e2e / 引擎间平价搬到 `benchmarks/qc/gpw/validation/`（acceptance 时跑，ctest 永不跑）。
- 留在 `tests/`（单元）：`cp2k_io`、`basis`、`heavy_elements`（→ 改为对树内表断言）、
  `real_harmonics`、`gth_integrals` 等纯数学项，以及 `forces` 的**小体系** FD 断言。
- 迁往 `benchmarks/`：`h2_cp2k`、`ar_cp2k`、`h2o_tutorial`、`zno_cp2k`、`nacl_cp2k`、
  `edl_frag`、`geoopt`、`aimd`、`multigrid_scf`（跑到收敛的完整 SCF/AIMD = e2e）。
  改名去掉 `_cp2k` 后缀（数值是钉死的锚点，不是一次外部调用）。

**M1–M4 — 已完成**（ExtXYZ molrs 桥、重元素摄入、任意-l 实谐函数与 d 投影子、
DIIS、解析力、FixCSVR）。M0 完成后这些的验证不再依赖外部数据目录。

**M5 — 1214 原子规模。** 已实测定位真正的墙：`gth_matrix_element` 对每个 (μ,ν) 遍历
所有原子无筛选 → O(N³)，N=24 时占 81%。已做：S/T/Vpp + ionic_density OpenMP → 11.8×。
待做：**原子循环筛选**（GTH PP 局域 → 跳过远原子，~50–100×）+ 对循环重叠筛选（~3×）。
（此前"稠密对角化是瓶颈"的判断已被实测推翻；稠密 OT 是性能倒退但正确，门控保留；
稀疏 LA 因重叠密度 24–34% 被否决 —— 均为已记录的实测结论，不重做。）

**M6 — DFT-D3(BJ) 接线。** 参数表已在树内；补解析力的 CN 导数第二遍 + FD 力自检。

约束（贯穿）：TDD；公开符号 Doxygen（`@brief`/`@param` 带单位/`@return`）；
无工厂函数；无回退值；**树内不出现第三方 package，构建与测试期不读任何第三方文件**。

## Files to create or modify

M0：
- `include/atomiverse/cpu/gpw/basis_data/gth_pbe.h` (new)
- `include/atomiverse/cpu/gpw/basis_data/molopt_dzvp.h` (new)
- `include/atomiverse/cpu/gpw/basis_data/gpw_data_table.h` (new)
- `tests/cpu/gpw/test_cpu_gpw_data_table.cpp` (new) — 树内表 vs 文献值 + 未支持元素抛异常
- `tests/CMakeLists.txt` — 删 `ATV_CP2K_DATA_DIR` 门控
- `CMakeLists.txt` — 删 `ATV_CP2K_DATA_DIR` 选项
- `CLAUDE.md` — 删该选项行

M0b：
- `benchmarks/qc/gpw/validation/` (new) — 迁入 9 个 e2e/锚点测试
- `benchmarks/CMakeLists.txt` — 注册 `val_gpw_*`
- `tests/cpu/gpw/test_cpu_gpw_{forces,heavy_elements}.cpp` — 改走树内表

M5：
- `src/cpu/gpw/gpw_solver.cpp` — `gth_matrix_element` 原子循环筛选

M6：
- `src/cpu/gpw/internal/dftd3.cpp` — CN 导数力项

## Tasks

- [x] M1 ExtXYZ molrs 桥 + `cpu::load_extxyz`
- [x] M2 重元素摄入（Zn/Na/Cl GTH+MOLOPT）；任意-l 实谐函数 + l=2 d 投影子
- [x] M2b GPW DIIS/Pulay（H2 49→5、ZnO 157→16 迭代）
- [x] M3 解析 GPW 力（W 矩阵 + 五项分解），GPWEvaluator FD→解析
- [x] M4 `cpu::md::FixCSVR`（Bussi α² 传播子 + H̃ 记账）
- [x] M0 树内 GTH-PBE + DZVP-MOLOPT-SR-GTH + D3 数据表（`cpu/gpw/basis_data/`；
      逐 token 对齐已发表参数表；未支持元素抛异常）
- [x] M0b 删 `ATV_CP2K_DATA_DIR` 选项与门控；14 个测试按单元/e2e 归位
      （4 个单元留 `tests/`，14 个 e2e 迁 `benchmarks/qc/gpw/validation/`）；
      examples 09/10 也改走树内表。全树无第三方数据读取。
- [x] **M3b 解析力 FD oracle 重新绿 —— 并揪出一个 17 个 commit 前就存在的真 bug**：
      性能攻坚把局域赝势挪到格点，力却仍取逐对解析导数 → 在求一个和 SCF 收敛的
      泛函不同的泛函的梯度（ZnO 0.534 Ha/Bohr）。修：共享内核
      `internal/local_pp.{h,cpp}`（Veff 带 V_loc^SR + 补 HF 项 + 去掉重复解析项），
      隐藏环境开关 `ATV_GPW_ANALYTIC_LOC` 升为 `GpwConfig.local_pp_on_grid`。
      修后 H2 <1e-5、周期 H2O 6.8e-5、ZnO 1.7e-5（残差是 FD oracle 自身的 O(h²)
      截断，实测 4×/倍步长）。
- [x] Doxygen（extxyz / gpw_forces / FixCSVR / 新数据表 / local_pp）
- [ ] **M5 `gth_matrix_element` 原子筛选 → 1214 原子单点可行**（阻塞收官）
- [ ] **ZnO 对角化 SCF 不收敛**（ac-014）：400 步撞上限，且比 OT 的能量高 ~4 mHa；
      mix=0.1 无效。OT 本身正常（H2 对角化一致到 2.6e-14）。这是新暴露的既有缺陷 ——
      该测试此前挂在 CP2K 门控后，从未在 CI 跑过。
- [ ] M5b EDL 收官：`final_frame.xyz` + GPW-CSVR NVT 短轨迹（能量守恒 + 无漂移诊断）
- [ ] M5c 全系统 wall-time / 峰值内存画像（benchmarks）
- [ ] M6 D3(BJ) CN 导数力项 + FD 自检

## Testing strategy

- **单元（`tests/`，ms 级）**：树内数据表 vs 文献值、未支持元素抛异常、实谐函数正交性、
  GTH 积分、collocate/FFT/Poisson/XC 格点、小体系解析力 vs FD 逐分量。
- **e2e / 锚点（`benchmarks/qc/gpw/validation/`，仅 acceptance 跑）**：H2/Ar/H2O/ZnO/NaCl
  单点能对**钉死的 constexpr 锚点**（一次性跑出、带出处注释，构建与测试期不调用任何外部
  引擎）；geoopt / AIMD 短轨迹；EDL 收官。
- **无参考的物理判据（最高价值）**：解析力 vs 6N+1 中心 FD oracle；AIMD 能量守恒 /
  CSVR 的 H̃ 无漂移。这些不需要任何外部引擎，失败即为本仓库缺陷。
- 完整 `cmake --build` + `ctest`（`ATV_BUILD_GPW`，含 RHF+PW+FFTW3），**且在
  `ATV_CP2K_DATA_DIR` 不存在的情况下全绿**（这是本次重构的验收信号）。

## Out of scope

- 完整 2000 步生产运行（下游使用，非代码交付物）。
- GAPW、GPW 中的杂化/HFX、GPW k 点（仅 Γ）、CUDA GPW 后端、金属/展宽。
- 通用 CP2K `.inp` 解析器。
- 稀疏 LA（重叠密度 24–34%，实测不划算 —— 已否决，不重开）。
- D3 之外的色散变体。
- 手写 C++ XYZ 解析器（所有 I/O 走 molrs）。

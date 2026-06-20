---
title: GPW production bridge — CP2K ZnO electric-double-layer AIMD
slug: gpw-edl-aimd
status: approved
science_required: true
created: 2026-06-20
acceptance: gpw-edl-aimd.acceptance.md
builds_on: gpw-gaussian-plane-wave
---

# GPW 生产桥接 — CP2K ZnO 电双层 AIMD

## Summary

本 spec 把已合并的 GPW 路径（[[gpw-gaussian-plane-wave]]，Phase 0–8 完成：GPW SCF、native PBE、MT Poisson、有限差分力、AIMD 接线、`examples/09_gpw`、`tests/cpu/gpw/`，但仅在 H/O/Ar ≤3 原子的微系统上验证）推到**生产可用**，目标是把 CP2K 的 `ZnO_MD.inp` 电双层 AIMD（1214 原子周期性 NVT：ZnO slab 192 Zn + 192 O、270 H₂O、10 NaCl 对；PBE；GTH-PBE 赝势；DZVP-MOLOPT-SR-GTH 基组；dt 0.5 fs；T 330 K；CSVR 恒温器；盒子 13.128×10.582×47.054 Å，PBC XYZ）在 Atomiverse 内启动。本 spec 兑现兄弟 spec 中延期的 Phase 8「解析力」和 Phase 9「benchmark 迁移」两个 TODO，并补齐剩余的生产缺口。最终用户可见结果：把 `final_frame.xyz` 读入一个 `cpu::Frame`，配置 GPW + CSVR NVT，跑出一条短轨迹，逐步能量与受力在容差内与 CP2K 参考一致，并产出全系统的 wall-time/内存画像。本 spec **不**复制兄弟 spec 已完成的工作。

整体按 M1–M6 里程碑组织，依赖顺序排列，镜像兄弟 spec 的 M/Phase 房屋风格。M6（DFT-D3）显式置于关键路径之外。

## Domain basis

科学引用（`science_required: true`）。公开单位 eV/Å/fs；内部 Hartree/Bohr。

**解析 GPW 力（M3 的核心）。** 总力按 CP2K Quickstep GPW 分解：
`F_A = F_HF(loc) + F_Pulay + F_nl + F_grid(H/XC) + F_Ewald`。
其中 Pulay 项 `F_Pulay = -Tr(P · dH_core/dR_A) + Tr(W · dS/dR_A)`，能量加权密度矩阵 `W = Σ_i f_i ε_i c_i c_i†`（由收敛后的 `C, ε, f` 构造）。`dT`、`dS`、`d⟨φ|V_loc|φ⟩` 经 Gaussian 递推解析求出。Hellmann–Feynman 项来自移动的 Gaussian 核电荷 `n_c^A` 与格点上的局域短程 `V_loc`。非局域 GTH 投影子导数需同时取 `⟨∂φ|p⟩` 与 `⟨φ|∂p⟩`，按牛顿第三定律 ±F 对称累加。
参考：VandeVondele, Krack, Mohamed, Parrinello, Chassaing, Hutter, *Comput. Phys. Commun.* **167**, 103 (2005), DOI 10.1016/j.cpc.2004.12.014。验证基准：与中心有限差分（步长 1e-4 Bohr）逐分量比对。

**GTH 赝势（M2 重元素路径）。** 局域部分
`V_loc(r) = -(Z_ion/r)·erf(r̄/√2) + exp(-r̄²/2)·Σ_k C_k r̄^{2(k-1)}`，`r̄ = r/r_loc`；
投影子 `p_i^l(r) = N_i^l · r^{l+2(i-1)} · exp(-r²/2r_l²)`（HGH 归一化）；
非局域 `V_nl = Σ_l Σ_ij h^l_ij |p_i^l⟩⟨p_j^l|`。Zn-q12 需要 `l = 0,1,2`（含 d 通道）。投影子数目必须从实际 CP2K `GTH_POTENTIALS` 条目读出，**不得硬编码**。XC=PBE 时使用 Krack PBE 优化参数集（非 HGH-LDA），逐字使用文件中的 `r_loc / C_i / {r_l, h^l_ij}` 以求 bit-parity。
参考：Goedecker-Teter-Hutter, *PRB* **54**, 1703 (1996)；Hartwigsen-Goedecker-Hutter, *PRB* **58**, 3641 (1998)；Krack, *TCA* **114**, 145 (2005)。

**CSVR / Bussi 恒温器（M4）。** 随机速度重标 SDE
`dK = (K̄-K)dt/τ + 2√(K·K̄/N_f)·dW/√τ`，离散传播子
`α² = e^{-Δt/τ} + (K̄/(N_f K))(1-e^{-Δt/τ})(R₁² + Σ_{i≥2}R_i²) + 2 e^{-Δt/2τ} √((K̄/(N_f K))(1-e^{-Δt/τ})) R₁`，
其中 `Σ_{i≥2}R_i²` 抽自 `Gamma((N_f-1)/2, 2)`；`v ← α v`；并累积守恒量记账项（有效能量 H̃，Eq.15，作无漂移诊断）。`τ = TIMECON`。`N_f` 记账须扣除 COM/约束自由度。
参考：Bussi, Donadio, Parrinello, *JCP* **126**, 014101 (2007), DOI 10.1063/1.2408420。

**DFT-D3(BJ)（M6，关键路径外）。**
`E_disp = -Σ_{A<B} Σ_{n=6,8} s_n · C_n^AB / (R_AB^n + (a1·R0^AB + a2)^n)`，`R0^AB = √(C8/C6)`；
配位数依赖的 `C6` 插值（k1=16, k2=4/3, k3=4）。解析力**必须**包含非局域 CN 导数的第二遍 `(∂E/∂C6)·(∂C6/∂CN)·(∂CN/∂R)`，否则受力静默错误。CP2K PBE-D3(BJ) 参数 `s6=1.0, a1=0.4289, s8=0.7875, a2=4.4407`（锁定前对照 CP2K `dftd3.dat` 核验）。
参考：Grimme 等, *JCP* **132**, 154104 (2010)；Grimme-Ehrlich-Goerigk, *JCC* **32**, 1456 (2011)。

## Design

里程碑（依赖顺序），每个里程碑是可独立测试/合并的实现单元：

**M1 — 经 molrs 的 ExtXYZ 加载（不手写 C++ 解析器）。** 架构约束（用户绑定）：**所有 I/O 走 molrs**。molrs 核心已解析 ExtXYZ（`molrs::io::data::xyz::read_xyz_frame` / `read_xyz_traj` / `XYZReader`，`parse_simbox` 从 `Lattice="..."` 头构造 SimBox）。缺口仅在 CXX 桥未暴露读取：`molrs-cxxapi/src/lib.rs` 当前桥接写（`xyz_write`/`xyz_append`）与 Zarr 读（`molrec_read_zarr_first_frame -> Box<FrameRef>`）。交付：(a) 在 molrs-cxxapi 新增 `xyz_read_first_frame(path) -> Box<FrameRef>`（可选轨迹变体）包裹 `read_xyz_frame`；(b) atomiverse `cpu::io` 薄包装 `load_extxyz(path) -> Frame`，镜像 `cpu/io/lammps_data.h` 风格，复用既有 `cpu::materialize(const molrs::FrameRef&) -> cpu::Frame` 路径（`include/atomiverse/cpu/frame/materialize.h`），把 `Lattice` → `frame.box`（`Mat3`）、每原子 Z → `type`。molrs 是位于 `/nobackup/proj/disk/teoroo/personal/jicli594/work/molcrafts/molrs` 的独立仓库（corrosion 构建，CXX 头是 `include/atomiverse/molrs.hpp`）；桥接编辑落在 molrs，atomiverse 重新对其编译。

**M2 — Zn/Na/Cl 重元素验证（GTH + MOLOPT）。** GTH/MOLOPT 解析器（`cpu/gpw/internal/cp2k_io`、`cp2k_basis`）元素无关但仅 H/O/Ar 验证过。摄入并验证 GTH-PBE-q12（Zn，需 l=0,1,2，d 通道；从实际 CP2K `GTH_POTENTIALS` 条目读取投影子数）、q9（Na，l=0,1 半芯）、q7（Cl，l=0,1），以及五元素的 DZVP-MOLOPT-SR-GTH 基组。门槛：用逐元素阶梯（ZnO → NaCl 离子对 → H₂O，每级单独对 CP2K 单点能）验证。

> **DISCOVERED (ac-004 rung-1 ZnO 验证，2026-06-20)：GPW core 角动量上限缺口（两处既有 TODO）。** 解析器无需改（ac-003 已 GREEN），但 ZnO 单点 **段错误**。根因:
> 1. **基组 f 壳层缺失**:`gpw_solver.cpp` 的 `harmonics()` 实-球谐→笛卡尔表只到 d（l=2）。Zn DZVP-MOLOPT-SR-GTH 含 f（l=3），`harmonics()[3]` 越界。→ 需补 l=3 的 7 个实-球谐笛卡尔展开项。
> 2. **d 非局域投影子缺失**:`gth_integrals.cpp` 的 `gth_projector_overlap`（`return 0.0; // TODO: l>=2`）与 `gth_nonlocal`（`n_m=(l==0)?1:(l==1?3:0)`）只支持 s/p 投影子。Zn-q12 有 l=2（d）投影子，当前被静默丢为 0 → 物理错误。→ 需补 l=2 投影子（5 个 m 分量的实 d 球谐 ⟨φ|p⟩，仿 l=1 路径）。
>
> 范围：五元素中**只有 Zn** 触发这两处（O/Cl 是 d 基组走 `harmonics()[2]`，已支持;Na/Cl 投影子仅 s/p）。这是运行整个 ZnO 目标体系(192 Zn)的硬前置。
>
> **下游已确认通用**(2026-06-20)：`overlap`/`kinetic`(Hermite 递推 `hermite_E`/`overlap_1d(i,j)`)、`collocate`(显式单项式循环)、`moment3d`(通用 `extra_dim`+`k`)均对任意笛卡尔角动量幂通用,**无需改动**。因此任意-l 支持只需泛化上述两处角动量生成器。
>
> **决策(用户)：实现任意-l 通用**(非仅 f/d)。Use A=未归一化实-固体谐(基组,span 决定,常数抵消);Use B=归一化实球谐(投影子,Σ_m 不变性 ⇒ 任一正交约定皆可)。两者均锚定复现现有 l≤2 表项作单元校验。

**M3 — 解析 GPW 力（阻塞项）。** 现有 `gpw_forces_fd`（`cpu/gpw/gpw_solver.h`）为有限差分：每步 6N+1 次 SCF，在 1214 原子下不可行。实现解析力 `cpu/gpw/internal/gpw_forces.{h,cpp}`（+ src），一次 SCF + 一次梯度遍历。按 Domain basis 的五项分解；W 矩阵在 `gpw_solver` 中由收敛 `C,ε,f` 构造；`dT/dS/d⟨φ|V_loc|φ⟩` 复用 `cpu/scf` 单电子/Boys/Rys 原语与 `gth_integrals`。验证：对既有 6N+1 FD oracle ≤ 1e-5 eV/Å 且对 CP2K 力参考 ≤ 1e-4 Ha/bohr 逐分量。经 `GPWEvaluator`（已有 `Capability::Force`）从 FD 切换到解析。

**M4 — CSVR / Bussi 恒温器。** 新增 `cpu::md::FixCSVR`（镜像 `cpu/md/fix_langevin.h`、`nose_hoover.h`；Thermostat 角色）。按 Domain basis 的 α² 传播子；`v ← αv`；累积守恒量记账（H̃，无漂移诊断）；`τ = TIMECON`；`N_f` 记账扣除 COM/约束自由度。先在便宜的经典 MD 系统上验证正则动能分布，再接入 GPW-AIMD。

**M5 — 周期性大系统规模/性能 + 收官运行。** GPW 最大测过 N=3；1214 原子周期性 CPU 的正确性、性能、内存均未测。**前置(M2 发现)**：GPW 求解器只有线性混合,Zn 体系已需 mix=0.1/~157 迭代;1214 原子靠线性混合不会收敛 —— 收官前必须先给 GPWSolver 加 DIIS/Pulay 混合(`cpu::scf` 已有 DIIS 可参考)。用 OpenMP（在既有 GPW/PW MPI 模式允许处加 MPI）剖析并并行 GPW 热点（collocation/integration、FFT、KS 构建、新解析力遍历）。收官验收：读 `final_frame.xyz`，配置 GPW（PBE、五个 GTH-PBE 赝势、DZVP-MOLOPT-SR-GTH、合适 CUTOFF）、CSVR NVT、dt 0.5 fs，跑几步短轨迹，逐步能量+受力在容差内匹配 CP2K 参考，并产出全系统 wall-time/内存画像。

**M6 — DFT-D3 色散（PBE-D3，后续/低优先）。** 加 D3(BJ)（CP2K PBE 默认）。作为力场式修正，置于可复用组件（`potentials/` 或小型 `cpu/dispersion` 模块），供 MD evaluators 复用。明确关键路径外：AIMD 可不带它启动（物理略有差异）。

约束（贯穿）：TDD（RED 测试引用 ac-NNN 在 GREEN 前）；公开符号 Doxygen（`@brief/@param`（带单位）/`@return`）；公开单位 eV/Å/fs，内部 Hartree/Bohr；无工厂函数；无回退值（查表对不支持输入抛异常）。

## Files to create or modify

M1：
- `/nobackup/proj/disk/teoroo/personal/jicli594/work/molcrafts/molrs/molrs-cxxapi/src/lib.rs` — 新增 `xyz_read_first_frame` CXX 桥接 fn（molrs 仓库）
- `include/atomiverse/cpu/io/extxyz.h` (new) — `cpu::load_extxyz` 声明
- `src/cpu/io/extxyz.cpp` (new) — 实现，复用 `cpu::materialize`
- `tests/cpu/gpw/test_cpu_io_extxyz.cpp` (new)

M2：
- `tests/cpu/gpw/test_cpu_gpw_heavy_elements.cpp` (new) — Zn/Na/Cl GTH+MOLOPT 摄入 + ZnO+水团簇单点
- `include/atomiverse/cpu/gpw/internal/cp2k_io.h` / `cp2k_basis.h` — 仅在 d 通道/半芯投影子计数需要修正时改

M3：
- `include/atomiverse/cpu/gpw/internal/gpw_forces.h` (new)
- `src/cpu/gpw/internal/gpw_forces.cpp` (new)
- `include/atomiverse/cpu/gpw/gpw_solver.h` — 暴露 W 矩阵构造 / 解析力入口
- `include/atomiverse/cpu/gpw/gpw_evaluator.h` — FD→解析切换开关
- `tests/cpu/gpw/test_cpu_gpw_forces.cpp` — 扩展：解析 vs FD oracle + CP2K 参考

M4：
- `include/atomiverse/cpu/md/fix_csvr.h` (new)
- `src/cpu/md/fix_csvr.cpp` (new)
- `tests/cpu/md/test_cpu_md_csvr.cpp` (new)

M5：
- `include/atomiverse/cpu/gpw/internal/collocate.h`、`grid_fft.h` — OpenMP/MPI 热点并行
- `benchmarks/gpw/znoedl/validation/` (new) — 收官参考夹具 + 运行（镜像 `benchmarks/pw/si/validation`）
- `examples/10_gpw_znoedl.cpp` (new) — 加载 + 启动短轨迹示例
- `tests/cpu/gpw/test_cpu_gpw_znoedl.cpp` (new) — 短轨迹 CP2K 逐步 parity

M6（关键路径外）：
- `include/atomiverse/cpu/dispersion/dftd3.h` (new)
- `src/cpu/dispersion/dftd3.cpp` (new)
- `tests/cpu/test_cpu_dftd3.cpp` (new)

## Tasks

- [x] Write failing tests for ExtXYZ molrs bridge + Frame load round-trip (tests/cpu/io/test_cpu_io_extxyz.cpp)
- [x] Implement xyz_read_first_frame CXX bridge in molrs-cxxapi and cpu::load_extxyz wrapper reusing cpu::materialize (include/atomiverse/cpu/io/extxyz.h, src/cpu/io/extxyz.cpp) — also fixed truncated symbol_for_z (1-18 → 1-118)
- [x] Write failing tests for Zn/Na/Cl GTH+MOLOPT ingest and per-element-ladder GPW single point vs CP2K — ingest (ac-003) GREEN; ladder (ac-004) GREEN: ZnO 0.10 mHa, NaCl 0.03 mHa, H2O 0.43 mHa (existing tutorial) vs CP2K. All 5 EDL elements {Zn,O,H,Na,Cl} validated. (test_cpu_gpw_zno_cp2k, test_cpu_gpw_nacl_cp2k)
- [x] Implement heavy-element projector-count and d-channel ingest in cp2k_io / cp2k_basis — NO change needed; parsers already element-generic (Zn d/f + SR-MOLOPT verified by ac-003)
- [x] DISCOVERED: arbitrary-l real harmonics — new cpu/gpw/internal/real_harmonics.{h,cpp} (`real_solid_harmonics`=Use A, `real_spherical_harmonics_normalized`=Use B via exact sphere monomial integral); wired into harmonics() + gth projector; l≤4 orthonormality + l≤3 forms unit-tested (test_cpu_gpw_real_harmonics); fixes Zn f-shell segfault + d-projector
- [x] DISCOVERED: l=2 (d) GTH nonlocal projector — folded into the above (gth_projector_overlap/gth_nonlocal now n_m=2l+1 any l; TODO lines gone)
- [x] Downstream kernels confirmed general (Hermite/monomial); ZnO rung GREEN vs CP2K (|dE|=0.10 mHa, mix=0.1, ~157 iters); NaCl + H2O rungs still pending
- [x] DISCOVERED (ZnO convergence): added GPW Pulay/DIIS — new cpu/gpw/internal/diis.{h,cpp} (single-channel, GPW-local, no scf::internal include; X=S^{-1/2} via cpu::linalg, computed once; KS-matrix extrapolation + commutator-residual convergence). GpwConfig.diis_depth=6 default (0=legacy linear). H2 49→5, ZnO 157→16 iters; energies bit-unchanged; suite 60s→8.7s. Removes the M5 linear-mixing blocker (still to prove at 1214-atom scale).
- [ ] Write failing tests for analytic GPW forces vs 6N+1 FD oracle and CP2K reference (tests/cpu/gpw/test_cpu_gpw_forces.cpp)
- [ ] Implement analytic GPW forces with W-matrix, HF/Pulay/nonlocal/grid/Ewald terms (include/atomiverse/cpu/gpw/internal/gpw_forces.h, src/cpu/gpw/internal/gpw_forces.cpp)
- [ ] Wire FD-to-analytic force switch into GPWEvaluator (include/atomiverse/cpu/gpw/gpw_evaluator.h)
- [ ] Write failing tests for FixCSVR canonical kinetic-energy distribution (tests/cpu/md/test_cpu_md_csvr.cpp)
- [ ] Implement cpu::md::FixCSVR thermostat with H-tilde bookkeeping (include/atomiverse/cpu/md/fix_csvr.h, src/cpu/md/fix_csvr.cpp)
- [ ] Add Doxygen with units to all new public symbols (extxyz, gpw_forces, FixCSVR)
- [ ] Write failing capstone test loading final_frame.xyz, running short GPW-CSVR NVT trajectory vs CP2K per-step parity (tests/cpu/gpw/test_cpu_gpw_znoedl.cpp)
- [ ] Profile and parallelize GPW hot paths (collocate/grid_fft/KS/force) with OpenMP/MPI; add benchmark + example (benchmarks/gpw/znoedl/validation, examples/10_gpw_znoedl.cpp)
- [ ] Implement DFT-D3(BJ) dispersion component with CN-derivative force term (include/atomiverse/cpu/dispersion/dftd3.h, src/cpu/dispersion/dftd3.cpp) — out of critical path
- [ ] Run full check + test suite

## Testing strategy

- **Happy path（M1）**：`load_extxyz(final_frame.xyz)` 返回 1214 原子 Frame，盒子 = 13.128×10.582×47.054 Å，五种元素的 Z/type 正确，往返 atom count 一致。
- **Happy path（M4）**：FixCSVR 接入经典 NVT，平均温度收敛到目标 T。
- **Edge cases**：M1 缺 `Lattice` 头 → 命名异常（无回退）；M2 未支持元素/投影子计数缺失 → 查表抛异常；M4 `N_f` 在有 COM 移除时正确扣除；M3 单原子（无 Pulay 交叉项）退化情形受力为零。
- **Domain validation（scientific）**：
  - M2：ZnO+水团簇 GPW 单点能 ≤ ~1e-4 Ha（匹配 grid+basis+PP+XC）vs CP2K。
  - M3：解析力对内部 6N+1 中心 FD ≤ 1e-5 eV/Å 逐分量，且对 CP2K 力参考 ≤ ~1e-4 Ha/bohr 逐分量。
  - M4：CSVR 在便宜系统上的动能分布通过正则分布卡方/方差检验。
  - M5：`final_frame.xyz` 短轨迹逐步能量与受力在容差内匹配 CP2K 参考。
  - M6：D3(BJ) 能量与解析力（含 CN 导数）对 CP2K D3 数值 FD 一致（关键路径外）。
- **Performance（M5）**：全系统单步 wall-time 与峰值内存记录入 benchmark 报告。
- 完整 `cmake --build` + `ctest`（ATV_BUILD_GPW，含 ATV_BUILD_RHF+ATV_BUILD_PW+FFTW3）。

## Out of scope

- 完整 2000 步生产运行（下游使用，非代码交付物；门槛是「启动 + 短轨迹 parity + 全系统画像」）。
- GAPW、GPW 中的杂化/HFX、GPW k 点（仅 Γ）、CUDA GPW 后端、金属/展宽。
- 通用 CP2K `.inp` 解析器（配置经 `GPWConfig` + 摄入的数据文件）。
- 波函数重启/外推调优（超出 gpw-gaussian-plane-wave 已有的 ASPC）。
- D3 之外的色散变体（仅 D3(BJ)，且 M6 置于关键路径外）。
- 手写 C++ XYZ 解析器（绑定约束：所有 I/O 走 molrs）。

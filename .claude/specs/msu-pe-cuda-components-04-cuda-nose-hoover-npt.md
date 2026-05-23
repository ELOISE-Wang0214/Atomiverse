---
title: msu-pe-cuda-components-04-cuda-nose-hoover-npt — CUDA NoseHooverNPT 集成 NPT 积分器
status: approved
created: 2026-05-23
parent: msu-pe-validation.md
chain: msu-pe-cuda-components-* (link 4/8)
---

# msu-pe-cuda-components-04-cuda-nose-hoover-npt — CUDA NoseHooverNPT 集成 NPT 积分器

## Summary

本规范是 `msu-pe-cuda-components-*` 链的第 4/8 环（链中最复杂的一环），承担 `msu-pe-validation.md` Phase G.7 的 CUDA 端 NPT 缺口。当前 CUDA 后端仅有分立的 `cuda::md::NoseHoover`（NVT 链式恒温器）+ `cuda::md::Barostat`（Berendsen 各向同性压浴），无法跑通 MSU 非晶态 PE 拉伸协议的 stages 2-5（各向同性 NPT 退火）和 stage 5（y/z 各向异性 NPT + x 自由轴）。本环把 CPU 端已与 LAMMPS `fix npt` 对齐（commit `1e2260f`：`mtk_term1 = Σ_{α active} 2·KE_α/(pdim·N)` 与 `mtk_term2 = Σ_{α active} ω̇/(pdim·N)`，来源 `fix_nh.cpp:2541-2569`）的 `cpu::md::NoseHooverNPT` 完整移植到 CUDA：单一组件同时承担 Integrator + Thermostat + Barostat 三重职责，支持 `couple_axis[3]` 逐轴启停的 M 链 NHC（mtchain/mpchain）+ MTK 桶活塞、`drag` 经验阻尼、以及由 `set_run_span(nfirst, nlast)` 驱动的 `T_start → T_stop` 线性温度斜坡（stage 3 的 500→100 K）。CUDA 实现保留 AV 特有的 `W = 3·omega_mass`、`G_eps = 3·f_omega` 约定使全部三轴耦合时位级（algebraically）回退到原各向同性标量整合器；逐轴动能 `KE_α` 与压力对角分量 `σ_αα` 通过专门的归约核读出，确保 Trotter 分裂对称、`mtk_term1`/`mtk_term2` 与 LAMMPS 一致。完成后链环 06（CUDA 路径基准）即可在 stages 2-5 调用 `cuda::md::MD::set_integrator(NoseHooverNPT{…})`。

## Domain basis

算法引用 Martyna-Tobias-Klein NPT 集成器（Tuckerman 2006, *J. Phys. A* 39, 5629, DOI `10.1088/0305-4470/39/19/S18`；Shinoda 2004, *Phys. Rev. B* 69, 134103, DOI `10.1103/PhysRevB.69.134103`）和 LAMMPS `fix_nh.cpp` 参考实现（GitHub `lammps/lammps` release `22 Jul 2025`，本机已克隆于 `/cfs/klemming/scratch/j/jichenli/salamander_2025/lammps-src/src/fix_nh.cpp`）。规范权威定义见 `include/atomiverse/cpu/md/nose_hoover_npt.h` 中的算法文档块；CUDA 实现必须按位语义复刻其 Trotter 分裂 P1-P11（外层 NHC 在前半步 baro-then-atoms、后半步 atoms-then-baro 以保持时间反演对称）和以下数值公式：

- 桶活塞驱动力（每个启用轴 α，对应 `fix_nh.cpp::nh_omega_dot` 中 `pdim=3` 时 `couple()` 已将 `p_current[i]` 平均到 `trace(σ)/(3V)`）：
  ```
  G_{ε,α} · W = 3·V·(P_α − P_target) + mtk_term,
  mtk_term     = 3·mtk_term1 = 6·Σ_{α active} KE_α / (pdim · N)
  ```
  其中 `P_α = (σ_αα + 2·KE_α) / V`，`pdim = Σ_α couple_axis[α]`，`N = (n_dof+3)/3`。
- 桶活塞速度耦合 `fix_nh.cpp::nh_v_press`（`fix_nh.cpp:2564-2569`）：
  ```
  v_α *= exp(−(ε̇_α + mtk_term2) · dt/2),
  mtk_term2 = Σ_{α active} ω̇ / (pdim · N).
  ```
  非启用轴（`couple_axis[α]==false` ⇒ `ε̇_α==0`）仍承担全局 `mtk_term2` 项，恰如 LAMMPS。
- 位置 + 盒子传播采用精确解析漂移因子 `(e^{ε̇_α·dt}−1)/(ε̇_α·dt)`（CPU 端 `expm1c`，零附近用 Taylor `1 + x + 2x²/3` 保证 O(h²) 平滑），`ε̇·dt` 被 clamp 在 ±0.5%/半步避免初始过压态发散。
- 温度斜坡（`fix_nh.cpp::compute_temp_target`）：`δ = clamp((step − nfirst)/(nlast − nfirst), 0, 1)`，`T_live = T_target + δ·(T_stop − T_target)`，每步更新 `kT_` 并按 `Q_t[0]=n_dof·kT·tau_T²; Q_t[k>=1]=kT·tau_T²` 重算原子链质量（baro 链 `Q_b` 不依赖 `n_dof`，但每步亦应按 `kT_` 重算以保持斜坡期间频率不变）。
- 经验阻尼 `tdrag_factor = 1 − drag·dt/tau_T`、`pdrag_factor = 1 − drag·dt/tau_P` 在 NHC 链反向扫描的每个链节都乘入，正向扫描不施加，匹配 LAMMPS `fix_nh.cpp::nhc_temp_integrate` / `nhc_press_integrate`。

单位：`eV / Å / fs / amu / K`，压力 `eV·Å⁻³`，与 `core/units.h` 中常量 `kBoltzmann`、`kMvv2e`、`kFtm2v` 配套。压力符号约定为 LAMMPS `compute pressure` 的 `P = +trace(σ)/(3V)`（正值⇔压缩），与 `cuda::VirialStress` / `frame.scalars->stress` 一致——注意 `cuda::md::Barostat` 当前用的是相反符号 `P = −trace(σ)/3` 仅出于 Berendsen 历史原因，本组件**不**沿用该反号路径，而是与 `cpu::md::NoseHooverNPT` 同号。

## Design

新增公共组件 `ATV::cuda::md::NoseHooverNPT : public cuda::ComponentBase`，结构和职责对齐 `cpu::md::NoseHooverNPT`：

- **Config 结构体**（公开聚合体，逐字段对应 LAMMPS `fix npt` 关键字）：`T_target, tau_T, P_target, tau_P, dt, drag=0.0, mtchain=3, mpchain=3, couple_axis={true,true,true}, T_stop=-1.0`——与 CPU 头部位级一致（`tools/cmp_signatures` 双后端审查会强制要求）。构造函数把 `T_stop<0` 折叠到 `T_target`（无斜坡哨兵）。
- **API 接口**：`name()="nose_hoover_npt"`，`role()=Integrator`，`capabilities()=None`（不直接写 `Energy`/`Force`），`reads()` 暂不在 CUDA `ComponentBase` 中暴露（参见 `barostat.cuh` 中已记录的延后说明），因此 `VirialStress` 自动注入暂走链环 02 的显式 `set_pair`/`set_neighbor` 注册路径，**本规范不引入 reads() 虚函数**；用户/链环 06 须自行 `on(Phase::post_force).add(cuda::VirialStress{})` 或依赖链环 02 的 categorical setter 自动注入。`default_bindings()` 返回 `{initial_integrate@0, final_integrate@0, end_of_step@-40}`，与 CPU 完全一致以保证 `ReduceEnergy @ -50` 之后再跑 P8-P11。`set_run_span(nfirst, nlast)` 返回 `*this`，链式可挂在 `set_integrator` 之前。
- **设备状态（`struct State`，PImpl 风格，匹配 `cuda::md::NoseHoover` 与 `cuda::md::Barostat` 模式）**：
  - `work_t* d_eta_t / d_eta_t_dot / d_eta_b / d_eta_b_dot` — NHC 链上的 η/η̇，长度分别 `mtchain` 和 `mpchain`（`cudaMalloc` 一次，`launch()` 内零拷贝）。
  - `work_t* d_eps_dot` — 三元素，每轴 `ε̇_α`，非启用轴在 `build()` 中初始化为 0 且永不写。
  - `work_t* d_ke_global` — 全局动能标量（重新归约用），`d_ke_axis` — 三元素逐轴动能 `KE_x/KE_y/KE_z`。
  - `work_t* d_diag_axis` — 三元素逐轴压力对角分量缓冲（合并 `σ_αα + 2·KE_α`，应用 `pdim==3` 平均），用于把 `update_eps_dot` 折成一个 host-side 标量装配 + 一个写 `d_eps_dot` 的 `<<<1,1>>>` 更新核。
  - `work_t* d_scale` — 三元素逐轴速度缩放因子 `sx/sy/sz`（`scale_v_by_baro`）。
  - `unsigned int* d_block_count` — last-block 归约同步计数器，沿用 `nose_hoover.cu` 已建立的 fused-reduce 模式。
  - 主机标量：`work_t Q_t[mtchain], Q_b[mpchain], W, W_axis=W/3, n_dof, kT_live, T_live, tdrag, pdrag, KE_global_cache, nfirst, nlast, pdim, couple_axis[3]`。Q/W 在 `build()` 计算，斜坡期间每步在 `compute_temp_target_host` 中重算 `Q_t` 并 `cudaMemcpyAsync` 上传（baro 链 Q_b 同理）。`d_eps_dot[a]` 用主机端 `couple_axis[a]==false` 守卫，避免写 GPU 上独立 mask。
- **核函数布局**（每个内部函数对应一个或若干 kernel，命名空间内 `anonymous namespace`）：
  | helper（对齐 CPU 同名函数 & LAMMPS 方法）   | kernel(s)                                          |
  |--------------------------------------------|----------------------------------------------------|
  | `compute_global_KE` / `compute_axis_KE`    | `reduce_ke_full_k`（fused 全局 + 三轴 4-into-1 归约，复用 last-block 模式） |
  | `nhc_atoms_half`                           | 主机端跑链节循环 → `cudaMemcpyAsync(d_scale, …)` → `scale_vel_global_k`（与 `nose_hoover.cu` 中 `scale_vel` 同形） |
  | `nhc_baro_half`                            | 纯主机标量（链节数 M ≤ 几个），仅写 `d_eps_dot[a] *= baro_scale` 经 `<<<1,3>>>` 或 `cudaMemcpyAsync` |
  | `update_eps_dot`                           | `reduce_stress_diag_k`（从 `frame.scalars->stress` + 逐轴 KE 装配 diag[3]，归约前 `pdim==3` 时取三轴均值） → 主机算 `G_eps_a`、`mtk_term` → `cudaMemcpyAsync(d_eps_dot, …)` |
  | `scale_v_by_baro`                          | 主机算 `sx/sy/sz` → `cudaMemcpyAsync(d_scale, …)` → `scale_vel_per_axis_k`；若 `sx==sy==sz`，复用快路径标量乘 `KE_global_cache` 而不重新归约 |
  | `kick_v_by_force`                          | `kick_v_k`（与现有 VV 风格一致） |
  | `propagate_x_and_box`                      | 主机算 `a[3], drift[3]`（含 `expm1c` Taylor 分支与 ±0.5% clamp）→ `propagate_x_k`（逐轴），随后主机端 `f.box` 缩放（`box[0/4/8] *= a[d]`，与 `barostat.cu` 同模式）+ 必要的 PBC wrap（沿用 `cuda::md::Barostat` 的 inline wrap 等价物） |
  | `extended_hamiltonian`                     | 主机端组装：`KE_global_cache + 已存 f.scalars->energy（D2H 一次）+ ½Σ W_axis·ε̇² + Σ ½Q·η̇² + n_dof·kT·η_t[0] + Σ kT·η + P·V` → `cudaMemcpyAsync(frame.scalars->conserved, …)` |
  非启用轴在主机端短路守卫，相应核启动不会触及该轴的 `d_eps_dot` 或盒子边。

- **生命周期**：`build()` 一次性 `cudaMalloc` 上述设备数组、零初始化 NHC 链、按 `KE_global_cache = 0.5·n_dof·kT` 设置等分初值，重算 `tdrag_factor` / `pdrag_factor`。`launch()` 按 `ctx.stage` 三分支调用上述 helper 序列（初始 P1-P6、final P7、end_of_step P8-P11，**顺序与 CPU 一致**：`nhc_baro_half → nhc_atoms_half → update_eps_dot → scale_v_by_baro → kick_v_by_force → propagate_x_and_box` 在 initial，`kick_v_by_force` 单独在 final，`compute_global_KE → scale_v_by_baro → update_eps_dot → nhc_atoms_half → nhc_baro_half`（反序）+ 末尾把 `KE_global_cache` 与 `extended_hamiltonian` 写回 `frame.scalars->kinetic` / `conserved` 在 end_of_step）。`destroy()` 释放全部 `cudaFree`。
- **精度策略**：与 `cuda::md::NoseHoover` 一致——配置标量与所有缩进 NHC 链上变量用 `work_t`（高精度内部累加），原子坐标/速度按 `real_t`（`ATV_REAL`，通常 `float`），核内显式 `static_cast` 在两者间过渡，避免单精度桶活塞驱动力的灾难性误差累积。
- **MPI**：CUDA 后端目前不跨 rank（单 GPU 假设，沿用 `nose_hoover.cu` / `barostat.cu` 现有姿态）。`compute_global_KE` 在 CUDA 端无需 `MPI_Allreduce`，与 CPU 的 `#ifdef ATV_HAS_MPI` 分支对应位置只写一条注释“CUDA single-GPU; MPI reduce N/A”。

## Files to create or modify

- `include/atomiverse/cuda/md/nose_hoover_npt.cuh` (new) — 公共头：`Config` 聚合体、`NoseHooverNPT` 类声明、`set_run_span`、默认绑定、`role()=Integrator`、PImpl `State` 前向声明。Doxygen 注释要：(1) 同 CPU 头部完整算法块；(2) 标明本组件**替换**而非叠加 `cuda::md::VelocityVerlet + NoseHoover + Barostat`；(3) `pressure-sign` 与 `frame.scalars->stress` / `cuda::VirialStress` 同号说明；(4) LAMMPS `fix_nh.cpp:2541-2569` 行号引用。
- `src/cuda/md/nose_hoover_npt.cu` (new) — 实现：`Config` 构造折叠 `T_stop`，`build()/launch()/destroy()`，全部 kernel 与主机 helper，`expm1c` 的 `__host__` 标记小工具。
- `tests/cuda/md/test_cuda_nose_hoover_npt.cu` (new) — 详细组件级 + 跨后端对照测试，独立翻译单元以隔离 NPT 的初始化噪声不污染 `test_thermostat.cu`。

## Tasks

- [ ] Write failing tests for `cuda::md::NoseHooverNPT` covering iso NPT relaxation, aniso y/z-only coupling, 500→100 K ramp, drag damping, CPU-vs-CUDA parity, and conserved-Hamiltonian flatness (tests/cuda/md/test_cuda_nose_hoover_npt.cu)
- [ ] Implement `cuda::md::NoseHooverNPT` `Config`, class shell, `set_run_span`, `default_bindings`, role tag, and PImpl `State` in include/atomiverse/cuda/md/nose_hoover_npt.cuh
- [ ] Implement `build()` / `destroy()` lifecycle, NHC chain mass init, `tdrag/pdrag` factors, equipartition KE seeding, and all `cudaMalloc`/`cudaFree` in src/cuda/md/nose_hoover_npt.cu
- [ ] Implement device kernels `reduce_ke_full_k` (global + per-axis KE in one pass), `scale_vel_global_k`, `scale_vel_per_axis_k`, `kick_v_k`, `propagate_x_k`, and `reduce_stress_diag_k` in src/cuda/md/nose_hoover_npt.cu
- [ ] Implement host-side Trotter helpers `nhc_atoms_half`, `nhc_baro_half`, `update_eps_dot`, `scale_v_by_baro`, `propagate_x_and_box`, `compute_temp_target`, and `extended_hamiltonian` per LAMMPS-aligned mtk_term1/mtk_term2 formulas with `couple_axis` per-axis guards in src/cuda/md/nose_hoover_npt.cu
- [ ] Implement `launch()` dispatch for `initial_integrate` / `final_integrate` / `end_of_step` stages reproducing the CPU P1–P11 ordering and writing `frame.scalars->kinetic` + `frame.scalars->conserved` at end_of_step in src/cuda/md/nose_hoover_npt.cu
- [ ] Add Doxygen docstrings per project style with `@brief` / `@param`/`@return` (units in fs/K/eV·Å⁻³) plus LAMMPS `fix_nh.cpp:2541-2569` cross-references in include/atomiverse/cuda/md/nose_hoover_npt.cuh and src/cuda/md/nose_hoover_npt.cu
- [ ] Verify CPU↔CUDA parity on an 8-atom LJ block-averaged ⟨T⟩, ⟨P⟩, ⟨box_x⟩ within 3% relative (single-precision noise band) over 2000 steps in tests/cuda/md/test_cuda_nose_hoover_npt.cu
- [ ] Run full check + test suite

## Testing strategy

- **Happy path** — small LJ block (8 atoms cubic, cf. `test_thermostat.cu::make_cubic`), `dt=0.001 fs`, `T_target=100 K`, `tau_T=100 fs`, `P_target=0`, `tau_P=1000 fs`, `couple_axis={true,true,true}`: 2000 步后 `frame.box(0,0)` 收敛到非零有限值，⟨P⟩ 块平均 < 5×10⁻⁵ eV·Å⁻³（≈ 80 bar，stat. tol for 8 atoms），所有坐标/速度有限。
- **Aniso path** — 同样系统，`couple_axis={false,true,true}`，初始 x 方向施加 +1% 应变：x 方向 `box(0,0)` 在整个 run 中保持初值（`ε̇_x≡0`），y/z 各自向 `P_target=0` 弛豫且块平均 `⟨Pyy⟩, ⟨Pzz⟩` 同样 < 5×10⁻⁵ eV·Å⁻³，`⟨Pxx⟩` 可以与之分离（≥ 2× tolerance）证明 x 轴并未被被动耦合。
- **Temperature ramp** — `T_target=500, T_stop=100`，`set_run_span(0, 999)`，run 1000 步：每 100 步抽样 `⟨KE⟩/(½·n_dof·k_B)`，与解析 `T_live(step)=500-400·step/999` 相比偏差 < 8%（在小系统、弱耦合下的允许涨落）。
- **Drag damping** — 同 happy-path 但人为给 `eps_dot[0]=1e-3 /fs` 初值，对比 `drag=0` vs `drag=2`：第一峰后 `box(0,0)` 的振荡振幅在 `drag=2` 路径上至少衰减 2×，确认 LAMMPS-style per-link drag 起效。
- **Edge cases** — `couple_axis={false,false,false}`（全轴刚性，pdim=0）退化为 NVT-only：盒子整步无变化、`ε̇_α≡0`、`mtk_term2 → 0`、`frame.box(d,d)` bit-identical 到 step 0；`T_stop == T_target`（无斜坡哨兵）+ 默认 `set_run_span` 不调用：`T_live ≡ T_target` 整步不变；`mtchain=1, mpchain=1` 仍可运行不 NaN（链长为 1 的退化分支只跑反向扫描 k=0 单步）。
- **Domain validation** — (1) 与 CPU `cpu::md::NoseHooverNPT` 同输入对照：8 原子 LJ，相同 FF/初速度种子/dt/2000 步，块平均 ⟨T⟩、⟨P⟩、⟨box_x⟩ 三量分别在 3% 相对误差以内（单精度噪声允许带，与链 06 基准对齐一致）；(2) 扩展 Hamiltonian 守恒：`frame.scalars->conserved` 在 500-step NVE 等效配置（`tau_T=tau_P=1e6 fs`，弱耦合）下的相对漂移 < 5×10⁻³（单精度 MTK 集成器的合理上界）。

## Out of scope

- **CUDA `cuda::md::MD::set_integrator` 接受 `NoseHooverNPT`** — 走链环 02 (`msu-pe-cuda-components-02-categorical-setter-api`) 的 categorical setter；本环只验证 `sim.add(cuda::md::NoseHooverNPT{…})` 显式注册路径，并在 doctstring 中注明 categorical 调用待链环 02。
- **`reads()==Stress` 虚函数** — `cuda::ComponentBase` 尚未暴露 `reads()`（见 `cuda/md/barostat.cuh:36-42` 已记录的延后），因此 `cuda::VirialStress` 的自动注入仍由用户/链环 02 显式 `on(Phase::post_force).add` 提供；引入 `reads()` 是更大的基础设施改动，单独立项。
- **多 GPU / MPI 归约** — CUDA 后端单 GPU 假设沿用，跨 rank 归约公式与 CPU 已经对齐，但本环不实现 NCCL/MPI 桥接（链环 02-08 全部维持单 GPU）。
- **三斜（triclinic）支持** — CPU 端 `propagate_x_and_box` 走 `box.scale_diagonal(ax, ay, az)`，正交盒子假设。三斜 NPT（off-diagonal `ε̇_xy/yz/xz` 链上变量）是独立的 LAMMPS `fix npt tri` 移植，单独立项。
- **Async / `defer_stream`** — `extended_hamiltonian` 的主机装配涉及一次 D2H 读 `frame.scalars->energy`；Phase H 异步化将统一处理 driver hot-path 的 stream sync，**本环遵循当前同步驱动模型**，不引入 AsyncBridge。
- **`benchmarks/md/polyethylene_tension.cpp` 的 CUDA 分支** — 该 wiring 归 Phase G.9 / 链环 06；本环只保证组件本身可被 `cuda::md::MD` 注册并通过单元测试。
- **替代温度斜坡形态**（非线性、阶梯）— LAMMPS `fix npt` 只支持线性 `T_start → T_stop`，本环保持相同最小契约。

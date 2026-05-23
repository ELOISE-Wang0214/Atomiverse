---
title: MSU PE CUDA components 03 — stage-1 CUDA integrators (NVELimit + FixLangevin)
status: approved
created: 2026-05-23
parent: msu-pe-validation.md (Phase G.5 + G.6)
chain: msu-pe-cuda-components-* (link 3/8)
---

# MSU PE CUDA components 03 — stage-1 CUDA integrators (NVELimit + FixLangevin)

## 综述 (Summary)

本规范是 `msu-pe-cuda-components-*` 链的第 3 / 8 环,目标是把 MSU 非晶 PE 5 阶段协议的 stage 1 所需的两个 CUDA 组件——位移截断积分器 `cuda::md::NVELimit` 与 Langevin 恒温器 `cuda::md::FixLangevin`——从已有的 CPU 等价物 (`cpu::md::NVELimit` / `cpu::md::FixLangevin`) 移植到 CUDA 后端。两者直接对应 `msu-pe-validation.md` 中的 G.5 与 G.6,且是 stage 1 在 CUDA 上能跑通的硬性前提:`msu_pe.data` 是未弛豫的初始构型,初始 `E_pair ≈ 3.4·10⁸ kcal/mol`,普通 Velocity-Verlet 一步即炸——必须有 `nve/limit 0.05` 把每步位移卡死在 0.05 Å;同时 `fix langevin 500 500 10.0 904297` 既是温度控制,又作为唯一的力源驱动体系松弛。两个组件移植完成后,既能通过既有的 `cuda::md::MD` 驱动器 + `add()` 注册路径直接使用,也可在链 02 (categorical setter API parity) 落地后通过 `set_integrator` / `set_thermostat` 调用——后者是软依赖,本规范不阻塞地等待链 02。

## Domain basis

**`nve/limit`** — Velocity-Verlet,带每步欧几里得范数位移截断:每个半步 kick 后,若 `|v|·dt > xmax`,则 `v ← v · (xmax / (|v|·dt))`。等价定义 `vmax = xmax / dt`,`vmax_sq = vmax²`,当 `vsq = vx²+vy²+vz² > vmax_sq` 时以 `s = vmax / sqrt(vsq)` 重标 `v`。对应 LAMMPS `fix nve/limit xmax`(`lammps/src/fix_nve_limit.cpp` L66 `vlimitsq = (xlimit/dtv)^2`,L105–112、L168–175 两半步均做 cap),CPU 实现 `src/cpu/md/nve_limit.cpp` 已就该映射做完逐参数审计 (verdict: ALIGNED)。MSU 协议中 `xmax = 0.05 Å`,`dt = 1.0 fs`。

**`fix langevin`** — Schneider-Stoll Langevin 热浴的力-增量形式。每步在 `post_force` 阶段对每个原子加两项力:

1. 摩擦项 (Stokes drag): `F_drag = -m / (damp · ftm2v) · v`
2. 随机踢力 (FDT-consistent): `F_rand = sqrt(2 m k_B T_live / (damp · dt · mvv2e)) / ftm2v · ξ`,`ξ ~ N(0, 1)`

其中 `T_live = T_start + clamp((step − nfirst)/(nlast − nfirst), 0, 1) · (T_stop − T_start)` 是线性温度坡道 (LAMMPS `fix_langevin.cpp::compute_target()` L794–802、L661)。单位常数 `ftm2v = kFtm2v`,`mvv2e = kMvv2e` 取自 `core/units.h`,`k_B = kBoltzmann`。MSU 协议中 `T_start = T_stop = 500 K`,`damp = 10.0 fs`,`seed = 904297`,`dt = 1.0 fs`。

**RNG 设计 (CUDA-specific)** — CPU 端 `FixLangevin` 用 `(seed, step, atom-tag)` 三元组 splitmix 流生成每原子 IID 的 Box-Muller 高斯,确保跨 rank 数可重现 (同一全局 tag 在不同 rank 切分下得到同一随机力)。CUDA 移植保留同一 `(seed, step, global tag)` 决定性映射:用 counter-based / Philox 或在 kernel 内联展开的 splitmix(每原子从 `f.atoms.tag[i]` 取 tag)生成 Box-Muller 对,**不**使用 `curand_init` 的持久 RNG state(那会把 reproducibility 钉在原子的行索引而非全局 tag,且需要额外的 device 内存)。这样 CUDA 单 GPU 上的随机轨迹与 CPU 单 rank 上的随机轨迹在固定 seed + 固定步号下逐原子可比;`msu-pe-align-03` 规定的跨引擎对比仍按统计分块均值 (Flyvbjerg & Petersen 1989, *J. Chem. Phys.* 91, 461, doi:10.1063/1.457480) 走,不做逐步比对。

参考文献:Verlet (1967) *Phys. Rev.* 159, 98–103;Schneider & Stoll (1978) *Phys. Rev. B* 17, 1302 (doi:10.1103/PhysRevB.17.1302);Tschopp et al. (2010) *Polymer* — MSU 协议出处;LAMMPS 源 `src/fix_nve_limit.cpp`、`src/fix_langevin.cpp` (branch `stable`, fetched 2026-05-20)。

单位:eV, Å, fs, amu, K(`core/units.h`)。

## Design

两个新组件,均落在 `ATV::cuda::md` 命名空间,继承 `cuda::ComponentBase`,匹配 `cuda::md::VelocityVerlet` / `cuda::md::NoseHoover` 的模式:CUDA-free `.cuh` 头声明类,`.cu` 实现 kernel + `launch()`。两者都是无状态-或-轻状态 component,生命周期为 `build()` → `launch()` × N → `destroy()`;`launch()` 内不分配。

**`cuda::md::NVELimit`**

- 构造:`NVELimit(double dt, double xmax)`,镜像 CPU。
- `role() = ComponentRole::Integrator`,`capabilities()` 沿用基类默认 `None`。
- `default_bindings() = { (initial_integrate, 0), (final_integrate, 0) }`,与 `cuda::md::VelocityVerlet` 一致——也意味着它**替代** VV,不与 VV 并存。
- 无持久 device 状态,无需 `build()` / `destroy()` 覆盖(基类空实现)。
- Kernel:`initial_integrate`(半步 kick + cap + 全步 drift + 正交 PBC wrap,按 `cuda::md::velocity_verlet.cu` 的 box[0]/box[4]/box[8] 写法)与 `final_integrate`(半步 kick + cap),按 `ctx.stage` 分派。kick 用 `dtf = 0.5·dt·kFtm2v`,cap 用 `vmax = xmax / dt`、`vmax_sq = vmax²`、`s = vmax / sqrt(vsq)`(代数上等价 LAMMPS `sqrt(vlimitsq/vsq)`,与 CPU 一致)。
- Precision tier 与 VV 一致:per-atom SoA `real_t`,工作精度 `work_t`,box 是 `prec_t`。

**`cuda::md::FixLangevin`**

- 构造:`FixLangevin(double T_start, double T_stop, double damp, double dt, std::uint32_t seed)`。`T_stop < 0` 折叠到 `T_start`(no-ramp 哨兵,与 CPU 一致)。
- `set_run_span(nfirst, nlast)` 链式 setter,返回 `*this`;`T_target() const` 暴露当前 live 设定值供 observer 读取。
- `role() = ComponentRole::Thermostat`,`capabilities() = Capability::Force`(同 CPU——FixLangevin 写 fx/fy/fz,让一个**只有** Langevin 的 driver 也满足 `DynamicsDriver` 的 "≥ 1 个 Force-capable component" 校验)。
- `default_bindings() = { (post_force, +500) }`——在所有物理力评估器之后、用户自定义 post-force 钩子之前累加。
- 无持久 device 状态(per-atom RNG 在 kernel 内从 `(seed, step, tag[i])` 重建,无 curand_init / 无 state buffer)。Build / destroy 无开销。
- Kernel `apply_langevin`:每原子计算 `gfactor1 = -m / (damp·kFtm2v)`、`gfactor2 = sqrt(2·kBoltzmann·T_live·m / (damp·dt·kMvv2e)) / kFtm2v`,从 `(base, tag[i])` 派生 64-bit seed,设备端 splitmix 生 6 个 u01,Box-Muller 产 3 个高斯 `rx, ry, rz`,累加 `fx[i] += gfactor1·vx[i] + gfactor2·rx`(y/z 同)。`T_live` 在 host 侧按 `(ctx.step, nfirst_, nlast_)` 计算一次后作为 kernel 标量参数传入,避免在 kernel 内重复条件分支。`T_live_` 成员同步更新,供 `T_target()` 查询。

无新增 frame 字段,无新增 BuildContext / LaunchContext 字段——两个组件都从 `ctx.frame.atoms.{x,y,z,vx,vy,vz,fx,fy,fz,mass,tag}`、`ctx.frame.box`、`ctx.step`、`ctx.runtime.stream` 读现有数据。生命周期与 VV/NoseHoover 完全一致,故 `cuda::md::MD` 接受它们走 `add()` 路径无须改 driver。链 02 落地后再补 `set_integrator(NVELimit{...})` / `set_thermostat(FixLangevin{...})` 派发,本环不阻塞等待。

## Files to create or modify

- `include/atomiverse/cuda/md/nve_limit.cuh` (new) — `NVELimit` 类声明
- `src/cuda/md/nve_limit.cu` (new) — kernels + `launch()` 派发
- `include/atomiverse/cuda/md/fix_langevin.cuh` (new) — `FixLangevin` 类声明 + `set_run_span` + `T_target`
- `src/cuda/md/fix_langevin.cu` (new) — kernels + RNG inline + `launch()`
- `tests/cuda/md/test_nve_limit.cu` (new) — 单组件单测
- `tests/cuda/md/test_fix_langevin.cu` (new) — 单组件单测
- `tests/CMakeLists.txt` — 注册 `test_cuda_md_nve_limit` 与 `test_cuda_md_fix_langevin` 两个 CTest 目标
- `src/cuda/CMakeLists.txt` — 在 atomiverse CUDA 源列表中加入两个新 `.cu` 文件
- `.claude/specs/INDEX.md` — 增加本规范条目

## Tasks

- [ ] Write failing tests for `cuda::md::NVELimit` (tests/cuda/md/test_nve_limit.cu): cap clamps a synthetic single-atom step to ≤ xmax; NaN-free after 100 steps on a configuration with 1e8 eV/Å seed forces; matches CPU `nve_limit` to single-precision tolerance on a 64-atom no-force drift over 50 steps
- [ ] Implement `cuda::md::NVELimit` in include/atomiverse/cuda/md/nve_limit.cuh + src/cuda/md/nve_limit.cu (initial/final kernels with cap, PBC wrap mirroring VelocityVerlet)
- [ ] Write failing tests for `cuda::md::FixLangevin` (tests/cuda/md/test_fix_langevin.cu): kinetic temperature converges to T_target within ±10% over 5000 steps on a 256-atom synthetic system at fixed seed; linear T ramp `T_start=200 → T_stop=500` honoured at `step = nlast/2`; per-atom force finite; component contract (`role()==Thermostat`, `capabilities()==Force`, default_bindings on `post_force @ +500`)
- [ ] Implement `cuda::md::FixLangevin` in include/atomiverse/cuda/md/fix_langevin.cuh + src/cuda/md/fix_langevin.cu (host-side T_live update + apply_langevin kernel with inline splitmix Box-Muller seeded by `(seed, step, tag)`)
- [ ] Add doxygen `@brief / @param (with units) / @return` per CLAUDE.md doc.style on every new public symbol in both new headers, mirroring the LAMMPS-correspondence audit tables from the CPU headers verbatim
- [ ] Register the two new `.cu` files in src/cuda/CMakeLists.txt and add `test_cuda_md_nve_limit` + `test_cuda_md_fix_langevin` to tests/CMakeLists.txt; add an INDEX.md entry for `msu-pe-cuda-components-03-cuda-integrators-stage1`
- [ ] Run full check + test suite (cmake build + ctest -R 'cuda_md_(nve_limit|fix_langevin)' then the full ctest run)

## Testing strategy

Happy path:

- `NVELimit` smoke: 64-atom lattice, dt = 1.0 fs, xmax = 0.05 Å, zero forces — every atom drifts at its initial v; no atom moves more than 0.05 Å in any step.
- `FixLangevin` smoke: 64-atom system at T = 0, no LJ, `T_start = T_stop = 300 K`, fixed seed; after 5000 steps the per-atom kinetic temperature mean is within ±10% of 300 K (Langevin equilibration tolerance, large-N statistical limit).

Edge cases:

- `NVELimit` cap engages: synthetic single atom with `f` chosen so that the half-step kick gives `|v|·dt = 10·xmax`; after `initial_integrate`, the resulting per-step displacement equals xmax exactly to single-precision tolerance.
- `NVELimit` NaN guard: 16-atom configuration with `fx[0] = 1e8 eV/Å` (mimicking msu_pe.data overlap force); 100 steps; no atom's position/velocity goes non-finite.
- `FixLangevin` ramp: with `T_start = 200, T_stop = 500, nfirst = 0, nlast = 1000`, at `ctx.step = 500` the host-computed `T_live_` equals 350 K (read via `T_target()`); at `step ≥ 1000` it stays at 500 K (clamped).
- `FixLangevin` Force-capability advertised: `comp.capabilities() & Capability::Force` is non-zero so a `cuda::md::MD` with only `NVELimit + FixLangevin` (no LJ) passes `validate_components()`.

Domain validation (per `$META.science.required`):

- `NVELimit` algebraic equivalence to LAMMPS `fix nve/limit`: with a fixed initial state and a fixed `f`, one step of `cuda::md::NVELimit{dt, xmax}` produces positions/velocities matching `cpu::md::NVELimit{dt, xmax}` (which already audits ALIGNED vs LAMMPS) to `1e-5` in single-precision (CUDA `real_t = float`).
- `FixLangevin` FDT check: ensemble-averaged ⟨T⟩ over the last half of a 10 000-step run at `T_target = 300 K` matches 300 K within ±5% (target tolerance — half-precision Langevin equilibration on 256 atoms saturates at this band per `msu-pe-align-03`).

Standard CUDA component contract (both):

- `role()`, `capabilities()`, `default_bindings()` return the documented values.
- `build()` → `launch()` × N → `destroy()` round-trips without `cudaGetLastError()` raising.

## Out of scope

- LAMMPS `fix langevin` keywords `scale` (per-type friction), `tally` (reservoir-energy bookkeeping), `zero` (COM-drift removal), `angmom`/`omega` (rotational DOF). CPU `FixLangevin` already lists these as deferred and they remain deferred on CUDA — the MSU deck does not use them, and `zero` is handled by the upcoming `cuda::md::ZeroMomentum` (G.1, separate spec).
- `set_integrator` / `set_thermostat` categorical setters for these two components — hard-dependency on link 02 (`msu-pe-cuda-components-02-cuda-categorical-setters`); registration via raw `add()` is functionally complete without it.
- Multi-GPU / CUDA+MPI variants. Per `msu-pe-validation` fixed scope decision 1 the CUDA arm is single-GPU.
- Stage 2–5 CUDA components (`NoseHooverNPT`, special-bonds exclusions, `fix momentum`) — those are G.7/G.8/G.1 and live in later links of this chain.

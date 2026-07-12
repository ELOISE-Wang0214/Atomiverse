---
title: MD constraints (SETTLE, then RATTLE) — stress-test the phase model, the write-set, and Ndof
slug: md-constraints
status: approved
science_required: true
created: 2026-07-12
acceptance: md-constraints.acceptance.md
depends_on: regressions-harness
conflicts_with: manybody-potentials
---

# 约束 — 压测相位模型 / 写集 / 自由度

## Summary

广度优先探针 #2。目标是撞 MD 相位模型和 Component 标签系统。**它撞出的三个缺陷全都
不是"约束的问题"，而是架构的问题** —— 这正是探针的价值。

> `conflicts_with: manybody-potentials` —— 两者都要动 `comm/halo`。**串行做，不要并行改。**

## 侦察结论（已核实）

**零命中**：SHAKE/RATTLE/LINCS/SETTLE/刚体/虚位点/冻结原子，全树没有任何实现，也没有
任何 group / 原子子集选择机制。可复用的是 `ResolveView`（`cpu/potentials/resolve_view.h:40-78`，
全局 tag → 本地行）与 `BondedTupleCache`（`cpu/potentials/bonded_cache.h:43-85`，
带 migration/halo 失效代的元组缓存）—— 二者形状正好，但**住在 `potentials/` 里**，
约束不是力场求值器，跨域 include 违反 CLAUDE.md 分层。**要先提升到共享位置。**

**相位模型够用，不需要新相位。** `post_integrate`（phase 1）是**空的**，且
`facade/stage.h:32` 与 `docs/developer/md-execution-model.md:65` 早就写明它就是留给约束的。
位置约束 → `post_integrate`；速度约束 → `final_integrate @ +200`（在 NoseHoover +100 之后）。

**三个真缺陷：**

- **(A) `Capability` 枚举有洞 —— 没有 Position / Velocity 位。**
  `facade/component_role.h:63-74` 只有 `Energy|Force|Stress|Charge|Dipole|DensityMatrix|Orbitals|Hessian`，
  而其文档注释（`:7-9`）说 Capability 就是"这个组件能写 Frame 的哪些字段"。
  于是 **VelocityVerlet / NoseHoover / Deform / ZeroLinearMomentum 这些真正在写 x/v 的
  组件，全部声明 `Capability::None`** —— 文档与行为矛盾。
  **必须先裁决：Capability 是"Evaluator 产出什么"（那就修文档，约束用 None），
  还是真正的写集（那就加 Position/Velocity 位，并修正上述四个组件）。本 spec 选后者**
  —— 因为约束是时间步里最具侵入性的组件，让它宣称"什么都不写"是个谎言，而写集是
  未来做组件冲突检测（两个组件同时写 v）的唯一基础。
- **(B) `ComponentRole` 没有 `Constraint`。** 若拿 `Integrator` 冒充，
  `DynamicsDriver::validate_components`（`cpu/driver/dynamics_driver.h:10-13`，要求
  ≥1 个 Integrator）会被骗过 —— 一个**只有约束、没有积分器**的 Plan 能通过 build 校验。
  静默物理错误。加一个枚举项。
- **(C) Ndof 硬编码 `3N−3`，散在 5 处，无发布通道**：
  `thermo.cpp:23,57`、`init.cpp:183`、`nose_hoover.cpp:52`、`fix_csvr.cpp:62`、
  `nose_hoover_npt.cpp:125`。约束会减少自由度（刚性水：每分子 −3）。不修则刚性水 NVT
  的温度报低约 33%，且恒温器照着错的设定点驱动。**且有构建顺序陷阱**：恒温器在自己的
  `build()` 里缓存 Ndof，组件按 Plan 顺序 build —— 所以自由度必须在
  `install_defaults()`/`setup_pass` 阶段发布，不能在约束的 `build()` 里减。

**MPI 的出路（不必新增通信轮次）**：走 LAMMPS `FixShake` 的**簇复制**路线 ——
每个 cluster 在所有拥有其任一原子的 rank 上冗余求解，用
`HaloParticipant::comm_forward_fields` 把自有的 `x_unc` / `v_unc` 影子数组挂到**现有**
halo 消息上（`halo_participant.h:61-74`）。SETTLE 是解析的、**无迭代**，完美契合；
迭代版 RATTLE 在完全复制的 cluster 上**本地**迭代，循环内不通信。
（注意：`CommField::data` 必须是组件自有、`build()` 时按 `atoms.nmax` 定尺寸、全程稳定的
指针 —— 不能指向 `frame.atoms.vx.data()`，因为 `Atom::grow()` 会重分配。）

**顺序陷阱**：VelocityVerlet 把 half-kick 和 drift **融合在一个 launch 里**
（`velocity_verlet.cpp:31-39`），且其后立刻 `wrap_orthorhombic`。所以位置约束看到的
`x` 已被 wrap、`v` 已被 NH 缩放。约束必须自己在 wrap **之前**留下未约束坐标的影子副本，
或对每个 cluster 键做 min-image。

**语义分歧**：bonded 力对无法解析的 tag 是**静默丢弃**（`bond_harmonic.cpp:57`）。
约束**必须抛异常** —— 丢一个 cluster = 一根没被约束的键 = 物理污染，不是近似。

## Design

先 **SETTLE**（解析、三原子刚性水、无迭代 → 验证相位/写集/Ndof 三条链路），
再 **RATTLE**（迭代 → 压测"迭代 + 分布式"）。

## Files

- `include/atomiverse/facade/component_role.h` — `Capability::Position|Velocity`；`ComponentRole::Constraint`
- `cpu/md/velocity_verlet.h`、`nose_hoover.h`、`deform.h`、`init.h` — 补正确的 `capabilities()`
- Ndof 单一所有者（`ctx.md.n_dof_global`，在 `install_defaults`/`setup_pass` 发布）
  + 5 个消费点改为读取它
- `resolve_view.h` / `bonded_cache.h` 提升出 `potentials/`（共享位置，见 CLAUDE.md 分层规则）
- `include/atomiverse/cpu/md/constraint_cluster.h` — cluster 发现 + 元组缓存（抛，不丢）
- `include/atomiverse/cpu/md/fix_settle.h` + `src/…`
- `include/atomiverse/cpu/md/fix_rattle.h` + `src/…`
- `cpu/md/md.h` — `set_constraint(...)`（可多个：水用 SETTLE，溶质用 RATTLE）
- `regressions/md/settle_water/`

命名：LAMMPS 式（`FixLangevin`/`FixCSVR` 先例）→ `cpu::md::FixSettle`、`cpu::md::FixRattle`。
**不要叫 `FixShake`** —— roadmap 已把 SHAKE 标为 ⊘（被 RATTLE/P-LINCS 上位替代）。

## Tasks

- [ ] **裁决 Capability 的语义**并落地（加 Position/Velocity 位 + 修正 4 个既有组件的标签）
- [ ] `ComponentRole::Constraint` + `DynamicsDriver` 校验（有约束但无积分器 → 抛）
- [ ] Ndof 单一所有者 + 5 个消费点迁移 + 构建顺序保证（约束的 DOF 扣减先于恒温器 build）
- [ ] `resolve_view` / `bonded_cache` 提升出 `potentials/`
- [ ] cluster 发现（从拓扑 + 类型筛选）；无法解析 → **抛**
- [ ] `FixSettle`（解析三原子）+ 未约束坐标影子数组（绕开 wrap 陷阱）
- [ ] `FixRattle`（迭代到容差）
- [ ] `CommField` 携带 `x_unc`/`v_unc` → 多 rank 与单 rank 一致
- [ ] regressions 用例
- [ ] 架构结论写回 `.claude/NOTES.md`

## Testing

- 单测：cluster 发现；SETTLE 解析解 vs 手算；Ndof 扣减正确；无积分器时抛。
- **不变量测试（约束的真判据）**：约束键长在整条轨迹上恒定到 1e-10；速度约束后
  `d(r²)/dt = 0`（即 `(v_i − v_j)·(r_i − r_j) = 0`）。
- 刚性水 NVT 的**温度必须落在设定点上** —— 这是 Ndof 修复的直接判据（不修则低 ~33%）。

### 三层测试的分工（CLAUDE.md 硬规则，本 spec 一律遵守）

| 层 | 跑多久 | 测什么 | 硬禁 |
|---|---|---|---|
| `tests/` | **整套秒级**（快速迭代循环） | **单个函数 / 方法** | **严禁 e2e** |
| `regressions/` | 单例 ≤ 60 s | **库外 e2e，对着文献** | 不进默认 ctest |
| `benchmarks/` | 验收时 | 对外部引擎 parity + 测速 | 不进 ctest |

本 spec 的每个交付功能，**必须在 `regressions/` 留下恰好一个教科书级用例**；
e2e **绝不允许**塞进 `tests/`。

### 风格（CLAUDE.md 硬规则）

- **OOP，不是自由函数汤**：查表即构造函数，**禁止** `make_*` / `builtin_*` /
  `resolve_*` / `parse_*` 工厂函数。
- **只运行一次的代码不要抽成函数** —— 就地内联在它运行的地方。函数是为了被调用
  ≥ 2 次（或被单独测试）才存在，不是为了给一段话取个名字。**在第二个调用者出现时
  才抽取。**

## Out of scope

- LINCS/P-LINCS（roadmap：LINCS 被 P-LINCS 上位替代；先做 RATTLE/SETTLE）
- 刚体 / 虚位点 / 氢质量重分配
- CUDA 镜像

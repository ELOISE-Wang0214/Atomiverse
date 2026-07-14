---
title: GCMC + NPT-MC — stress-test the Frame's fixed-N assumption and the MC transaction
slug: mc-gcmc-npt
status: code-complete
science_required: true
created: 2026-07-12
acceptance: mc-gcmc-npt.acceptance.md
depends_on: regressions-harness
---

# GCMC + NPT-MC — 压测 Frame 的固定 N 假设

## Summary

广度优先探针 #4。假设是"Frame 假定 N 固定，GCMC 会撞破它"。

**侦察推翻了这个假设，而这正是最有价值的发现：`comm/` 早就把变长 N 解决了。
本 spec 的第一原则是——不要重新发明它，把它提升出来复用。**

## 侦察结论（已核实）

**变长 N 已经存在，只是被墙起来了：**

- `cpu::Atom` **本来就可变长**：`grow(int)`（`frame/atom.h:117`，LAMMPS DELTA 增长，
  零填充尾部）与 `set_capacity(int)`（`:131`，**可增可缩**）。**N 归 `Atom` 所有，不归
  `Frame`** —— `Frame` 只是个薄持有者，它的"固定 N"只存在于构造函数里。
- `comm::pack_exchange` / `unpack_exchange` / `max_exchange_doubles`
  （`cpu/comm/atom_pack.h:112,119,127`）**已经能把一个原子的完整行序列化进 double 数组
  ——含 x/v/f/mass/charge/type/tag/image + 键/角/二面角 + special 列表 —— 而且是公开的、
  零 `#ifdef ATV_HAS_MPI`。这三个函数直接就是 GCMC 事务的载荷。**
- MPI migration **本来就在做插入与删除**（`src/cpu/comm/atom_migration.cpp`）：
  插入 = grow → 在 `nlocal` 处写新行 → `++nlocal`（`:264-282`）；
  删除 = compact-survivors-forward（`:222-234`）配合 `move_row(frame, ni, oi)`（`:134-180`，
  逐字段拷贝整行含拓扑）。
- **邻居表也早就准备好了**：`NeighborList::launch` 有个 **非 LAMMPS 的 `size_stale`
  触发器**（`neighbor/neighbor_list.h:149`，注释明说是为 migration 的变长 N 加的），
  N 一变就无条件全重建。`CellList::bin_atoms` 也按需自增容量。**插入/删除不需要碰邻居表。**

⇒ **`move_row` 现在是 `static` + `#ifdef ATV_HAS_MPI`（`atom_migration.cpp:39,134`），
外部够不着。必须提升为公开的 `cpu::Atom::copy_row(int dst, int src)`，并让
`atom_migration.cpp` 改调它（单一定义，不许复制代码到 `mc/`）。**

**真正破的地方（2026-07-14 对照 driver-handler 重构后的树刷新）：**

- **(A) `mc::Transaction` 只是个位置差分**（`cpu/mc/transaction.h:32-38`：`moved` +
  `old_xyz` + `new_xyz` 七个平行数组），**没有动词**。拒绝一次插入 → `nlocal` 已加上
  收不回；拒绝一次删除 → 原子救不回；box 改动完全没记录。`MC::rollback_`
  （`cpu/mc/mc.h:190-197`）只会还原 x/y/z。

  **用户指令（2026-07-12）：「我们的输入输出尽可能是 Frame，所以不应该有 old_xyz
  *_xyz 这种东西。」**

  ⇒ **不要给这个错的设计加动词。把 `old_x/old_y/old_z/new_*` 七个平行 `std::vector`
  整体删掉，改成 Frame 语义。** `Transaction` 记录的是一个 **Frame 补丁**：被这次
  trial 触碰到的原子所构成的一个小 `cpu::Frame`（连同 `Box`）。回滚 = 把这个补丁
  Frame 拼回主 Frame。

  这一步把四个动词统一成**同一件事**——"这里有一个要还原的 Frame"：
  - Displace：补丁 = 被移动原子的原始行；
  - Delete：补丁 = 被删原子的完整行（**含拓扑**，这正是裸 xyz 数组永远做不到的）；
  - Insert：补丁为空 + 记录新增的行数（回滚 = 截断 `nlocal`）；
  - Box：补丁携带原 `Box`。

  行的搬进搬出复用**提升后的 `Atom::copy_row`**，序列化（若需要）复用已有的
  `comm::pack_exchange`。**`mc/` 里不再出现任何裸坐标数组。**
- **(B) ~~acceptance seam~~ — driver-handler 重构已经兑现了它。** 写死的
  `metropolis_accept_` 和整个 `SamplingDriver` 类都已删除；seam 现在叫
  **`OutcomePolicy`**（`cpu/driver/outcome_policy.h`），`Metropolis` 是它的 NVT 实现，
  换接受准则 = 一个 `set_outcome()` 调用。**剩余工作只是两个新 policy 实现**：
  μVT 需要 N/V/μ/Λ，NPT 需要 `P·ΔV − N k_BT ln(V'/V)` —— 二者都能从
  `LoopState.frame` 读到 N 和 V，不需要动 seam 本身。注意 `MC` preset 把
  `metropolis_&` 存成成员来回读 `current_energy()`（`mc/mc.h:212`）；新 preset
  同样各自持有自己 policy 的引用。RNG 走 `ctx.mc.rng`（与 Mover 共用一条流，
  可复现性 load-bearing —— 见 `outcome_policy.h:79-88` 的注释）。
- **(C) 没有全局 tag 分配器**：`grow()` / `set_capacity()` 给的是 identity tag
  （`atom.cpp:71-78`、`:151-154`）。删除压缩 + 插入之后，新行会继承落点行的
  **陈旧 tag** → tag 重复。需要一个"下一个可用 tag"计数器（LAMMPS 的 `tag_extend`）。
- **(D) `units.h` 里没有普朗克常数** → 热德布罗意波长 Λ 算不出来（GCMC 接受率必需）。
- **(E) `CellList::grid_stale()` 只比较 x 方向边长**（`cell_list.cpp:131`，
  `box_at_last_build_` 本身就只是一个标量）→ **各向异性**体积移动会静默丢 pair
  （该失效模式的注释就在 `cell_list.cpp:253-266`，曾导致 NPT 炸掉）。
  各向异性体积移动之前必须修。
- **(F) 没有任何单粒子能量通路**：`PairLJCut` 只写聚合标量，没有 `eatom[]`、没有
  `energy_of(i)`。⇒ 每次 trial 都是 **O(N) 全量重算**。这是 `MC` 的**刻意契约**
  （`mc/mc.h:96-99` 的 run() 文档：全量重算是"对所有势都成立的唯一形式，含长程/kspace"）。
  **v1 接受它**（正确但慢）；`PairEvaluator::energy_of_rows()` 是后续性能 spec。

**架构判断：GCMC 不需要新的 Driver category —— 这一点已被 driver-handler 重构
彻底坐实：category 本身都没了。** 终止语义没变（还是 trial-move 预算）。GCMC/NPT-MC
各是一个 **preset**（构造函数里装上自己的 mover + `OutcomePolicy`，与 `cpu::mc::MC`
同型，直接派生 `CPUDriver<sampling::Phase>`）；是派生 `MC` 还是与它并列，实现时定，
不在 spec 层面锁死。

**v1 单 rank。** `cpu::mc::MC::install_defaults`（`mc/mc.h:160-166`）从不安装 DD，
`ctx.md.dd` 在 MC 下是 null。多 rank GCMC 需要全局 tag 分配器 + DD，留给后续。

## Design

先 **NPT-MC**（体积移动，不改 N —— 先把 Transaction 的 verb 和 acceptance seam 跑通），
再 **GCMC**（插入/删除 —— 压变长 N）。

多个 mover 并存有个真实的 API 问题：`MC::validate_components`（`mc/mc.h:169-182`）
只要求 ≥1 个 Mover，而所有绑在 `propose` 槽的 mover **每个 sweep 都会全部触发**。
所以需要一个**选择器 mover**（按概率抽本轮做哪种移动），而不是把 `TranslateMove` 和
`InsertDeleteMove` 并排注册。这是本 spec 最大的 API 设计决策。

注意 Frame-patch 重写会连带改 `MC::rollback_` 与 run() 里
`frame.scalars.energy = metropolis_.current_energy()` 的恢复逻辑（`mc/mc.h:125-129`）
—— 回滚从"还原 x/y/z"变成"把补丁 Frame 拼回去"。

## Files

- `frame/atom.h` / `atom.cpp` — **提升** `copy_row(dst, src)`；`atom_migration.cpp` 改调它
- `frame/atom.h` — 全局 tag 分配器（下一个可用 tag）
- `core/units.h` — 普朗克常数 / `lambda_thermal(mass, T)`
- `cpu/mc/transaction.h` — **重写为 Frame 补丁**（删掉 `old_xyz`/`new_xyz` 平行数组）
- `cpu/driver/outcome_policy.h`（或 `cpu/mc/` 下新文件）— **μVT / NPT 两个
  `OutcomePolicy` 实现**（seam 已由 driver-handler 兑现，`Metropolis` 是现成范本）
- `cpu/mc/volume_move.h`、`cpu/mc/insert_delete_move.h`、`cpu/mc/move_selector.h` + `src/…`
- `cpu/mc/npt_mc.h`、`cpu/mc/gcmc.h`（preset：构造函数装 mover + outcome policy，
  与 `cpu::mc::MC` 同型，持 μ/T/P 状态）
- `cpu/neighbor/cell_list.cpp` — `grid_stale` 检查三个轴，不只 x
- `regressions/mc/gcmc_lj/`

## Tasks

- [x] `Atom::copy_row` 提升 + `atom_migration.cpp` 去重（**不许复制代码**）
- [x] `Transaction` 重写为 **Frame 补丁**（删除 `old_xyz`/`new_xyz`）+ `MC::rollback_` 全面支持
- [x] μVT / NPT `OutcomePolicy` 实现（seam 已存在 —— `set_outcome()` 即插即用）
- [x] 全局 tag 分配器（删除+插入后 tag 唯一性单测）
- [x] `units.h` 普朗克常数 + Λ
- [x] `grid_stale` 修成三轴
- [x] `VolumeMove` + `NptMC`（先跑通 Box 回滚）
- [x] `InsertDeleteMove` + `Gcmc`（变长 N）
- [x] `MoveSelector`（多 mover 的正确语义）
- [x] regressions 用例
- [x] 架构结论写回 `.claude/NOTES.md`

## Testing

- 单测：插入后拒绝 → Frame 逐字段还原（含 nlocal/tag/拓扑）；删除后拒绝 → 原子完整复活；
  体积移动拒绝 → box **与**坐标都还原；tag 唯一性；Λ 对着解析式。
- **细致平衡是 MC 的真判据**：接受率随 μ 的响应必须满足 μVT 的解析关系。

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

- 多 rank GCMC（需要全局 tag 分配器 + DD；v1 单 rank）
- 单粒子/子集能量优化（`energy_of_rows()` 是后续性能 spec）
- CBMC / 分子（非单原子）插入 —— 分子插入还需重跑 `add_bond`/`build_special_list`
- Gibbs 系综 / 反应系综

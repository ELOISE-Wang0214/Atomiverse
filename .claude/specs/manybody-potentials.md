---
title: Many-body potentials (EAM, then Stillinger-Weber/Tersoff) — stress-test the pair/neighbor abstraction
slug: manybody-potentials
status: approved
science_required: true
created: 2026-07-12
acceptance: manybody-potentials.acceptance.md
depends_on: regressions-harness
---

# 多体势 — 压测 pair/neighbor 抽象

## Summary

广度优先探针 #1。目标**不是** EAM 本身，而是用 EAM 去撞 `potentials/` + `neighbor/`
+ `comm/` 的抽象边界：现有的整个 pair 体系只被 **LJ 一个最温顺的两体势**验证过。

侦察结论（file:line 见下）：**抽象在 `potentials/` 层没破，破在 `neighbor/` 和
`comm/`。** 这两处的修复就是本 spec 的真正交付物。

## 侦察结论（已核实）

**没有任何多体/多遍力代码存在**（EAM/Tersoff/SW/嵌入密度，四个后端全零命中）。
不存在可复用的多体骨架 —— 但存在三个**关键的可复用/可提升资产**：

1. **`HaloParticipant::CommField`**（`cpu/component/halo_participant.h:61-74`）——
   声明式的 per-atom 额外 halo 载荷，`comm_forward_fields()` / `comm_reverse_fields()`，
   引擎侧打包已实现（`halo_exchange.h:178-181, :271`，绑定于 `cpu_driver.h:186-203`）。
   **这正是 EAM 的 ρᵢ 需要的机制，而它在生产代码里从未被使用过** —— 全树唯一实现者是
   测试桩 `tests/cpu/test_install_infra_hook.cpp:142`。EAM 是它的第一个真实用户。
2. **全邻居表已写好但从未被调用**：`CellList::build_full` / `refresh_full`
   （`cpu/neighbor/cell_list.h:100-104, :131-138`；实现 `cell_list.cpp:788, :947`）。
   `NeighborList` 只调 `refresh_half`。
3. `BondedForceScratch` / `Voigt6` / `ResolveView` / `BondedTupleCache` 均可直接复用。

**破的地方，精确定位：**

- **(A) 只发布半表。** Tersoff/SW 的三体角需要 atom i 的**全部**邻居；半表里有一半
  邻居挂在别人行上。CPU 侧没有请求开关 —— 而 **CUDA 侧已经有这个模式**
  （`cuda::DynamicsServices::flat_pairs_requested`，`cuda/component/context.cuh:51`），
  且 CUDA 的 PairList 本来就是全表（`cuda/neighbor/neighbor_list.cuh:50-83`）。
  **CPU 落后于 CUDA。**
- **(B) halo 每步只有一次 forward（`pre_neighbor`）+ 一次 reverse（`post_force`），
  且 reverse 载荷被显式审计死为 fx/fy/fz**（`cpu/comm/halo_exchange.h:85-109` 的注释
  逐字预言了这个场景："Any future kernel that writes to additional ghost-row fields
  must either extend this payload list … or operate on owners only"）。
  EAM 需要在 **force 相位中间**通信 ρᵢ —— **这个时间点不存在**。
- **(C) `PairKernel`/`SumKernel`（`potentials/pair_kernel.h`）只建模 `fpair(r²)`，
  结构上无法表达多体**（无 per-atom 状态、无三元组）。EAM 不是它的扩展，是它的兄弟。
- **(D) `Pair` 基类不合身**：`launch()` 与 `role()` 均为 `final`（`cpu/potentials/pair.h:30,38`），
  `set_pair` 是**单例槽**（`cpu_driver.h:265-269`，第二个 pair 风格会抛）。

## Design

### 分类：走 CoulLong 的路，不要塞进 Pair

`CoulLong`（`cpu/potentials/coul_long.h:14-17`）已经确立了范式：一个力风格如果不是
`fpair(r²)`，就做成裸的 `ComponentBase + NeighborConsumer`，用
`.on(md::Phase::force).add(comp, order)` 挂载，而不是 `Pair` 子类 + `set_pair`。
EAM/SW/Tersoff 全部照此办理。

### (A) 全邻居表：CPU 补齐 CUDA 已有的请求开关

`cpu::DynamicsServices`（`cpu/component/context.h:48-68`）增 `bool full_list_requested`
（镜像 CUDA 的 `flat_pairs_requested`）。消费者在 `build()` 里置位；`NeighborList`
据此改调 `refresh_full`。**这是把死代码 `refresh_full` 接上电，不是新写邻居算法。**

### (B) ρᵢ 的通信 —— 本 spec 的核心架构决策

EAM 的两遍：pass1 累加 ρᵢ（owner 上只是**部分和**，因为半表；ghost 行上也有部分和），
pass2 需要 `F'(ρ_j)`，而 j 可能是 ghost。三条路，**spec 选 (B2)**：

- **(B1) LAMMPS 经典**：ghost 的 ρ 部分和 reverse 回 owner → owner 补全 → `F'(ρ)`
  forward 出去。需要在 force 相位中间插入**两次**额外通信。破坏"每步一次 forward /
  一次 reverse"的不变量，且 `HaloSend::launch` 不是可重入的。
- **(B2) 全表 + 加宽 halo（选定）**：用 (A) 的全表，令 halo cutoff = `2 × r_cut`
  （ghost-of-ghost），则**每个 ghost 行的 ρ 都能在本地独立算全**，pass2 直接可用。
  **完全消除 force 中段通信**，代价是更厚的 halo。这条路只依赖 (A)，不动 comm 的
  时序不变量 —— 对一个探针而言这是正确的风险取舍。
- **(B3)** 让 EAM 声明 `CommField{rho}` 走 forward —— **错误**：forward 在
  `pre_neighbor` 触发，会把**上一步**的 ρ 送出去，静默错一个时间步。

**但 `CommField` 仍要在本 spec 里被点亮**：EAM 的 `F'(ρ_i)` 作为 forward field 走
声明式通道（在 (B2) 下它是可选的优化，在多 rank 校验里作为交叉验证）。这样这个
"设计好却从未被使用"的协议第一次拿到生产用户，其真实缺陷才会暴露。

### 数据表：静态 .h，禁止解析器

CLAUDE.md 的 GPW 先例是硬约束：**没有解析器，没有 `builtin_*`/`make_*` 工厂，
构建/测试/运行期不读任何数据文件**。因此：

- `cpu/potentials/eam_data/eam_tables.h`（静态 C++ 表，一次性转录，%.17g）
- 生成器 `tools/gen_eam_tables.cpp` 保留在树内（GPW 的 `tools/gen_gpw_tables.cpp` 先例）
- 查找构造函数 `EAMPotential(Element, name)`，未知元素/名字**抛异常**（无 fallback）
- **禁止**运行时 setfl/`.eam.alloy`/`.tersoff` 解析器

## Files

- `include/atomiverse/cpu/component/context.h` — `DynamicsServices.full_list_requested`
- `include/atomiverse/cpu/neighbor/neighbor_list.h` / `src/…` — 据 flag 调 `refresh_full`
- `include/atomiverse/cpu/potentials/eam.h` + `src/cpu/potentials/eam.cpp` — `cpu::EAM`
- `include/atomiverse/cpu/potentials/sw.h` + `src/…` — `cpu::StillingerWeber`（三体）
- `include/atomiverse/cpu/potentials/eam_data/eam_tables.h`（生成）
- `tools/gen_eam_tables.cpp`
- `regressions/md/eam_cu/`、`regressions/md/sw_si/`
- CUDA 镜像 **不在本 spec**（正交性规则要求最终补上，但探针先 CPU）

## Tasks

- [ ] `full_list_requested` 开关 + `NeighborList` 接上 `refresh_full`（含半表/全表一致性单测）
- [ ] halo cutoff 在有全表消费者时放宽到 `2 × r_cut`
- [ ] `EAMPotential(Element, name)` + 静态表 + 生成器
- [ ] `cpu::EAM`（两遍，单 launch 内；`ComponentBase + NeighborConsumer`）
- [ ] EAM 解析力 vs 有限差分 oracle（这是防止"对错泛函求导"重演的硬门）
- [ ] `cpu::StillingerWeber`（三体，验证全表足以支撑三元组）
- [ ] `CommField` 首个生产用户（`F'(ρ)` forward field）+ 多 rank 一致性
- [ ] regressions 用例 ×2
- [ ] 把撞出来的架构结论写回 `.claude/NOTES.md`

## Testing

- 单测（tests/，ms 级）：全表 vs 半表覆盖一致；EAM ρ 累加；SW 三元组枚举；表查找抛异常。
- **FD 力 oracle 是强制的**（见 CLAUDE.md 的 GPW 教训：能量对、力错，潜伏 17 个 commit）。
- regressions（库外，对文献）：见验收。

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

- Tersoff（SW 已足够证明三体通路；Tersoff 的键序求导量大，后续 spec）
- MEAM / ReaxFF / 反应势
- CUDA 镜像
- 每原子应力

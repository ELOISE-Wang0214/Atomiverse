---
title: Path driver + CI-NEB — stress-test the Driver taxonomy with N coupled replicas
slug: path-driver-neb
status: approved
science_required: true
created: 2026-07-12
acceptance: path-driver-neb.acceptance.md
depends_on: regressions-harness
---

# Path 驱动 + CI-NEB — 用 N 个耦合副本压测 Driver 分类学

## Summary

广度优先探针 #3，**架构价值最高的一个**。CLAUDE.md 自己预告了第六个 category
（"Path for NEB — terminate when every image of a coupled chain has converged"），
但从没人验证过这个预告是否成立。侦察的结论是：**它只对了一半。**

## 侦察结论 —— 最深的一条：分类学少了一个轴

CLAUDE.md 规定"category 由**终止语义**决定，绝不由算法或方法决定"。但
**"N 个耦合副本"是与终止语义正交的另一个轴**：

| | 全收敛才终止 | 步数预算终止 |
|---|---|---|
| **1 个 Frame** | `OptimizeDriver` | `DynamicsDriver` / `SamplingDriver` |
| **N 个 Frame** | **Path（NEB）** | **REMD / 并行回火 / PIMD / 伞形采样 —— 无家可归** |

照现有一轴分类，REMD/PIMD 将来会逼出**第七个 category**，而它与 Path 有 ~90% 的
多副本机器重复。

**必须在写 `path_driver.h` 之前做的决策（本 spec 的头号交付物）：**

- **(甲) 承认分类学是二维的**：把多副本机器抽成一个**可复用的 `cpu::ReplicaSet`**
  （新的基础设施簇，不是 Driver），由 Path 与将来的 Ensemble/REMD 共用。
- **(乙) 维持一轴**：`PathDriver` 私有持有 `vector<Frame>`，将来 REMD 复制一遍。

**本 spec 选 (甲)** —— 因为广度优先的全部意义就是提前暴露这种重复；等 PIMD 来的时候
再重构，代价是两倍。`ReplicaSet` 承载：副本存储、每副本一个 Backend/Context、
副本间耦合的挂载点、集体终止的归并。`PathDriver` 只提供"全收敛才终止"这一条语义。

## 用户指令（2026-07-12）：副本由**展开策略**生成，不由用户提供

> N 个 Frame 的 Driver 是不是还有 replicate 算法？最好的解决方案是输入一个 Frame，
> 然后有一个抽象算法给出如何将一个 Frame 展开成多个 Frame。因为即便是 REMD/PIMD
> 也不是要求用户提供多个不同的 Frame，而是通过某种策略展开。

**采纳，并且它修正了本 spec 原来的设计。** 原设计让 `ReplicaSet` 直接持有
`vector<Frame>`，隐含"用户交出 N 个 Frame"——既难用，又和门面
`Driver::run(frame, ...)`（**单个** Frame 入参）打架。

正确抽象是一个 **`ReplicaExpansion` 策略**：**一个 Frame 进 → N 个 Frame 出**。

| 方法 | 展开策略 | 每副本参数 |
|---|---|---|
| NEB | 端点间插值（线性 / IDPP） | — |
| PIMD | 恒等复制 × P 个 bead | bead 序号（ring-polymer 弹簧） |
| REMD / 并行回火 | 恒等复制 × N | **温度阶梯** |
| 伞形采样 | 恒等复制 × N | **bias 中心（沿 CV）** |

三个直接收益：

1. **门面签名一个字都不用改**：`PathDriver::run(frame, ...)` 仍然是 **Frame 进**；
   展开在内部发生。结束时把结果**折叠回那一个 Frame**（NEB → 攀爬 image = 鞍点构型），
   于是 **Frame 进、Frame 出**的货币不变式在新 category 上继续成立。
2. **per-replica 的参数（温度/bias 中心/bead 序号）天然归属于展开策略**，而不是散落在
   driver 字段里 —— 这正是把 category 焊死到 NEB 上的那种东西。
3. NEB 的"两个端点"不是"用户给 2 个 Frame"的特例，而是**展开策略自己的配置**：
   `LinearInterpolation{final_frame, n_images}`（或 IDPP）。初态是 `run()` 的入参，
   终态属于策略。

⇒ `ReplicaSet` 由 `(Frame, ReplicaExpansion&)` 构造，**不**由 `vector<Frame>` 构造。
`ReplicaExpansion` 是一个抽象基类，`expand(const Frame&) -> N frames + per-replica params`，
外加一个 `collapse(const ReplicaSet&) -> Frame`（把结果折叠回单 Frame）。

## 其它已核实的事实

**零多副本代码存在。** `cpu::Backend` 只持有一个 `Frame*`（`backend.h:175`）；
`Driver<>` 根本不持有 frame（frame 是 `run()` 的参数）；`Context` 里**没有 frame**
（`context.h:118-142`），组件只能通过 `ctx.frame` 看到**当前那一个**。

**`cpu::Frame` 是可拷贝的**（`frame.h:37-57`，全 xtensor SoA + std 容器，无 deleted 拷贝），
且已经在按值返回（`materialize()`）。所以 `std::vector<cpu::Frame>` **今天就能编译**。
代价是每个 image 深拷贝一份拓扑（浪费但正确）。

**不要复用 async 快照环**（`cpu/async/frame_snapshot.h:41-55`）：它是**有损**的 POD 子集
（无 mass/charge/拓扑/scalars），且是只写的发布通道，槽位会被回收。

**添加第六个 category 的完整改动面**（已逐文件追踪 Sampling 的加入过程）：
`core/binding.h:14`（`PhaseKind::Path` + `:19-33` 的 switch）、`facade/stage.h`（新 phase 命名空间）、
`cpu/driver/cpu_driver.h:114-119`（**必须**特化 `CategoryTraits<path::Phase>`，否则继承来的
categorical setter 无法实例化）、新建 `cpu/driver/path_driver.h`、
**`src/cpu/driver/driver.cpp` 的显式模板实例化**（`template class Driver<path::Phase, cpu::device_tag>;`）、
`CMakeLists.txt` 的 `ATOMIVERSE_SRCS` 显式列表（cpu/ 不是 glob）、
`tests/cpu/test_driver_binding_kind.cpp`（注释里写着 "five-driver"，要变六）。

**LBFGS 焊死在 Frame 上（真实技术债）：**
- 维度从 frame 推导：`lbfgs.cpp:100-102` `dim = owned_count(ctx)*3`，无注入维度的缝。
- 两环递推**内联在 `LBFGS::launch` 里**（`lbfgs.cpp:359-451`），不是自由函数 →
  band 级优化器只能**复制**这段算法，除非先把它抽成
  `lbfgs_two_loop(span x, span g, History&)`。
- **收敛写是无条件的全局写**：`lbfgs.cpp:201-216` 直接 `ctx.runtime.outer_status.converged = true`。
  共享一个 Context 的话，**image 0 收敛就会 break 掉整条链的循环**。这是最硬的一处。
- 可以干净复用的：`cpu/minimizer/line_search.h` 的 `armijo_backtrack` / strong-Wolfe /
  `zoom` **本来就是自由函数模板**（over 一个 `Phi` 仿函数）→ band 级 `phi(alpha)` 是 drop-in。

**MPI**：全树没有 replica communicator，DD 建在 `MPI_COMM_WORLD` 上。LAMMPS 式
"每副本一个 MPI 分区"会**死锁**。v1 只做：**band 在每个 rank 上复制，每个 image 内部
做空间 DD**。

**每副本一个 Backend ⇒ 每副本一个 Context ⇒ 每副本一个 `outer_status`，集体终止免费得到**
（`PathDriver` 把它们 AND 起来）。这也是绕开上面"LBFGS 无条件写收敛"的最省力办法。

## Design

- 新增基础设施簇 `cpu/replica/`：
  - `ReplicaExpansion`（抽象策略：`expand(const Frame&)` → N 个 Frame + 每副本参数；
    `collapse(...)` → 单个 Frame）。首个实现：`LinearInterpolation{final, n_images}`。
    `IdentityReplication{n}`（PIMD/REMD 将来直接用）作为第二个实现，**证明策略抽象
    不是为 NEB 量身定做的**。
  - `ReplicaSet`（由 `(Frame, ReplicaExpansion&)` 构造；每副本一个 Backend/Context；
    集体状态归并）。
- 新增 category 基类 `cpu/driver/path_driver.h`（`CPUDriver<ATV::path::Phase>`，
  终止 = 所有 image 的 `outer_status.converged` 全为真）。
- 新增域目录 **`cpu/path/`（不是 `cpu/neb/`）** —— 域目录按域命名，绝不按算法命名
  （`minimizer/` 装 LBFGS，`mc/` 装 TranslateMove）。叫 `neb/` 会把 REMD/PIMD 永远挡在门外。
- **相位枚举必须是通用的**：`force → couple → search_direction → line_search →
  convergence_check → end_of_step`。**不许出现 `tangent`/`spring`/`climb`** ——
  那样就把 category 焊死在 NEB 上了。`couple` 槽同样服务 PIMD 的 ring-polymer 弹簧
  （数学上就是同一个 k±1 谐振子）和 REMD 的交换算符。
- 切线 / 弹簧 / 攀爬像**组件**一样挂在 `couple` 槽（`role()==Evaluator`,
  `Capability::Force`），而不是写死在 driver 里 —— 遵守"算法是组件，不是 driver 子类"。

## Files

- `core/binding.h`、`facade/stage.h`、`cpu/driver/cpu_driver.h`（CategoryTraits 特化）
- `src/cpu/driver/driver.cpp`（显式实例化）
- `include/atomiverse/cpu/replica/replica_expansion.h`（策略基类 + `LinearInterpolation`
  + `IdentityReplication`）+ `replica_set.h` + `src/…`（新基础设施簇）
- `include/atomiverse/cpu/driver/path_driver.h`（新 category 基类）
- `include/atomiverse/cpu/path/{neb.h, tangent.h, spring.h, climbing_image.h}` + `src/…`
- `src/cpu/minimizer/lbfgs.cpp` → 抽出两环递推（**去焊接**）。
  **注意这与"只运行一次的代码不要抽成函数"不冲突 —— 恰恰相反，它是那条规则的正例：
  抽取发生在第二个调用者出现的时候**（单 Frame 的 `Minimize` + band 级优化器）。
  在此之前它内联在 `LBFGS::launch` 里是对的；现在有了第二个调用者，才抽。
  形式上优先做成可复用的类（`LbfgsHistory` 持有曲率历史 + 一个作用于 span 的方法），
  而不是一个裸自由函数 —— OOP 优先。
- `CMakeLists.txt`、`tests/cpu/test_driver_binding_kind.cpp`（five → six）
- `CLAUDE.md` 分类学表 + `docs/roadmap/optimization.md`
- `regressions/path/neb_muller_brown/`

## Tasks

- [ ] **裁决并记录分类学的二维性**（选 (甲)：抽出 `ReplicaSet`）
- [ ] `PhaseKind::Path` + `path::Phase`（**通用**相位名）+ `CategoryTraits` 特化 + 显式实例化
- [ ] **`ReplicaExpansion` 策略**：`expand(Frame) -> N`，`collapse(...) -> Frame`；
      实现 `LinearInterpolation`（NEB）**与** `IdentityReplication`（证明可复用）
- [ ] `ReplicaSet`：由 `(Frame, ReplicaExpansion&)` 构造 / N 个 Backend / 集体
      `outer_status` 归并（AND 收敛，OR 停止）
- [ ] `PathDriver`（category 基类）
- [ ] 把 LBFGS 的两环递推抽成自由函数（先去焊接，再复用）
- [ ] NEB 切线（Henkelman 改进切线）+ 弹簧 + **攀爬 image**，全部作为 `couple` 槽组件
- [ ] regressions 用例（Müller-Brown 解析势的鞍点）
- [ ] 架构结论写回 `.claude/NOTES.md` + 更新 CLAUDE.md 分类学表

## Testing

- 单测：`ReplicaSet` 的集体收敛归并（一个未收敛 ⇒ 整体未收敛）；跨 workflow 的
  PhaseKind 误绑定仍然抛；切线公式对着手算。
- regressions：Müller-Brown 势的鞍点位置与能垒对着文献值。

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

- REMD / PIMD / 伞形采样（`ReplicaSet` 为它们而抽，但本 spec 不实现它们）
- Dimer / string method
- 每副本一个 MPI 分区（v1：band 全 rank 复制）
- CUDA

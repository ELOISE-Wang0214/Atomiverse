---
title: Delete the Driver category — a Driver is a staged loop, termination is a Handler
slug: driver-handler
status: approved
science_required: false
created: 2026-07-12
acceptance: driver-handler.acceptance.md
supersedes: the CLAUDE.md "Driver taxonomy" section in its entirety
---

# 删掉 category —— Driver 是一个分阶段的循环，终止是一个 Handler

## Summary

三个探针连续撞出同一堵墙，墙的名字叫 **category**：

- probe 2：LBFGS **无条件**往 `ctx.outer_status.converged` 里写 —— 收敛判据（一个**策略**）被埋进了优化器（一个**算法**）。
- probe 3：因此第一个 image 一收敛就会停掉整条 NEB 链。为了绕开它，`ReplicaSet` 只好自己再维护一份 per-replica 的收敛表。
- probe 3 又发现：要表达"N 个耦合副本"，一维的 category 不够，于是被迫承认第二根轴。

**这三件事是同一个病：终止语义被焊死在类继承体系里，于是它既不能复用，也不能组合，还每来一个新方法就长一个新 category。**

本 spec 删掉 category 这个概念。

## 用户裁决（2026-07-12，binding）

1. **category 直接删掉**，改成**可复用的 handler**。
2. `MDDriver` 这类**方法相关的 pipeline**，改成**默认构建好的一套 driver + handler**（preset），不是继承体系的一层。
3. **stage 必须用 enum** —— 不许 string，不许 int。

## Driver 是什么

剥到底，只有三样：

1. **Stages** —— 带类型的相位序列（enum），Plan 把节点绑在上面。
2. **Iteration** —— 重复走 stages。
3. **Termination** —— 什么时候停。

**(3) 被焊进类继承里，这就是 category。** 六个类，彼此只差一个 `if`。把 (3) 拿走，category 不再是类。

## 新的三分

| | 是什么 | 数量 |
|---|---|---|
| **Runtime**（今 Backend+Context） | frame 绑定、节点存储、timing、MPI | **每次 run 一个**，整棵 driver 树共享 |
| **Driver** | **一个循环节点**：stages + plan + handlers | 任意多，**可嵌套** |
| **Handler** | 每轮迭代后被问一次 | 任意多，可组合，**用户可定义** |

两条关键结构决定：

- **Driver 不拥有 Runtime**，它被交给一个 Runtime 去执行 → 嵌套时 Backend 不翻倍、DD 不重建、MPI 不分裂。
- **Driver 本身是一个可执行节点**（build/launch/destroy，和组件同样的生命周期）→ Plan 的 slot 里既能放组件，也能放**子 Driver**。"driver 能不能调 driver"这个问题于是自己消失：Driver 是节点树里的节点。

## Handler

```cpp
enum class Verdict { Continue, Converged, Stop };

/// 每轮迭代之后被问一次。**不是模板** —— 这正是它能跨 workflow 复用的原因：
/// 它看不见 stage，只看见迭代状态。
class LoopHandler {
  public:
    virtual ~LoopHandler() = default;
    virtual Verdict after_iteration(const LoopState&) = 0;
};
```

`LoopState{ int iteration; const Frame&; const Context&; const LoopMetrics& }` —— **没有 phase**。

今天的六个 category，退化成六个 handler：

| 原 category | 现在只是 |
|---|---|
| Evaluate | `SinglePass{}` |
| Dynamics | `StepBudget{n}` |
| Optimize | `ForceConverged{ftol, rms_ftol, ...}` |
| Response | `FixedBudget{6N}` |
| Sampling | `TrialBudget{n}` |
| Path | `AllReplicasConverged{}` |

而用户第一次能自己写：`StopOnNaN`、`StopWhenTemperatureExceeds{2000}`、`StopOnWallClock{3600}`、`StopOnSentinelFile`（= LAMMPS `fix halt`）。**今天这些一个都写不出来**，因为"停不停"是 driver 类的私有 `if`。

**组合规则（定死）：Stop > Converged > Continue** —— 按严重性取第一个非 Continue。

**MPI 契约：handler 的裁决必须是集合的**（allreduce）。各 rank 停在不同迭代上 = 死锁。今天 LBFGS 是手工做这件事的；新设计里这是 `LoopHandler` 基类的契约，由 driver 统一 allreduce，handler 作者不用想。

### 这一步单独就修掉两个已知 bug

LBFGS 不再写收敛，它只**发布度量**（`max_force` / `rms_force` / `ΔE` → `LoopMetrics`），由 handler **裁决**。于是：

- 换收敛判据不必改优化器；
- **status 属于循环节点，不属于 Context** —— 子循环看不见、也够不着父亲的 status。probe 2 和 probe 3 那个"内层一收敛就停掉外层"的 bug，从"需要记住的纪律"变成**结构上写不出来的东西**。

## Outcome policy —— MC 的接受/拒绝不是终止

它是**迭代结果策略**：每轮之后决定 commit 还是 rollback。与 handler 正交的第二个策略位：

```cpp
enum class Outcome { Commit, Rollback };
class OutcomePolicy {
    virtual Outcome after_stages(LoopState&) = 0;   // Metropolis 是它的一个实现
};
```

分开之后 `SamplingDriver` 这个类没有存在的必要。

## stage 必须是 enum（用户裁决 3）

- **公开 API 里不许出现 int stage。** 要修的洞：`Backend::run_phase(int)` → `run_stage(PhaseEnum)`（模板化，或由 `Driver<PhaseEnum>` 在**唯一一处** `slot_of()` 里做 enum→槽位换算）。当前 `path_driver.h` 里的 `run_phase(static_cast<int>(p))` 是这个洞的现行受害者。
- **`DefaultBinding` 里的 `int phase` 保留，且只能经 `make_binding(PhaseEnum, order)` 构造。** 它是**擦除本身**：`ComponentBase` 不是模板，`default_bindings()` 必须返回同一个类型。`PhaseKind` 由 enum 类型经 ADL 推出，跨 workflow 误绑在注册时抛。**没有人在写代码时碰得到那个 int。**
- **嵌套不需要 int**：子 Driver 作为**节点**被持有，它自己的 stage 仍是它自己的 enum。**类型擦除发生在节点边界，不在 stage 边界。**

## Preset（用户裁决 2）

`MD` / `Minimize` / `MC` / `Evaluate` / `Hessian` / `NEB` 活下来，但**不再是继承体系的一层**：它们是**在构造函数里配好默认 handler + 默认基础设施注入**的 `Driver<Phase>`。

```cpp
class MD : public Driver<md::Phase> {
  public:
    MD() { set_handler(StepBudget{});  /* + 邻居表/halo 的自动注入规则 */ }
};
```

**是类，有构造函数，不是 `make_md()` 工厂函数**（CLAUDE.md 硬规则）。

继承深度：`Driver<P,Dev> → CPUDriver<P> → {6 个 category 基类} → {5 个具体 driver}` **四层**，塌成 `Driver<P> → {preset}` **两层**。

## 分类学之后剩下什么

Category **不再是类**。剩下的是 **stage 词汇表**（`md::Phase`、`minimizer::Phase`、`path::Phase`…）—— 这些是真的，它们命名了"一轮迭代里发生什么"。preset 属于**便利**，不属于架构。

probe 3 那个"二维分类学"的发现一并折叠：**副本数根本不是 category 的一根轴，它是一层嵌套循环。** REMD/PIMD 于是**不需要任何新 category** —— 它们是"副本循环 + 预算 handler + 不同的 couple 组件"。

## 诚实的未解项

**组件的每实例状态。** 同一个子 Driver body 在 N 个副本上跑，里面的组件就要 N 份独立状态（LBFGS 的曲率历史）—— 这**又回到 probe 3 撞到的"组件不能克隆"**。嵌套让这个问题**更尖锐**，没让它消失。三条出路都还没验：(i) 每副本一个子 Driver 实例（需要 `clone()`）；(ii) 组件状态按循环实例做键（脏）；(iii) 保持 band 级组件（NEB 今天的做法，能跑但没真正用上嵌套）。**本 spec 不裁决它**，只保证不把它变得更难。

## 落地顺序（每步独立可验，每步跑全套）

- **① Runtime / Driver 分离** —— Driver 不再拥有 Backend；`run(runtime, frame)`。纯机械，无行为变化。
- **② Handler + per-loop status** —— **最小、最独立、立刻验证设计是对的**：判据从 LBFGS 抽出来，probe 2/3 那个 bug 结构性消失。`Verdict` / `LoopHandler` / `LoopMetrics` / 组合规则 / MPI allreduce 契约。
- **③ `run_stage(PhaseEnum)`** —— 堵上 int 洞（裁决 3）。
- **④ Driver 成为可执行节点，允许嵌套** —— Plan slot 收 component | sub-driver。
- **⑤ 六个 category 基类删除，preset 化**（裁决 1 + 2）。CUDA 的 3 个镜像同步。

## Testing

- 单测：handler 组合规则（Stop > Converged > Continue）；per-loop status 隔离（子循环写 converged **不会**停掉父循环 —— 这条直接钉死 probe 2/3 的 bug）；跨 workflow 误绑仍抛；preset 的默认 handler 是对的那个。
- regressions：**现有 5 个用例全部原样通过，且不改一行**（它们是库外消费者，用的是公开 API —— 如果 preset 做对了，它们不该知道底下换了架构）。这是本 spec 最强的验收信号。

## Out of scope

- 组件 `clone()` / 每实例状态（见"诚实的未解项"）。
- REMD / PIMD（新架构为它们而生，但本 spec 不实现）。
- torch 后端。

---
title: Frame 的组织方式 — 后端分类 vs 后端×领域分类（决策 spec）
slug: frame-organization
status: draft
science_required: false
created: 2026-07-12
acceptance: frame-organization.acceptance.md
depends_on: manybody-potentials, md-constraints, path-driver-neb, mc-gcmc-npt
---

# Frame 的组织方式（决策 spec）

## Summary

用户提问（2026-07-12）：

> Frame 需要更好的解决方案，是按照后端分类，还是后端 × 领域分类？

**这是一个决策 spec，不是实现 spec。它刻意排在四个探针之后** —— 因为四个探针正是
回答这个问题所需要的证据。现在就重构 Frame 是投机；跑完探针再重构是有据可依。

## 现状

Frame 按**后端**分类，一个后端一个类型，语义字段集共享、表示各自拥有：

| Frame | 表示 | 拓扑布局 |
|---|---|---|
| `molrs.Frame` | 扁平 Block 列 —— **边界上的交换 schema** | 扁平 bonds/angles/dihedrals block |
| `cpu::Frame` | 主机端 xtensor SoA | 每原子 LAMMPS 布局（`num_bond`/`bond_atom`，全局 tag） |
| `cuda::Frame` | 裸设备指针 + 主机端缓存计数 | **双份**：单卡扁平 list / 多卡每原子 SoA |
| `torch::Frame` | `at::Tensor` block（可自动微分） | 扁平 sub-struct |

转换点各自拥有恰好一个边界：`materialize`/`write_back`、`to_cuda`/`to_cpu`、`torch::transfer`。

## 三个候选

### (甲) 维持现状：只按后端分类
- 优点：转换边界恰好 3 个，往返契约已被金标测试钉死（`test_materialize.cpp`）。
- 缺点：Frame 有变成 god-struct 的压力 —— 每个新领域都想往里加字段。

### (乙) 后端 × 领域分类（`cpu::md::Frame`、`cpu::scf::Frame` …）
- 优点：每个领域只看见自己要的字段。
- **致命缺点：转换边界从 3 个变成 3×D 个。** 每一个都要有自己的 `materialize`/
  `write_back`/`to_cuda` 与往返契约。而 AIMD 这类跨领域工作流（MD 驱动 + SCF 求值器）
  会**同时**需要两种 Frame —— 于是要么再加领域间转换，要么退化回一个大 Frame。
  **建议否决。**

### (丙) 一个后端一个 Frame，但可选载荷是**可组合的 block**（**倾向选项**）
核心字段（x/v/f/mass/type/charge/tag/image/box/scalars）恒在；领域按需挂载 block
（拓扑 block、per-atom 派生量 block、影子数组 block …）。

**这不是新发明 —— 树里已经在手搓它了：**
- `cpu::Atom` 的 `has_bonded_topology_` / `ensure_bonded_storage()`（`frame/atom.h:258,267`）
  就是一个**手写的、简陋的 optional block**：纯 LJ 体系不为键角二面角数组付内存。
- **molrs 的 `Frame` 本来就是 Block 列的组合** —— block 模型在边界的另一侧已经是既成事实，
  只是 C++ 侧没有把它表达出来。

⇒ 方向是 **composition，不是 multiplication**：转换边界仍是 3 个，领域各取所需。

## 四个探针会提供的证据（这就是本 spec 排在最后的原因）

| 探针 | 它会告诉我们 Frame 缺什么 |
|---|---|
| `manybody-potentials` | per-atom 派生量（EAM 的 ρᵢ、F′(ρ)）该住在哪：组件自有数组？还是一个 Frame block？（现在被迫是组件自有 + `CommField`） |
| `md-constraints` | 约束需要"未约束坐标"的**影子数组** `x_unc`/`v_unc`。这是 Frame block，还是组件私有？ |
| `mc-gcmc-npt` | **变长 N 已经在 `comm/` 里解决了**，且 `Transaction` 将被改写为 **Frame 补丁** —— 这会直接压测"Frame 能不能被切片/拼接" |
| `path-driver-neb` | `ReplicaExpansion` 要求 Frame 能被**展开与折叠**（一个 Frame → N 个 → 一个）。这是对 Frame 可组合性最强的一次压测 |

**四个探针合起来，会把"Frame 到底需要什么"从猜测变成观测。** 特别是 GCMC 的
Frame-补丁 与 NEB 的 expand/collapse —— 两者都要求 Frame 支持**部分 Frame / 子 Frame**
的语义，而这恰恰是 block 模型天然给出、而当前固定 struct 给不出的能力。

## Tasks

- [ ] 等四个探针 `code-complete`
- [ ] 汇总四份 `.claude/NOTES.md` 的架构结论，列出 Frame 被真实需求撑出的每一个洞
- [ ] 就 (甲)/(乙)/(丙) 给出带证据的裁决，交用户签字
- [ ] 若选 (丙)：设计 block 模型（与 molrs 的 Block 列对齐），保持转换边界仍为 3 个，
      并保留 `test_materialize.cpp` 的往返金标契约不破

## Out of scope

- 本 spec **不写任何实现代码**，直到裁决被签字。

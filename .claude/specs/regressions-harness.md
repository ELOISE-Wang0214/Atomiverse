---
title: regressions/ — out-of-library textbook end-to-end regression harness
slug: regressions-harness
status: approved
science_required: false
created: 2026-07-12
acceptance: regressions-harness.acceptance.md
---

# `regressions/` — 库外教科书级端到端回归

## Summary

用户指令（2026-07-12）：

> 每个功能交付的时候要在库外尝试运行一个端到端的、教科书级别的例子确保和参考文献
> 一致，作为回归测试的一部分，这个例子要尽可能小，不要和单元测试放在一起，可以专门
> 放一个 regressions/ 文件夹。

本 spec 建立这个 harness。它是后面四个广度探针 spec（`manybody-potentials`、
`md-constraints`、`path-driver-neb`、`mc-gcmc-npt`）的**共同前置**：每个探针交付时
必须在 `regressions/` 里留下一个通过的教科书用例。

## 为什么不是 tests/ 也不是 benchmarks/

树里已经有两个地方，都不合适，所以确实需要第三个：

| 目录 | 现有语义 | 为什么不合适 |
|---|---|---|
| `tests/` | 纯单元测试，ms 级，一个函数/方法。CLAUDE.md 硬规则：**不允许 e2e**、不允许计时、不允许与外部引擎比对，且"注册了就会跑" | 用户明确要求"不要和单元测试放在一起" |
| `benchmarks/<domain>/<case>/{validation,performance}` | 与外部**引擎**（LAMMPS/PySCF/CP2K）比对 + 性能测量，验收时手工跑，不进 ctest | 语义是"对外部引擎的 parity"与"测速"，且不是每次回归都跑 |
| **`regressions/`（本 spec）** | **对着教科书/文献的解析解或公认数值**，跑通一个最小 e2e，作为**回归**每次都能跑 | — |

判据差别一句话：`benchmarks/` 问"我们和 LAMMPS 一样快/一样准吗"；`regressions/`
问"我们还和**教科书**一样对吗"。后者不需要装任何第三方引擎，因此可以常态化跑。

## Design

### 布局

```
regressions/
  CMakeLists.txt              定义 atv_add_regression()
  README.md                   规则 + 用例清单表
  <domain>/                   md | mc | qc | path …（与 benchmarks 同名法）
    <case>/
      reg_<case>.cpp          单文件、main() 返回非零即失败
      REFERENCE.md            引用出处（书/论文 + 页码/式号 + 期望值 + 容差理由）
```

### 硬规则（写进 regressions/README.md，也是验收判据）

1. **库外消费者视角**：只允许 include `atomiverse/` 公开头文件，只 link `atomiverse`
   目标（照抄 `examples/CMakeLists.txt` 的 `atv_add_cpp_example` 模式）。**禁止**
   include 任何 `internal/` 头、禁止 friend、禁止访问私有实现。一个回归用例跑不通，
   要么是物理错了，要么是**公开 API 不足以完成这件事** —— 后者本身就是要暴露的缺陷。
2. **对着文献，不对着我们自己**。期望值来自教科书/论文的解析解或公认数值，写在
   `REFERENCE.md` 里，含**完整出处与容差的物理理由**。禁止"黄金值"式地把我们当前
   的输出固化成期望值。
3. **尽可能小**：单个用例 CPU 串行 **≤ 60 秒**（远严于验收 10 分钟线）。这是回归，
   不是生产。
4. **不进默认 ctest**：遵守 CLAUDE.md「tests/ 注册了就会跑、且不含 e2e」。回归用
   独立开关 `ATV_BUILD_REGRESSIONS=ON` 构建，独立目标 `regressions` 一次跑完。
5. **非零退出即失败**，并把 `got / want / tol` 打印出来。
6. **每个 case 一个文件，写成直白的 `main()`。** 用户指令：**只运行一次的代码不要抽出来
   做成单独的函数**。回归用例是一次性脚本 —— setup / run / check 就地写在 `main()` 里，
   不许为了"好看"切成 `setup_frame()` / `run_md()` / `check_energy()` 一串一次性 helper。
   一个读者应当从上往下读完就懂，不必在文件里跳来跳去。

### 三层测试的分工（同步写进 CLAUDE.md）

| 层 | 跑多久 | 测什么 | 硬禁 |
|---|---|---|---|
| `tests/` | **整套秒级** —— 快速迭代循环 | **单个函数 / 方法** | **严禁 e2e** |
| `regressions/` | 单例 ≤ 60 s | **库外 e2e，对着文献** | 不进默认 ctest |
| `benchmarks/` | 验收时 | 对**外部引擎** parity + 测速 | 不进 ctest |

**每交付一个功能，必须在 `regressions/` 留下恰好一个教科书级用例** —— 而不是把 e2e
塞进单元测试里。

### CMake

- 新增 `option(ATV_BUILD_REGRESSIONS "Build end-to-end textbook regression cases" OFF)`。
- `regressions/CMakeLists.txt` 提供 `atv_add_regression(<domain> <case>)`：建 `reg_<case>`
  可执行、link `atomiverse`、挂到聚合目标 `regressions` 上（`add_custom_target(regressions)`
  依赖并**运行**全部用例）。
- **不** `add_test()` 到默认 ctest（见规则 4）。

## Files

- `regressions/CMakeLists.txt`（新）
- `regressions/README.md`（新）
- `CMakeLists.txt`：加 `ATV_BUILD_REGRESSIONS` 选项 + `add_subdirectory(regressions)`
- `regressions/md/lj_argon/`：**首个种子用例**，见下
- `CLAUDE.md`：Repository layout 增 `regressions/` 一行；tests/benchmarks 的二分法改三分法
- `docs/developer/index.md`：同上

### 种子用例：`regressions/md/lj_argon/`

用已实现且已 ⚡ 的 LJ + NVE 做种子，证明 harness 本身可用（而不是同时调试 harness 和新物理）：
液氩 LJ NVE **能量守恒** —— 固定 dt 下总能量漂移 < 1e-4 ε/atom/ps。参考：
Rahman (1964) *Phys. Rev.* **136**, A405 的 LJ 氩参数（ε/k_B = 120 K, σ = 3.4 Å）；
守恒判据是 velocity-Verlet 的辛性质（Frenkel & Smit §4.3）。

## Tasks

- [ ] `ATV_BUILD_REGRESSIONS` 选项 + `add_subdirectory(regressions)`
- [ ] `regressions/CMakeLists.txt` + `atv_add_regression()` + 聚合目标 `regressions`
- [ ] `regressions/README.md`（五条硬规则 + 用例清单表）
- [ ] 种子用例 `regressions/md/lj_argon/reg_lj_argon.cpp` + `REFERENCE.md`
- [ ] CLAUDE.md / docs 的 tests-vs-benchmarks 二分法 → 三分法

## Testing

harness 自身无单元测试（它是构建设施）。验收 = 种子用例在干净 build 里绿。

## Out of scope

- CI 接线（谁在什么时候跑 `regressions` 目标）——留给后续。
- CUDA 回归用例。先 CPU。

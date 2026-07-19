---
slug: python-on-demand-install-01-cmake-modules
title: CMake product-module gates
status: code-complete
grilled: true
created: 2026-07-19
chain: python-on-demand-install
phase: 01-cmake-modules
---

# CMake product-module gates (01)

## Summary

在根 CMake 中为 Atomiverse 产品域建立可独立开关的模块选项与无环依赖 DAG：关闭的模块不再进入 `atomiverse` 的编译单元，并对外暴露与既有 `ATV_HAS_PW` 同形的 `ATV_HAS_*` 宏；nanobind 仅链接已开启模块的 `bind_*.cpp`。默认 monorepo / CI 配置仍打开全部经典域表面，保证现有一条 configure 行与 `tools/ci_cpu_unit.sh` 行为不回退。

## Domain basis

本子-spec 为 **构建/打包模块化**，无物理方程、无 DOI/arXiv 校验目标。模块边界对齐仓库后端目录（`md/`、`mc/`、`minimizer/`、`path/`、`vib/`、`semiempirical/`、`potentials/`、`neighbor/` 等），不改变力场、积分器或电子结构算法。

## Design

### Always-on core（无独立 `ATV_BUILD_*` 产品门）

始终编入 `atomiverse`：

| 集群 | 路径（示意） | 说明 |
|---|---|---|
| Core | `src/core/*` | log、molrs contract |
| Frame / Driver / Component | `src/cpu/frame/*`, `driver/*`, `component/*` | `cpu::Frame`、`CPUDriver`、`Evaluate` preset、LoopHandler / SubDriver |
| I/O + 通用 recorder/observer/async/comm/task | `io/*`, `recorder/*`, `observer/*`, `async/*`, `comm/*`, `task/*` | 跨 workflow 基础设施 |
| CUDA 设备骨架（`ATV_BUILD_CUDA`） | `cuda/device/*`, `cuda/frame/*`, `cuda/driver/*`, `cuda/async/*`, `cuda/comm/*` | 与现有 CUDA 门一致 |
| Toy SCF | `src/cpu/scf/solve.cpp` | 现有始终编入的 toy 参考，不并入 `ATV_BUILD_RHF` |
| Direct Ewald | `src/cpu/kspace/ewald.cpp` | 始终可编；PME 仍仅 `ATV_BUILD_KSPACE` |

`ATV_BUILD_EVALUATE`（默认 ON）仅用于发出 `ATV_HAS_EVALUATE` 与产品矩阵对称，不额外 gate 源文件。

### 新产品 / 基础设施门（默认 ON）

| Option | 默认 | 源目录 | Public def |
|---|---|---|---|
| `ATV_BUILD_MD` | ON | `src/cpu/md/*`；`src/cuda/md/*` | `ATV_HAS_MD` |
| `ATV_BUILD_MC` | ON | `src/cpu/mc/*` | `ATV_HAS_MC` |
| `ATV_BUILD_MINIMIZER` | ON | `src/cpu/minimizer/*`；`src/cuda/minimizer/*` | `ATV_HAS_MINIMIZER` |
| `ATV_BUILD_VIB` | ON | `src/cpu/vib/*` | `ATV_HAS_VIB` |
| `ATV_BUILD_PATH` | ON | `src/cpu/path/*`, `src/cpu/replica/*`, `src/cpu/bias/*` | `ATV_HAS_PATH` |
| `ATV_BUILD_SEMIEMPIRICAL` | ON | `src/cpu/semiempirical/*` | `ATV_HAS_SEMIEMPIRICAL` |
| `ATV_BUILD_POTENTIALS` | ON | `src/cpu/potentials/*`；`src/cuda/potentials/*` | `ATV_HAS_POTENTIALS` |
| `ATV_BUILD_NEIGHBOR` | ON | `src/cpu/neighbor/*`；`src/cuda/neighbor/*` | `ATV_HAS_NEIGHBOR` |
| `ATV_BUILD_EVALUATE` | ON | （无额外 `.cpp`） | `ATV_HAS_EVALUATE` |

既有门保持：`ATV_BUILD_{RHF,PW,GPW,KSPACE,CUDA,TORCH,PYTHON,…}`。

### Auto-enable DAG（单向、无环）

在 option 声明之后、组装源列表之前解析；子依赖 `CACHE BOOL FORCE` + `message(STATUS)`：

```
MD            → POTENTIALS, NEIGHBOR
MC            → POTENTIALS, NEIGHBOR
MINIMIZER     → POTENTIALS, NEIGHBOR
PATH          → MINIMIZER   # 再经 MINIMIZER 拉 POTENTIALS+NEIGHBOR
KSPACE        → POTENTIALS
GPW           → RHF ∧ PW    # 已有 FATAL，不 auto-enable 缺失端
TORCH         → CUDA
LAUNCHER      → MPI ∧ PYTHON
```

禁止环。无效组合 FATAL（保留 GPW/TORCH/LAUNCHER 等既有契约）。

### 源文件与 nanobind

1. 经典域按 `if(ATV_BUILD_*)` 追加 `target_sources`（对齐 RHF/PW 块）。
2. 开启门 → `target_compile_definitions(… PUBLIC ATV_HAS_*)`。
3. CUDA 与 CPU 同门正交：`if(ATV_BUILD_CUDA AND ATV_BUILD_MD)` 等。
4. 测试/回归注册按门 wrap。
5. `nanobind_add_module` 仅加入 ON 模块的 `bind_*.cpp`；`bind.cpp` 用 `#ifdef ATV_HAS_*` 调 init。

实现落点：`cmake/AtvProductModules.cmake`（新）+ 根 `CMakeLists.txt`。`python/pyproject.toml` **本子-spec 不修改**。

### Reuse decision

- **reuse** 既有 `option(ATV_BUILD_*)` / `ATV_HAS_*` / `FATAL_ERROR` 模式。
- **generalize** `ATOMIVERSE_SRCS` 始终编入的经典域 → always-on 核心 + 按门追加。
- **reuse** `nanobind_add_module(_atomiverse)` — 仅条件化源列表。
- **pattern** `bind_pw.cpp` 的 `#ifdef ATV_HAS_PW` — 推广到 `bind.cpp` 与其余 bind。
- **new — deferred** `pyproject.toml` cmake.define → 02/03。

## Files to create or modify

- `cmake/AtvProductModules.cmake` (new)
- `CMakeLists.txt`
- `python/src/bindings/bind.cpp`
- `tests/CMakeLists.txt`
- `tests/architecture/module_dag_auto_enable.cmake` (new)
- `tests/architecture/module_off_excludes_sources.cmake` (new)
- `regressions/CMakeLists.txt`
- `regressions/architecture/module_compile_defs/reg_module_compile_defs.cpp` (new)
- `regressions/architecture/module_compile_defs/REFERENCE.md` (new)

## Tasks

- [x] Write failing architecture cmake tests for product-module DAG and source gating (`tests/architecture/module_dag_auto_enable.cmake`, `tests/architecture/module_off_excludes_sources.cmake`; register in `tests/CMakeLists.txt`)
- [x] Write failing regression compile-def checks in `regressions/architecture/module_compile_defs/reg_module_compile_defs.cpp`
- [x] Implement `cmake/AtvProductModules.cmake` options, auto-enable DAG, FATAL messages, and STATUS summary
- [x] Generalize `ATOMIVERSE_SRCS` / `target_sources` in `CMakeLists.txt` into always-on core + per-module gated lists; emit public `ATV_HAS_*`
- [x] Gate `nanobind_add_module` sources and wrap `bind_*` init in `python/src/bindings/bind.cpp` with `#ifdef ATV_HAS_*`
- [x] Gate unit-test and regression registration in `tests/CMakeLists.txt` and `regressions/CMakeLists.txt` by `ATV_BUILD_*`
- [x] Add regression example `regressions/architecture/module_compile_defs` (public-header consumer of compile defs; ≤60 s)
- [x] Verify against default monorepo configure plus one modular OFF smoke (`-DATV_BUILD_MD=OFF`)
- [x] Run full check + test suite

## Testing strategy

- **Happy path**：默认经典门全 ON 时与现有 ctest 注册表一致；`ATV_HAS_*` 齐全。
- **DAG / OFF 源排除**：architecture cmake 脚本。
- **nanobind**：`MD=OFF` 则无 `bind_md.cpp`。
- **Regression**：`module_compile_defs` 编译期要求 full-classical `ATV_HAS_*`；不跑 SCF/MD e2e。

## Out of scope

- Python CLI / pure wheel / 源码拉取 / `--config`（02/03）
- 按模块裁剪安装头或拆多静态库
- 改变算法 / Frame / Driver 语义
- `tests/` 内 e2e SCF/MD
- 用户手册长文（仅 CMake 注释）

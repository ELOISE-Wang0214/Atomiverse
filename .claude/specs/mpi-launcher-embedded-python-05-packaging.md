---
title: MPI launcher packaging — wheel 内置原生 atomiverse、CMake install、运行时解释器发现
status: code-complete
created: 2026-07-15
---

# MPI launcher packaging — wheel 内置原生 atomiverse、CMake install、运行时解释器发现

## Summary

本篇是 mpi-launcher-embedded-python 五链的收官（05）：把 04 落地的原生 `atomiverse` 可执行程序（rank 0 内嵌 CPython 的 MPI 启动器）真正交付到用户手里。三条交付路径与一个核心设计面：（1）wheel 路径——scikit-build-core 构建的 wheel 可选地（默认 OFF、环境变量 `ATV_BUILD_LAUNCHER=ON` 显式开启）把原生 ELF 二进制装进 wheel 的 scripts 区，`pip install` 后它以真实原生二进制的身份出现在 venv/conda 的 `bin/` 即 PATH 上，明确不是 console-script 垫片；（2）纯 CMake 路径——`cmake --install` 为非 wheel 用户安装 `bin/atomiverse` 及最小的库/头文件布局；（3）运行时解释器发现——内嵌 CPython 不再依赖 04 假设的构建期 FindPython3 路径，而是在启动时解析用户当前环境（venv/conda/module-loaded Python）的 `python3` 并据此初始化，配 `ATV_PYTHON` 环境变量兜底覆盖，不支持的布局给出清晰致命错误、绝不静默回退到构建机路径；（4）一页 docs/user-guide 文档讲清 `python atv.py` 与 `mpirun -np N atomiverse atv.py` 两种运行方式及环境要求。验收标准是"可安装 + 可发现"，不是物理。

## Design

**不涉及物理**（无 Domain basis）；工程依据为 CPython PEP 587（`Py_InitializeFromConfig` / `PyConfig` 路径计算）与 scikit-build-core ≥ 0.10 的 wheel 布局机制。

### 1. Wheel 内的原生二进制（scikit-build-core scripts 机制，精确表述）

scikit-build-core（本仓已 pin `>=0.10`）在 CMake configure 时定义 `SKBUILD=2` 并注入安装目录变量族，其中 **`SKBUILD_SCRIPTS_DIR`** 映射到 wheel 的 `atomiverse-<ver>.data/scripts/`。放进该目录的文件在 `pip install` 时被安装到目标环境的 `bin/`（Windows 为 `Scripts/`）并置可执行位；pip 只对以 `#!python` 开头的文本脚本重写 shebang，**ELF 二进制原样拷贝**——这正是 patchelf 等 PyPI 包发行原生可执行文件的既有机制。因此 04 的 launcher target 只需：

```cmake
# app/atomiverse/CMakeLists.txt
include(GNUInstallDirs)
if(SKBUILD)
    install(TARGETS atomiverse_launcher RUNTIME DESTINATION "${SKBUILD_SCRIPTS_DIR}")
else()
    install(TARGETS atomiverse_launcher RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
endif()
```

（`atomiverse_launcher` 为 04 的可执行 CMake target 名：库 target 已占用 `atomiverse` 名，可执行 target 以 `OUTPUT_NAME atomiverse` 落盘。）

**决策：launcher-in-wheel 默认 OFF、opt-in。** `python/pyproject.toml` 的 `[tool.scikit-build.cmake.define]` 新增一行，沿用现有 env-override 模式：

```toml
ATV_BUILD_LAUNCHER = {env = "ATV_BUILD_LAUNCHER", default = "OFF"}
```

理由：launcher 依赖 `ATV_USE_MPI`（MPI wheel 天然绑定环境，不可移植）且直接链接宿主 `libpythonX.Y`（非 manylinux 可分发）；HPC 用户本就 `pip install .` 源码构建。CMake 侧加守卫：`ATV_BUILD_LAUNCHER=ON` 而 MPI 未找到 → `FATAL_ERROR`（不回退、点名开关，与 143–148 行 `ATV_USE_MPI` 守卫同风格）。`ATV_USE_MPI` 保持根 CMakeLists 的自动探测默认，不进 pyproject。

### 2. 纯 CMake install（现状为零，加最小规则）

根 CMakeLists 今天唯一的 install 是 nanobind 模块（759 行）。新增（`include(GNUInstallDirs)`，全部 `if(NOT SKBUILD)` 守卫——wheel 不携带头文件/静态库）：

- `install(TARGETS atomiverse ARCHIVE/LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})`
- `install(DIRECTORY include/atomiverse TYPE INCLUDE)`
- 可执行文件走上面 app/ 的规则（launcher 开启时）。

已知局限（文档化、不在本篇解决）：`atomiverse` 是静态库且传递依赖 molrs_cxx/xtensor，纯 install 的库对下游可用性有限；完整的 export set / `find_package(Atomiverse)` 配置**不做**。本篇的硬要求是 `bin/atomiverse` 落地。

### 3. 运行时解释器发现（本篇真正的设计面）

**选定的单一机制：解析 PATH 上的 `python3` 并把它作为 `PyConfig.executable`。** 顺序：

1. 环境变量 **`ATV_PYTHON`**（指向某个 python 可执行文件的路径）——PYTHONHOME 式兜底，服务奇异布局；
2. 否则在 `PATH` 上搜 `python3`（再 `python`）——venv 激活后即 venv 的解释器，conda 同理，Lmod module-loaded Python 同理;
3. 都没有 → 致命错误（列出查找顺序 + 提示 `ATV_PYTHON`），**绝不**回退构建期路径（no-fallback 规则）。

初始化用 `PyConfig_InitPythonConfig`（非 Isolated——要尊重用户环境），只设 `config.program_name` 与 `config.executable` 为解析到的路径，**不手工设 `config.home`**：CPython 的 getpath（PEP 587）会按该解释器自身的语义计算 `sys.prefix`/`sys.exec_prefix`，包括 venv 的 `pyvenv.cfg` 解析（venv 的 site-packages + base 的 stdlib）——手工设 home 到 venv 前缀反而因 venv 无 stdlib 而挂在 `encodings` 导入上。这替换掉 04 的构建期 FindPython3 路径假设。

**ABI 守卫**：launcher 链接的是构建期的 `libpython3.X`，只能承载同 minor 版本的环境。启动时（rank 0、`Py_InitializeFromConfig` 之前、一次性 `popen`）查询解析到的解释器 `sys.version_info[:2]`，与编译期 `PY_VERSION_HEX` 比对；不匹配 → 致命错误，同时报出两个版本号与 `ATV_PYTHON` 提示。

**libpython 链接**：正常动态链接 `Python3::Python`（Development.Embed），target 加 `INSTALL_RPATH "$ORIGIN/../lib"`——覆盖 conda（`env/bin/atomiverse` 旁即 `env/lib/libpython3.X.so`）与常见前缀安装；module-loaded Python（本机 Lmod 场景）由模块自身的 `LD_LIBRARY_PATH` 供给，在文档中写明。不用 dlopen 方案（复杂度不换收益）。

**可测性拆分**：发现逻辑落成独立类 `InterpreterLocator`（与 04 的 `PythonHost` 同命名空间；查找即构造函数，遵守 no-factory 规则；构造参数注入 override 值与 PATH 字符串），**头/源不含 `Python.h`、不执行任何解释器**——纯路径解析 + 版本字符串比对，因此可进 `tests/`（秒级、零第三方读取，符合 tests/ 硬规则）。实际的 `popen` 查询与 `PyConfig` 调用留在 `PythonHost`，属验收级运行时行为。单元测试 target 直接编译 `app/atomiverse/interpreter_locator.cpp` 这个 TU（app/ 不是库的一部分）。

### 4. 文档

`docs/user-guide/launcher.md`（新）：单机 `python atv.py` vs `mpirun -np N atomiverse atv.py` / `srun -n N atomiverse atv.py` 对照；环境要求（PATH 上有匹配 minor 版本的 `python3`、libpython 可见性 / `LD_LIBRARY_PATH`、Lmod 场景）；`ATV_PYTHON` 覆盖；两条安装路径（`ATV_BUILD_LAUNCHER=ON pip install .` 与 `cmake --install`）；wheel 不可移植（MPI + libpython 绑定环境）的说明。`docs/user-guide/index.md` 加链接。

## Files to create or modify

- `python/pyproject.toml` — 新增 `ATV_BUILD_LAUNCHER` env-override define（默认 OFF）
- `CMakeLists.txt` — `GNUInstallDirs` + `NOT SKBUILD` 守卫的库/头文件 install 规则；launcher-需-MPI 的 `FATAL_ERROR` 守卫（若 04 未含）
- `app/atomiverse/CMakeLists.txt` (04-embed) — launcher 的 SKBUILD/BINDIR 双路 install、`INSTALL_RPATH "$ORIGIN/../lib"`、locator 源文件
- `app/atomiverse/interpreter_locator.h` (new) — `InterpreterLocator`（Python.h-free）
- `app/atomiverse/interpreter_locator.cpp` (new)
- `app/atomiverse/python_host.h` (04-embed) — 接入运行时发现的接口调整
- `app/atomiverse/python_host.cpp` (04-embed) — `PyConfig_InitPythonConfig` + `config.executable`；删除构建期路径假设；ABI minor 版本守卫
- `tests/app/test_interpreter_locator.cpp` (new) — 单元测试（注入式，不执行 Python）
- `tests/CMakeLists.txt` — 注册上述测试（直接编译 locator TU）
- `docs/user-guide/launcher.md` (new)
- `docs/user-guide/index.md` — 加链接

（`app/atomiverse/*` 路径按 04 链节落地布局；本 dev 快照中尚不可见，已用 Glob 核实并标注。）

## Tasks

- [x] Write failing tests for InterpreterLocator (tests/app/test_interpreter_locator.cpp, registered in tests/CMakeLists.txt): ATV_PYTHON override 优先、PATH 逐目录搜 python3/python、找不到即 throw、版本字符串不匹配即 throw——全部注入字符串，绝不执行解释器
- [x] Implement InterpreterLocator in app/atomiverse/interpreter_locator.h/.cpp（查找即构造函数，Python.h-free，no-fallback：不支持布局 throw 且报文含 ATV_PYTHON 提示）
- [x] Rewire PythonHost (app/atomiverse/python_host.h/.cpp) to runtime discovery：PyConfig_InitPythonConfig + program_name/executable 取自 locator，一次性 popen 查询 sys.version_info 做 ABI minor 守卫，删除 04 的构建期 FindPython3 路径假设
- [x] Implement install rules in app/atomiverse/CMakeLists.txt：SKBUILD → `${SKBUILD_SCRIPTS_DIR}`、否则 `${CMAKE_INSTALL_BINDIR}`；`INSTALL_RPATH "$ORIGIN/../lib"`
- [x] Implement plain-install rules in CMakeLists.txt：`NOT SKBUILD` 守卫的 libatomiverse + include/atomiverse 安装；launcher-ON-而-MPI-缺失 的 FATAL_ERROR 守卫
- [x] Add ATV_BUILD_LAUNCHER env-override define (default OFF) to python/pyproject.toml
- [x] Add docstring per doxygen（InterpreterLocator / PythonHost 公开符号的 @brief/@param/@return；路径/环境变量语义写进 @param）
- [x] Write docs/user-guide/launcher.md and link it from docs/user-guide/index.md（python vs mpirun/srun、环境要求、ATV_PYTHON、LD_LIBRARY_PATH/Lmod、两条安装路径、wheel 不可移植说明）
- [x] Hygiene cleanup（rides the fix(neighbor) cycle; ruff clean on touched scripts）
- [x] Verify installability in a fresh venv：`ATV_BUILD_LAUNCHER=ON pip install .` → wheel 的 `.data/scripts/atomiverse` 为 ELF、PATH 上可执行、`atomiverse script.py` 单 rank 导入 venv 的 site-packages；`cmake --install --prefix <scratch>` 落出 bin/atomiverse（+ lib/include）——验收级手动运行，不进 ctest
- [x] Run full check + test suite（168 tests green incl. 13/13 locator; parity re-verified post neighbor-fix）

## Testing strategy

- **单元（tests/，秒级，唯一进 ctest 的部分）**：`InterpreterLocator` 的解析次序（override > PATH）、PATH 逐目录语义、缺失 → throw、版本比对 → throw。全部通过构造参数注入伪 PATH/伪版本字符串完成；**测试不 popen、不 import、不读任何第三方包**（tests/ 硬规则）。
- **边界**：空 PATH；`ATV_PYTHON` 指向不存在的文件（throw，报文点名该变量）；版本串 "3.11" vs 编译期 "3.12"（throw，报文含两个版本号）。
- **验收级运行时（手动，需网络 + fresh venv，非 ctest）**：wheel 安装后 `command -v atomiverse` + `file` 判 ELF；`atomiverse script.py` 中脚本打印 `sys.prefix` 等于 venv 前缀并成功 `import atomiverse`；`PATH` 剥离 python3 后启动 → 非零退出 + 清晰报错；`cmake --install` 布局检查。这些是 acceptance.md 中 runtime 型判据的执行方式。
- 无新 regressions/ 用例：链的物理回归归 04（parity case）；本篇交付面是安装与发现，无可计算观测量。

## Out of scope

- CMake export set / `find_package(Atomiverse)` 包配置、pkg-config——纯 install 只做最小布局（已在 Design §2 声明局限）。
- manylinux / PyPI 可分发的 launcher wheel（libpython + MPI 绑定环境，本篇明确默认 OFF）。
- Windows/macOS 启动器打包；dlopen-libpython 方案。
- 04 已定的协议/广播/嵌入行为本身（TaskSpec、Frame 广播、dispatch hook）——本篇不改语义，只改发现与安装。
- `python` 前缀命令包装、conda 包（conda-forge recipe）。

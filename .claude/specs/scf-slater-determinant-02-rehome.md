---
title: SCF Slater 行列式重构 · 02 疏散 scf_result.h 同住者（类型各归其家，SCFResult 留守待删）
slug: scf-slater-determinant-02-rehome
status: approved
chain: scf-slater-determinant
link: 02/04
created: 2026-07-20
---

# SCF Slater 行列式重构 · 02 疏散 scf_result.h 同住者

## Summary

`scf_result.h` 里与 `SCFResult` 同住的五名住户（`MeanFieldType`/`SpinType`/`InitGuess` 三枚举、`SCFInitialGuess`、`GradientResult`、`IterationRecord`）各自搬去语义正确的家，`SCFResult` 本体这一轮留在原文件不动，使文件末态只剩 `SCFResult` 一个类型——这是 03-unbag 删除该文件的前置条件。其中 `IterationRecord` 经消费面核实为**零读者的写死数据**（criterion-unify 之后每迭代判定量已由 `ScfConverged` + driver 自有迭代状态承载），判定冗余、直接删除。行为零变化：纯搬迁 + include 修正，编译绿 + 全测绿即完成。豁免点名：`scf::toy::SCFResult`（`solve.h` 独立 toy 命名空间）不动；NDDO semiempirical 自有类型不动。

## Domain basis

无物理。领域锚点是**行为保持**：搬迁前后所有 SCF/gradient/MP2/FCI/cube-writer 路径逐位一致；唯一的语义变动是删除一个全树零读者的诊断向量（`SCFResult::trace`），不被任何测试、示例、回归或文档观测。

## Design

**验证内部（已逐一读文件核实）：**

- `scf_result.h:16-71` 现有六个类型：`IterationRecord`(:16)、三枚举(:23-25)、`SCFInitialGuess`(:28)、`SCFResult`(:36，含 `std::vector<IterationRecord> trace` :61)、`GradientResult`(:68)。
- **直接 include 面恰好八处**（librarian 结论已复核，全树 grep 一致）：`mean_field_method.h:42`（六类型全用）、`internal/hamiltonian_builder.h:16`（**仅** `MeanFieldType`/`SpinType`）、`internal/gradient.h:20`（枚举 + `GradientResult` + `SCFResult`）、`internal/integral_derivatives.h:21`（枚举 + `SCFResult`）、`internal/mp2_pack.h:7` 与 `internal/fci_pack.h:9`（仅 `SCFResult`，**不动**）、`src/cpu/recorder/density_cube_writer.cpp:9`（`SCFResult` + `SpinType`）、`src/cpu/recorder/orbital_cube_writer.cpp:9`（仅 `SCFResult`，**不动**）。其余消费者（`rhf/uhf/rks/uks.cpp`、`orbital_info_observer.h`、`test_method_driver_parity`、`delta_scf_reference.h`、`examples/13_mp2`、`14_fci`）全部经 `mean_field_method.h`/`methods.h` 传递可达，include 面零改动。

**1. 三枚举 → 新公共头 `scf_enums.h`（不复活 `scf_types.h`）。** 请求首选的公共 `scf_types.h` 经核实**不适合**：它现在是一块活的 `#error` 墓碑（"SCFConfig was deleted…"），其存在意义就是给树外陈旧 include 一条清晰死路；复活它会把失败模式从明确报错变成名字查找错误，且与 `internal/scf_types.h`（`SpinResolution`）同名互扰——请求自己都特意提醒区分。新家 `include/atomiverse/cpu/scf/scf_enums.h`（new）：只装三枚举，零依赖，`internal/hamiltonian_builder.h` 这类内部头可以只拿枚举而不再拖入结果类型。

**2. `SCFInitialGuess` → `scf_guess.h`。** 已读现文件：`scf_guess.h` 是初猜公共头（`apply_homo_lumo_mixing`），主题吻合；搬入后补 `#include <xtensor/containers/xtensor.hpp>`。消费者（`mean_field_method.h/.cpp`、`tests/cpu/scf/delta_scf_reference.h` 经 `methods.h`）传递可达。

**3. `GradientResult` → 新公共头 `gradient_result.h`（公共可见性保持）。** 公共 `gradient.h` 同样是 `#error` 墓碑（"not public — use methods.h"），复活理由与 1 同样不成立；`internal/gradient.h` 是 internal，不满足"公共可见性必须保持"。新家 `include/atomiverse/cpu/scf/gradient_result.h`（new）：只装 `GradientResult`（forces 公共 API 返回类型，eV/Å），`mean_field_method.h` 与 `internal/gradient.h` 都 include 它。

**4. `IterationRecord` 判定冗余，删除。** 证据：全树（include/src/tests/examples/docs）唯一写点是 `mean_field_method.cpp:450` 的 `res.trace.push_back(...)`，**读者为零**——没有任何测试、示例、回归、文档读 `SCFResult::trace` 或 `IterationRecord`。criterion-unify（06 环）后，每迭代的 ΔE/DIIS 误差由 `ScfConverged` 自持成员测量（`scf_converged.h:83-96`）、迭代状态由 CPUDriver 拥有，逐迭代历史无人消费。处置：删 `IterationRecord` 结构体、删 `SCFResult::trace` 成员（其唯一消费者）、删 `mean_field_method.cpp:450-451` 的 push_back。这是请求"冗余则删除并迁移其消费者"分支的直接执行，`SCFResult` 仍留在原文件（成员删除≠搬家，03 前置目标"文件只剩 SCFResult"照常达成）。

**5. `scf_result.h` 末态 + include 修正。** 文件只剩 `SCFResult`（去掉 `trace` 后不再用 `<vector>`，一并删该 include；`@file` 注释同步）。include 面机械修正：`internal/hamiltonian_builder.h:16` 改指 `scf_enums.h`（唯一"只为枚举"的整替换点）；`mean_field_method.h`、`internal/gradient.h`、`internal/integral_derivatives.h`、`density_cube_writer.cpp` 在保留 `scf_result.h`（仍用 `SCFResult`）之外**补** `scf_enums.h`（及各自所需的 `gradient_result.h`/`scf_guess.h`），不再靠 `scf_result.h` 传递提供枚举。

**Reuse decision：**
- `reuse` `scf_guess.h` —— `SCFInitialGuess` 的现成主题头，不新建。
- `new` `scf_enums.h`、`gradient_result.h` —— 两块墓碑（公共 `scf_types.h`/`gradient.h`）保持墓碑，不复活。
- `delete` `IterationRecord` + `SCFResult::trace` + 唯一写点 —— 零读者。
- **不动**：`internal/mp2_pack.h`、`internal/fci_pack.h`、`orbital_cube_writer.cpp`、`orbital_info_observer.h`、tests/examples 的传递消费者、`scf::toy::SCFResult`（`solve.h`）、NDDO semiempirical、docs 链接（03 删文件时一并处理）。

## Files to create or modify

- `include/atomiverse/cpu/scf/scf_enums.h` (new) — `MeanFieldType`/`SpinType`/`InitGuess` 逐字搬入。
- `include/atomiverse/cpu/scf/gradient_result.h` (new) — `GradientResult` 逐字搬入。
- `include/atomiverse/cpu/scf/scf_guess.h` — 接纳 `SCFInitialGuess`（+xtensor include）。
- `include/atomiverse/cpu/scf/scf_result.h` — 删五同住者，末态只剩 `SCFResult`（去 `trace` 成员与 `<vector>`）。
- `include/atomiverse/cpu/scf/mean_field_method.h` — 补 include `scf_enums.h`/`scf_guess.h`/`gradient_result.h`。
- `src/cpu/scf/mean_field_method.cpp` — 删 `res.trace.push_back`（:450-451）。
- `include/atomiverse/cpu/scf/internal/hamiltonian_builder.h` — `:16` 的 `scf_result.h` → `scf_enums.h`。
- `include/atomiverse/cpu/scf/internal/gradient.h` — 补 `scf_enums.h` + `gradient_result.h`（保留 `scf_result.h`）。
- `include/atomiverse/cpu/scf/internal/integral_derivatives.h` — 补 `scf_enums.h`（保留 `scf_result.h`，注释同步）。
- `src/cpu/recorder/density_cube_writer.cpp` — 补 `scf_enums.h`（`SpinType`）。
- `tests/cpu/scf/test_scf_type_homes.cpp` (new) — 新家 include 面单元测试。
- `tests/CMakeLists.txt` — 注册 `test_scf_type_homes`。

## Tasks

- [ ] Write failing include-surface test tests/cpu/scf/test_scf_type_homes.cpp (includes ONLY scf_enums.h / scf_guess.h / gradient_result.h; uses all three enums, constructs SCFInitialGuess and GradientResult) and register it in tests/CMakeLists.txt
- [ ] Implement include/atomiverse/cpu/scf/scf_enums.h with MeanFieldType/SpinType/InitGuess moved verbatim from scf_result.h, doxygen per doc.style
- [ ] Implement include/atomiverse/cpu/scf/gradient_result.h with GradientResult moved verbatim (forces [eV/Å], energy [eV] doc retained)
- [ ] Move SCFInitialGuess verbatim into include/atomiverse/cpu/scf/scf_guess.h with the xtensor include added
- [ ] Delete IterationRecord, the SCFResult::trace member, and the res.trace.push_back at src/cpu/scf/mean_field_method.cpp:450-451 (redundant: zero readers tree-wide)
- [ ] Trim include/atomiverse/cpu/scf/scf_result.h to SCFResult-only (drop <vector>, update @file doc) and fix the direct includers: mean_field_method.h adds scf_enums.h/scf_guess.h/gradient_result.h; internal/hamiltonian_builder.h:16 switches to scf_enums.h; internal/gradient.h and internal/integral_derivatives.h add scf_enums.h (+gradient_result.h where used); src/cpu/recorder/density_cube_writer.cpp adds scf_enums.h
- [ ] Verify end state by grep: scf_result.h defines only SCFResult; IterationRecord count is zero tree-wide; scf::toy::SCFResult and semiempirical types untouched
- [ ] Run full check + test suite

## Testing strategy

- **新家 include 面（单元，秒级）**：`test_scf_type_homes` 只 include 三个新/改后的公共头即可用全部搬迁类型——枚举可取值、`SCFInitialGuess`/`GradientResult` 可构造——不 include `scf_result.h`。这正是 03 删文件所依赖的性质。
- **行为保持（既有套件即守卫）**：`test_method_driver_parity`、`delta_scf_reference.h` 驱动的 ΔSCF 测试、SCF stages 测试、cube writer 测试、`examples/13_mp2`/`14_fci` 编译运行全部不改一行断言、全绿；`trace` 删除无任何断言观测（零读者已核实）。
- **边界情况**：`internal/hamiltonian_builder.h` 换头后独立可编译（枚举不再经 `scf_result.h` 传递）；`density_cube_writer.cpp` 的 `SpinType` 分支照常。
- **Domain validation**：无物理——编译绿 + 全测绿 + 终态 grep 即验证（行为零变化规格）。
- **Regression example**：本环不含新 `regressions/` 案例——纯搬迁，非交付特性；既有 SCF 回归即行为守卫。

## Out of scope

- `SCFResult` 本体搬迁与其消费者迁移、`scf_result.h` 文件删除 —— 03-unbag。
- `build_*` 及方法侧重构 —— 04-methods。
- `SlaterDeterminant` 类型 —— 01-type（链内顺序前置，无代码依赖，02 不消费它）。
- 两块 `#error` 墓碑（公共 `scf_types.h`/`gradient.h`）的最终删除，及 docs 里指向 `scf_types.h`/`scf_result.h` 的陈旧链接 —— 随 03 文件删除一并处理。
- `scf::toy::SCFResult`（`solve.h`）与 NDDO semiempirical 自有类型 —— 点名豁免，不动。

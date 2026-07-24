---
slug: scf-review-05-api-hygiene
title: "SCF API 卫生整顿 — 网格构造器、rys_roots 抛异常、π 常量归一、内部 .cpp 归位"
status: done
created: 2026-07-24
revised: 2026-07-24
domain: scf / api-hygiene
queue: scf-review
depends_on: [scf-review-01-linalg]
grilled: true
---

# SCF API 卫生整顿(scf-review-05-api-hygiene)

## Summary

本规格收拢 SCF 审阅中五处 API 卫生问题,除一处外均为行为保持(behavior-preserving)的整顿。将 `cpu/scf/grid.h` 里作为唯一构造路径的工厂 `build_molecular_grid` 转为 `MolecularGrid` 的构造函数(遵循 `GthPotential(Element,name)` 查表构造器的既有先例),并同步移除 `no_factory_scrub.sh` 对它的白名单、刷新 `compatibility_surface_scrub.sh` 的隐式-SG0 守卫;把 `gf2.cpp` 里触发工厂命名约定的 `build_se_part` 重命名;把 `rys_roots.cpp` 中域内唯一的 `abort()` 改为抛异常并新增一条断言测试;将 7 处本地 `kPi/PI` 常量与 5 处 `M_PI` 归一到 `units::pi`;并把 11 个错放在 `src/cpu/scf/` 根、实则实现内部头的 `.cpp` 归位到镜像的 `src/cpu/scf/internal/`。对用户可见的结果是:SCF 域不再有工厂式主构造路径、不再有不可捕获的 `abort`、π 常量单一来源、源码树布局与头文件镜像一致;所有既有 SCF/DFT 数值与回归保持不变。

## Domain basis(数值等价性,非新物理)

本规格不引入新物理,不改变任何科学数值。相关的数值不变量已逐位核对:

- π 常量:所有被清扫的本地常量(`xc_gga.cpp:35`、`grid.cpp:26`、`rys_eri.cpp:43`、`native_integrals.cpp:23`、`internal/pcm.cpp:26`、`dftd4.cpp:54`、`xc_functional.cpp:23`)字面均为 `3.14159265358979323846`,与 `include/atomiverse/core/units.h:118` 的 `units::pi` **逐位相同**;`M_PI`(glibc `math.h`)同样是该字面值。因此归一为纯文本重构,数值逐位不变。
- 错误路径:`rys_roots.cpp` 的 `abort` 改 `throw` 携带**同一条消息**(`nroots=%d not supported (max 7)` / `rys_roots fails: nroots=%d x=%g`),对所有 `nroots∈{1..7}` 的成功路径零改动 —— 仅把不可达/不可捕获的失败终止变为可捕获、可测试的异常。
- 网格构造器:`n_radial`/`n_angular` 仅对 `GridLevel::SG0` 生效的语义**原样保留**(required-then-defaulted 形参序 `MolecularGrid(atoms, natm, level, n_radial=50, n_angular=110)`)。

## Design

### 依赖与链顺序

本规格属 8 段链的第 5 段,依赖 `scf-review-01-linalg`(链序保证其先落地)。原因仅在于:item 13 用 `git mv` 搬动的 11 个 `.cpp` 中,`orbital_hessian.cpp`/`cphf_pack.cpp`/`ucphf_pack.cpp` 三个文件的**内容**由 01 段修改;本段只搬动这三个已被 01 修改过的文件,故须在 01 之后执行。除此之外无跨段耦合。

### item 10 — 网格工厂转构造器(🔴 arch)

`include/atomiverse/cpu/scf/grid.h:38` 的 `MolecularGrid` 是无构造函数的聚合体,`grid.h:79` 的自由函数 `build_molecular_grid(atoms, natm, level, n_radial=50, n_angular=110)` 是它**唯一**的构造路径 —— 这正是"工厂即主构造器"的违规。改为在 `struct MolecularGrid` 上声明构造函数:

```cpp
MolecularGrid(const native::Atom* atoms, int natm, GridLevel level,
              int n_radial = 50, int n_angular = 110);
```

保持 required-then-defaulted 形参序与 SG0-only 的 `n_radial`/`n_angular` 语义。这与 `GthPotential(molrs::Element, name)` 查表构造器先例一致(查表构造器是允许的,不是工厂)。`grid.cpp:687` 的定义随之从自由函数改为成员构造器。

**门禁 lockstep(铁律)**:`build_molecular_grid` 目前被两处架构 scrub 显式豁免,必须在同一次变更中同步收紧,否则门禁与实现漂移:
- `tests/architecture/no_factory_scrub.sh:116` 的白名单行(在 `grep -vE 'build_molecular_grid|…'` 中)删除其 token —— 删除后该门禁转为**强制** SCF 树中不存在 `build_molecular_grid` 自由函数。
- `tests/cpu/compatibility_surface_scrub.sh:167-177` 的 awk 守卫当前匹配 `MolecularGrid build_molecular_grid(` 且次行为 `int n_radial`(拦截隐式-SG0 重载)。刷新为构造器形态:匹配构造器声明并断言 `natm` 之后是显式 `GridLevel`(次行不得直接是 `int n_radial`),**保留其原意**(不得出现省略 `GridLevel` 的隐式-SG0 构造)。

### item 11 — 重命名 `build_se_part` 与 adc 的三个 `build_m*_`(🟢，audit 扩充)

`gf2.cpp:94` 的 `AuxSlice build_se_part(...)` 抽取本身是合理的(两处调用 `gf2.cpp:418-419`),但 `build_` 前缀触发 `no_factory_scrub.sh` 的 `cpu/scf` 树 `build_*` 通配门禁(该门禁对 `cpu/scf` 递归 grep `build_*(`,且 `build_se_part` **不在**白名单内)。**选定新名 `compute_se_part`**(动词 `compute`,非工厂前缀,保留 `_se_part` 后缀以维持与 PySCF `ragf2_slow.build_se_part` 的可读映射)。改定义 + 两处调用点。

**Audit 扩充（grill 2026-07-24，实测两条 scrub 在 dev 上均为红）**：`no_factory_scrub` 的现存违规不止 `build_se_part`——`adc.cpp` 的三个文件局部助手 `build_mss_`(:178, 调用 :1274)、`build_mab_ea_`(:880, 调用 :1539)、`build_mij_ip_`(:1312, 调用 :1489) 同被点名。同一机械处理：重命名为 `compute_mss_` / `compute_mab_ea_` / `compute_mij_ip_`(定义 + 各调用点，纯文本改名，无语义变化)。(`build_imds` 的 3 处违规由 scf-review-03 以 `EomImds` 计算构造函数偿清。) 此项与 item 10 的白名单删除一并,是 `no_factory_scrub` 在本段之后**恢复绿色**的必要条件(见测试策略中的门禁说明)。

**Audit 追加 — item 15：compat scrub 的 gpw `molrs::Element` 违规修复**。`compatibility_surface_scrub.sh` 在 dev 上亦为红：`using molrs::Element;` 别名违反"Element 必须以 `molrs::Element` 使用、不得设 ATV 别名/适配器"门禁,现存 4 个代码位点(`tests/cpu/gpw/test_cpu_gpw_dftd3.cpp:209,258`、`regressions/gpw/gpw-pbe-d3bj-ice/reg_gpw-pbe-d3bj-ice.cpp:23`、`regressions/gpw/gpw-d4-water/reg_gpw-d4-water.cpp:23`)+ 2 处文档(`docs/tutorial/methods/gpw-pbe-d3bj-ice.md:93`、`gpw-d4-water.md:109`)。修复为逐点删除 `using` 别名、改用限定名 `molrs::Element`(机械单行改动 ×6,跨入 gpw 域但属 iron-law 就地偿债——不修则本段 ac-003 的"exit 0"永不可验证)。

### item 14 — `abort` 改 `throw`(🔴 janitor)

`src/cpu/scf/rys_roots.cpp:96-105` 的 `default:`(不支持的 `nroots`)与通用错误分支调用 `std::fprintf(stderr,…)+std::abort()`,是整个 SCF 域**唯一**的 `abort`(数百处均 `throw`);且经 g/h 壳层(`nroots>7`)可达,不可捕获、不可测。修复:两处 `abort` 均替换为携带同一消息的 `throw std::invalid_argument`(不支持的 `nroots`)/ `throw std::runtime_error`(迭代求根返回非零 `err`),删去紧邻的两个 `fprintf`(消息随异常传递)。`rys_roots_iterative.cpp` 的 `err` 返回契约保持不变 —— 其唯一调用者(`rys_roots.cpp`)现在把 `err` 表面化为异常而非 `abort`,原本"死胡同"的 `err` 管路就此贯通。

**交互守卫**:延后规格 `rsh-eri-cam-b3lyp` 将扩展 `rys_roots`/`rys_eri` 以支持 erf-衰减 ERI;本段**只改错误处理,不重绘任何签名**。

### item 9 — π 常量归一(🟡+🟢)

7 处本地 `constexpr kPi/PI`(见 Domain basis 列举)+ 5 处 `M_PI`(`one_electron.cpp:261`、`integral_derivatives.cpp:2138`、`internal/cart2sph.cpp:115`、`ecp_integrals.cpp:40,50`)全部改用 `units::pi`(`#include "atomiverse/core/units.h"`),删除本地常量。纯文本重构,逐位等价(已在 Domain basis 核对字面)。

### item 13 — 内部 .cpp 归位(🟢 arch)

11 个 `.cpp` 实现的是 `include/atomiverse/cpu/scf/internal/` 下的头文件(11 个内部头均已存在于 `internal/`),但 `.cpp` 错放在 `src/cpu/scf/` 根。对每个 `git mv` 到 `src/cpu/scf/internal/`。**构建须同步**:`CMakeLists.txt` 采用**显式源清单**(非递归 glob),第 493/495/497/502/529/530/546/550/551/552/553 行分别列出这 11 个文件;每行路径前缀改为 `src/cpu/scf/internal/`。搬动后若某 `.cpp` 以同目录相对 `#include` 引用其内部兄弟头,须相应修正(实现时逐一核对)。注意 `one_electron.cpp`/`integral_derivatives.cpp`/`ecp_integrals.cpp` 同时被 item 9 编辑、`gradient.cpp` 同时被 item 10 编辑 —— 均为同一文件的不同处改动,先内容后 `git mv`(或反之,`git mv` 保留历史,顺序不影响结果)。

### Reuse decision

未提供 `librarian_report`(caller 标注 blueprint 刷新延后),无复用候选需裁决。设计直接沿用两条 house 先例:构造器路径遵 `GthPotential(Element,name)` 查表构造器;π 单一来源为既有 `units::pi`。不引入任何新符号,不新建类型。

## Files to create or modify

item 10(网格构造器 + 门禁 lockstep):
- `include/atomiverse/cpu/scf/grid.h` — 工厂声明改为 `MolecularGrid` 构造器
- `src/cpu/scf/grid.cpp` — 工厂定义(:687)改为构造器(另见 item 9 :26)
- `src/cpu/scf/tdhf.cpp` — 4 处调用点(:547/838/909/1339)
- `src/cpu/scf/mean_field_method.cpp` — 1 处调用点(:792)
- `src/cpu/scf/internal/hamiltonian_builder.cpp` — 1 处调用点(:88)
- `tests/cpu/scf/test_xc_reference.cpp` — 调用点(:483/507/527)
- `tests/cpu/scf/test_init_guess.cpp` — 调用点(:152/180/213/233)
- `tests/cpu/scf/test_full_grid_alda.cpp` — 调用点(:41/86)
- `tests/architecture/no_factory_scrub.sh` — 删除 `build_molecular_grid` 白名单 token(:116)
- `tests/cpu/compatibility_surface_scrub.sh` — 刷新 grid 守卫为构造器形态(:167-177)

item 11(audit 扩充):
- `src/cpu/scf/gf2.cpp` — `build_se_part`→`compute_se_part`(定义 :94 + 调用 :418/:419)
- `src/cpu/scf/adc.cpp` — `build_mss_`→`compute_mss_`、`build_mab_ea_`→`compute_mab_ea_`、`build_mij_ip_`→`compute_mij_ip_`(定义 :178/:880/:1312 + 调用 :1274/:1539/:1489)

item 15(audit 追加,compat scrub 偿债):
- `tests/cpu/gpw/test_cpu_gpw_dftd3.cpp` — 删 `using molrs::Element;`(:209/:258),改限定名
- `regressions/gpw/gpw-pbe-d3bj-ice/reg_gpw-pbe-d3bj-ice.cpp` — 同上(:23)
- `regressions/gpw/gpw-d4-water/reg_gpw-d4-water.cpp` — 同上(:23)
- `docs/tutorial/methods/gpw-pbe-d3bj-ice.md`、`docs/tutorial/methods/gpw-d4-water.md` — 文档代码块同步(:93/:109)

item 14:
- `src/cpu/scf/rys_roots.cpp` — 两处 `abort` 改 `throw`(:96-105)
- `tests/cpu/scf/test_rys_roots.cpp` — 新增 `nroots>7` 抛异常断言测试

item 9(π 归一,尚未在上文列出的文件):
- `src/cpu/scf/xc_gga.cpp`、`src/cpu/scf/rys_eri.cpp`、`src/cpu/scf/native_integrals.cpp`、`src/cpu/scf/internal/pcm.cpp`、`src/cpu/scf/dftd4.cpp`、`src/cpu/scf/xc_functional.cpp`、`src/cpu/scf/internal/cart2sph.cpp`

item 13(归位 11 个 .cpp,`src/cpu/scf/X.cpp → src/cpu/scf/internal/X.cpp`):
- `cphf_pack.cpp`、`ucphf_pack.cpp`、`orbital_hessian.cpp`、`stability.cpp`、`ecp_integrals.cpp`(另见 item 9)、`minao.cpp`、`one_electron.cpp`(另见 item 9)、`shell_pair.cpp`、`integrals.cpp`、`integral_derivatives.cpp`(另见 item 9)、`gradient.cpp`(另见 item 10)
- `CMakeLists.txt` — 上述 11 个显式源清单条目改前缀为 `internal/`

回归:
- `regressions/layout/scf_api_hygiene/reg_scf_api_hygiene.cpp` (new) — 库外见证:`MolecularGrid` 构造器可解析/链接 + `rys_roots(8,…)` 抛可捕获异常
- `regressions/CMakeLists.txt` — 新增 `atv_add_regression(layout scf_api_hygiene)`(HF 守卫块内)

## Tasks

- [x] Write failing unit test asserting rys_roots throws for unsupported nroots (tests/cpu/scf/test_rys_roots.cpp → RysRoots)
- [x] Replace both abort sites with throws carrying the same message in src/cpu/scf/rys_roots.cpp
- [x] Convert build_molecular_grid factory to a MolecularGrid constructor (include/atomiverse/cpu/scf/grid.h + src/cpu/scf/grid.cpp), preserving SG0-only n_radial/n_angular semantics
- [x] Port the 13 build_molecular_grid call sites to the MolecularGrid constructor (tdhf.cpp, mean_field_method.cpp, internal/hamiltonian_builder.cpp, test_xc_reference.cpp, test_init_guess.cpp, test_full_grid_alda.cpp, gradient.cpp)
- [x] Update architecture scrubs in lockstep: remove build_molecular_grid whitelist in tests/architecture/no_factory_scrub.sh, refresh the implicit-SG0 guard for the constructor shape in tests/cpu/compatibility_surface_scrub.sh, and fix the 6 pre-existing `using molrs::Element` compat-scrub violations (tests/cpu/gpw/test_cpu_gpw_dftd3.cpp, regressions/gpw/gpw-pbe-d3bj-ice, regressions/gpw/gpw-d4-water, + 2 tutorial docs) so both scrubs restore to green
- [x] Rename build_se_part to compute_se_part in src/cpu/scf/gf2.cpp and build_mss_/build_mab_ea_/build_mij_ip_ to compute_mss_/compute_mab_ea_/compute_mij_ip_ in src/cpu/scf/adc.cpp (definitions + all call sites)
- [x] Sweep the 7 local kPi/PI constants and 5 M_PI uses to units::pi across the enumerated cpu/scf files, deleting the local constants
- [x] Relocate the 11 misplaced internal .cpp files under src/cpu/scf/internal/ (git mv) and repath their 11 explicit entries in CMakeLists.txt
- [x] Add regression example regressions/layout/scf_api_hygiene/reg_scf_api_hygiene.cpp (out-of-library witness: MolecularGrid ctor resolves/links, rys_roots(8,…) throws catchably; register via atv_add_regression(layout scf_api_hygiene))
- [x] Run full check + test suite (unit suite, no_factory_scrub with whitelist removed, compatibility_surface_scrub, regressions incl. existing SCF/DFT cases)

## Testing strategy

- 单元(唯一新行为,item 14):在 `tests/cpu/scf/test_rys_roots.cpp`(镜像 `src/cpu/scf/rys_roots.cpp`,`RysRoots` 套件,单函数 `native::rys_roots`)新增断言:
  - happy:`nroots∈{1..7}` 仍正常返回(既有测试已覆盖,不回退)。
  - edge/域校验:`rys_roots(8, 1.0, u, w)` 抛 `std::invalid_argument`(硬编码期望:`EXPECT_THROW(..., std::invalid_argument)`);`nroots<1` 同样抛出。异常消息含 `nroots`。
- 行为保持项(item 9/10/11/13)不新增数值单元测试 —— 由**既有 SCF/DFT 单元套件** + 两条架构 scrub + 编译共同守护:
  - `no_factory_scrub.sh`:删除 `build_molecular_grid` 白名单后必须**通过**(即现在强制 `build_molecular_grid` 不存在);item 11 的重命名同时确保 `cpu/scf` 树 `build_*` 通配门禁绿色(`build_se_part` 目前即触发该门禁,重命名清除之)。
  - `compatibility_surface_scrub.sh`:刷新后的 grid 守卫通过,且守卫指向构造器(强制 `GridLevel`)形态。
  - π 归一逐位等价(Domain basis 已核对),既有 XC/ERI 单元数值不变。
- 回归示例(架构见证型,遵 `regressions/layout/architecture_gates/reg_architecture_gates.cpp` 先例;无数值 golden、无第三方运行时):`regressions/layout/scf_api_hygiene/reg_scf_api_hygiene.cpp` 作为库外消费者 (1) 用新构造器构造一个单原子 `MolecularGrid`(见证 item 10 形态可解析并链接),(2) 调用 `native::rys_roots(8, 1.0, u, w)` 并断言其抛出、可在库外 `catch`(见证 item 14 错误路径可捕获),运行返回 0。
- 既有 SCF/DFT 回归(`regressions/scf/*`,经 RKS/UKS XC 路径间接行使网格)须复现其引用文献值**不变**,作为行为保持的端到端护栏。

## Out of scope

- `rys_roots`/`rys_eri` 的任何**签名**改动(erf-衰减 ERI 扩展保留给延后规格 `rsh-eri-cam-b3lyp`;本段只改错误处理)。
- 其余仍在 `no_factory_scrub` 白名单中的 `build_*`(`build_minao_density`、`build_minao_ao_list`、`build_shell_list_from_table`、`build_shell_to_atom_from_integrals`)—— 各属自己的 IntegralPack-删除/gradient 轮次,本段不动。
- 命名空间与目录对齐(`ATV::cpu::scf` 未与 `internal/` 目录对齐)—— 全仓延后项,非本段。
- 任何 grid/XC/ERI 的数值行为变化 —— 本段全部行为保持(仅 item 14 是"更好的错误路径")。
- `scf-review-01-linalg` 对 `orbital_hessian/cphf_pack/ucphf_pack` 的**内容**改动 —— 由 01 段负责;本段仅在 01 落地后搬动这些文件。
- 观察到但不处理:`tests/cpu/compatibility_surface_scrub.sh:28` 的 `CPU_DIRS` 含一条陈旧路径 `${ROOT}/tests/test_xc_reference.cpp`(实际文件在 `tests/cpu/scf/`,已被 `tests/cpu` 目录覆盖);无害(grep 命中空),不在本段五项范围内,记录以备后续清理。

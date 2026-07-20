---
title: SCF Slater 行列式 · 04 cpu/scf build_* 自由函数收编为成员方法 + no_factory 通配门禁
slug: scf-slater-determinant-04-methods
status: code-complete
chain: scf-slater-determinant
link: 04/04
depends: scf-slater-determinant-01-type, scf-slater-determinant-02-rehome, scf-slater-determinant-03-unbag
created: 2026-07-20
---

# scf-slater-determinant-04-methods — cpu/scf build_* 自由函数收编为成员方法 + no_factory 通配门禁

## Summary

把 cpu/scf 剩余的全部 `build_*` 自由函数收编为数据所属类的成员方法，消灭 SCF 域最后一批 factory 形自由函数：operators.h 六个 Fock 构建算子归入 `ShellPairList` / `ShellBlockPairList` /（stored-ERI 路径）`IntegralPack`，pair-list 两个 builder 变构造函数，一/二电子积分四个 builder 变 `IntegralPack` 私有成员，`build_libcint_context` 变 `LibcintContext` 构造函数。kernel 本体逐字搬移，单遍 ERI 融合、8 重对称、Schwarz 筛选一个 token 不动，零性能回归。全部转换落地后，`tests/architecture/no_factory_scrub.sh` 开闸 `build_*` 通配门禁（作用域为 cpu/scf 树）并对 13 个旧函数名保留全树黑名单防回潮；脚本在本 spec 完成态全绿。本 spec 是 scf-slater-determinant 四链（01-type → 02-rehome → 03-unbag → 04-methods）的收口。

## Domain basis

本 spec 不引入任何新物理；约束是**数值逐位不变**。被搬移 kernel 实现的物理量保持原样：

- 闭壳层 Fock 二电子部分 G = J − ½K，G_ij = Σ_kl P_kl [(ij|kl) − ½(ik|jl)]（Szabo & Ostlund, *Modern Quantum Chemistry*, Dover 1996, §3.4.4）。
- direct SCF：ERI 现算不存储、canonical 四重组单遍求值 + 8 重置换对称散射（Almlöf, Faegri, Korsell, *J. Comput. Chem.* **3** (1982) 385, DOI 10.1002/jcc.540030314）。
- Cauchy–Schwarz 筛选 |(ij|kl)| ≤ √(ij|ij)·√(kl|kl)（Häser & Ahlrichs, *J. Comput. Chem.* **10** (1989) 104, DOI 10.1002/jcc.540100111）。

单位不变：积分与能量 Hartree（内部），坐标 Bohr（私有 scratch），Frame 边界 eV/Å（`core/units.h`）。

## Design

### 归属映射（librarian 已穷举调用面，本轮逐一核实）

| 旧自由函数（声明位置） | 新成员方法 | 实现文件 |
|---|---|---|
| `build_g_rhf_direct` (operators.h:55) | `ShellPairList::g_rhf(P_total, G) const` | operators.cpp（原地改签名） |
| `build_jk_uhf_direct` (operators.h:72) | `ShellPairList::jk_uhf(P_total, P_alpha, P_beta, J, K_alpha, K_beta) const` | 同上 |
| `build_j_direct` (operators.h:82) | `ShellPairList::j(P_total, J) const` | 同上 |
| `build_g_rhf_direct_block` (operators.h:123) | `ShellBlockPairList::g_rhf(P_total, G) const` | 同上 |
| `build_k_direct_block` (operators.h:102) | `ShellBlockPairList::k(P, K) const` | 同上 |
| `build_coulomb_exchange_fused` (operators.h:39) | `IntegralPack::coulomb_exchange_fused(P_total, G) const` | 同上 |
| `build_shell_pair_list` (internal/shell_pair.h:52) | `ShellPairList::ShellPairList(const IntegralPack&, double thr = 1e-10)` | shell_pair.cpp |
| `build_shell_block_pair_list` (internal/shell_pair.h:93) | `ShellBlockPairList::ShellBlockPairList(const IntegralPack&, double thr = 1e-10)` | shell_pair.cpp |
| `build_overlap` / `build_kinetic` (internal/one_electron.h:33,38) | `IntegralPack` **私有**成员 `overlap(shells)` / `kinetic(shells)` | one_electron.cpp |
| `build_nuclear_attraction` (internal/one_electron.h:45) | `IntegralPack` **私有**成员 `nuclear_attraction(shells, atoms)` | one_electron.cpp |
| `build_eri_dense` (internal/two_electron.h:33) | `IntegralPack` **私有**成员 `eri_dense(shells)` | two_electron.cpp |
| `build_libcint_context` (libcint_bridge.h:42) | `LibcintContext::LibcintContext(const Frame&)`（`ATV_HAS_LIBCINT` 下） | libcint_bridge.cpp |

**归属拍板（读实现后）：**

- `build_coulomb_exchange_fused` 输入为 P + ERI + n，无任何 shell-pair 数据；ERI 与 n_basis 恰是 `IntegralPack` 的成员 → 归 `IntegralPack` 公有 const 成员（唯一生产调用者 hamiltonian_builder.cpp:64 手上正是 `integrals`，调用点改为 `integrals.coulomb_exchange_fused(P_alpha, G)`）。候选 ShellPairList 被否：该 kernel 是 stored-ERI 的 BLAS 收缩路径，与 pair list 零耦合。
- 一/二电子 builder 签名保留 `ShellList` / `AtomList` 参数（kernel 逐字搬移），标 `private` + `static`（不读对象状态，唯一调用者为 `IntegralPack` 构造函数 integrals.cpp:49-57）。注明：`IntegralPack` 是排期删除的过渡债（CLAUDE.md banned 表），这些方法随其替代者一同迁移，本轮是务实一跳，不扩大其公有表面。

### 逐字搬移保证（零性能回归的机制）

成员化只改签名行，不改循环体。为使 kernel 体零 token 改动可被 diff 审查验证：

- pair-list 方法体首行加 `const ShellPairList& sp = *this;`（block 版同理），原体内所有 `sp.pairs` / `sp.shells` 拼写原样保留；
- `IntegralPack::coulomb_exchange_fused` 体首行加 `const int n = n_basis;`，体内 `ERI` 直接解析到同名成员；
- 构造函数版沿用原函数体的局部 `ShellPairList out; … ` 逻辑，收尾 `*this = std::move(out);`（默认构造不递归）；`LibcintContext` 构造函数同理。

禁止顺手改动：不引入表达式模板、算符代数、任何中间物化；OpenMP thread-local 归约、canonical 序、mirror-permutation 修正逐字保留。`build_coulomb_exchange_fused` 中现有的 `xt::adapt` 零拷贝 + DGEMV 结构原样。

### 文件处置（最小 diff 方案，读文件后定）

- `include/atomiverse/cpu/scf/operators.h`：六个声明分别移入 `internal/shell_pair.h`（5 个）与 `internal/integrals.h`（1 个）后，该文件改为 `#error` 毒丸 shim——与既有 `scf/shell_pair.h` / `one_electron.h` / `two_electron.h` 毒丸同款先例。includer 仅 3 处：hamiltonian_builder.cpp、operators.cpp、test_integrals.cpp，全部改 include。
- `src/cpu/scf/operators.cpp`：**保留文件名**（CMakeLists.txt:401 不动，git 历史零 rename 噪音），改 include 为 `internal/shell_pair.h` + `internal/integrals.h`，定义改成员限定名。
- `internal/one_electron.h`：三个 build 声明迁走后仅剩 `nuclear_repulsion_energy(AtomList)`（非 build_*，不在本轮范围）——文件保留，头注释里的 `build_overlap(...)` 示例代码更新为成员拼写。
- `internal/two_electron.h`：唯一声明迁走后为空 → 删除；two_electron.cpp 改 include `internal/integrals.h`（若 CMake 有列该头则同步移除）。integrals.cpp 顶部 doc 示例同步更新。
- `internal/shell_pair.h` 的 `blocks_from_shells` 不是 build_* 名、非 factory 形（返回视图分组），不动。

### 调用面（全部仅换拼写）

- 生产：hamiltonian_builder.cpp:22-23（prepare 内 `sp_block_ = ShellBlockPairList{integrals}; sp_ = ShellPairList{integrals};`）、:64,66,73,127,129,141-142,158,169（成员持有 `sp_`/`sp_block_`，改 `sp_.g_rhf(…)` 等拼写）；integral_derivatives.cpp:689（`LibcintContext cint_ctx{frame};`）、:1001（`const ShellBlockPairList bpl{integrals, 1e-10};`）；integrals.cpp:49-57（私有成员拼写）。
- 测试：test_integrals.cpp:389-394（构造函数 + `g_rhf` 成员拼写，去掉 operators.h include；AO 级 vs 块级 parity 断言不变，无需新 fixture）；test_rys_eri.cpp:187,197,211,222,246,261（构造函数拼写）。

### 门禁设计（no_factory_scrub.sh）

1. **旧名全树黑名单**（`FACTORY_PATS` 逐个追加，作用域沿用现有六目录）：`build_g_rhf_direct`、`build_g_rhf_direct_block`、`build_jk_uhf_direct`、`build_j_direct`、`build_k_direct_block`、`build_coulomb_exchange_fused`、`build_shell_pair_list`、`build_shell_block_pair_list`、`build_overlap`、`build_kinetic`、`build_nuclear_attraction`、`build_eri_dense`、`build_libcint_context`（各 `'name[[:space:]]*\('`）。新成员名均去掉 `build_` 前缀，不与黑名单冲突。
2. **`build_*` 通配开闸**：新增 pattern `'build_[a-z_0-9]*[[:space:]]*\('`，过滤器：注释行（沿用现有 regex）、成员/限定访问（`grep -vE '(\.|->|::)build_'`，即请求白名单 `.build(`/`->build(` 的一般化——同时覆盖 `linalg::build_density_matrix(` 这类跨域限定调用）。
3. **作用域裁定（对请求的显式偏离）**：通配 pattern 作用域为 **cpu/scf 树**（`include/atomiverse/cpu/scf` + `src/cpu/scf`），而非全部六目录。依据：全树审计发现约 40 个合法 `build_*` **成员方法声明**（`Frame::build_special_list`、`Backend::build_components/build_plan`、`CellList::build_half/full/stencil`、CUDA `build_halo/build_tag_maps`、NDDO `build_packed/build_one_center` 等）——声明行无法用 grep 与自由函数区分，全树通配需要约 40 项白名单，门禁即被打穿。scf 域通配 + 旧名全树黑名单同时保住两个约束：cpu/scf 不能再进任何新 build_* 自由函数；已转换旧名任何地方不得回潮。`include/atomiverse/cpu/linalg/density.h` 的 `build_density_matrix`（GPW 六调用点共用，随 GPW 轮次收编）因此天然在通配作用域外——脚本中保留带注释的豁免说明，防止后续把通配扩域时误伤。
4. **scf 域残留豁免白名单**（每项必须带注释说明 + 归属轮次）：`build_molecular_grid`（grid.h:79，公有网格 builder，随 IntegralPack 删除轮收编）、`build_minao_density` / `build_minao_ao_list`（minao 初猜机器，同轮）、`build_shell_list_from_table`（basis_set.cpp 文件局部 helper，同轮）、`build_shell_to_atom_from_integrals`（integral_derivatives.cpp 内部，gradient 轮）。这些不在本 spec 转换清单（librarian 清单如实为界），豁免防止门禁把本轮范围逼爆。
5. 脚本在本 spec 完成态**全绿**，且黑名单可证伪（注入旧名临时文件应使脚本非零退出）。

## Files to create or modify

- include/atomiverse/cpu/scf/internal/shell_pair.h
- src/cpu/scf/shell_pair.cpp
- include/atomiverse/cpu/scf/operators.h（改 #error 毒丸 shim）
- src/cpu/scf/operators.cpp
- include/atomiverse/cpu/scf/internal/integrals.h
- src/cpu/scf/integrals.cpp
- include/atomiverse/cpu/scf/internal/one_electron.h
- src/cpu/scf/one_electron.cpp
- include/atomiverse/cpu/scf/internal/two_electron.h（删除）
- src/cpu/scf/two_electron.cpp
- include/atomiverse/cpu/scf/libcint_bridge.h
- src/cpu/scf/libcint_bridge.cpp
- src/cpu/scf/integral_derivatives.cpp
- src/cpu/scf/internal/hamiltonian_builder.cpp
- tests/cpu/scf/test_integrals.cpp
- tests/cpu/scf/test_rys_eri.cpp
- tests/architecture/no_factory_scrub.sh

## Tasks

- [x] Write the RED state: update tests/cpu/scf/test_integrals.cpp:389-394 and tests/cpu/scf/test_rys_eri.cpp:187,197,211,222,246,261 to member/constructor spellings and drop the operators.h include (build fails until members exist; parity assertions unchanged, no new fixture)
- [x] Implement ShellPairList / ShellBlockPairList constructors from (IntegralPack, thr) in include/atomiverse/cpu/scf/internal/shell_pair.h + src/cpu/scf/shell_pair.cpp, keeping default ctors and moving kernel bodies verbatim via the `out` + move-assign idiom
- [x] Implement the five pair-list Fock members (g_rhf / jk_uhf / j on ShellPairList; g_rhf / k on ShellBlockPairList) plus IntegralPack::coulomb_exchange_fused in src/cpu/scf/operators.cpp with `*this`-alias verbatim bodies; declarations into internal/shell_pair.h and internal/integrals.h; turn include/atomiverse/cpu/scf/operators.h into the #error shim
- [x] Implement IntegralPack private static members overlap / kinetic / nuclear_attraction / eri_dense (declarations in internal/integrals.h; definitions stay in src/cpu/scf/one_electron.cpp and src/cpu/scf/two_electron.cpp); update src/cpu/scf/integrals.cpp:49-57 call spellings; delete include/atomiverse/cpu/scf/internal/two_electron.h; refresh doc snippets in internal/one_electron.h and integrals.cpp
- [x] Implement LibcintContext::LibcintContext(const Frame&) in include/atomiverse/cpu/scf/libcint_bridge.h + src/cpu/scf/libcint_bridge.cpp under ATV_HAS_LIBCINT, keeping the default ctor
- [x] Update production call sites: src/cpu/scf/internal/hamiltonian_builder.cpp:22-23,64,66,73,127,129,141-142,158,169 and src/cpu/scf/integral_derivatives.cpp:689,1001 to member/constructor spellings only
- [x] Add doxygen (@brief/@param with units/@return) on every new member method and constructor per doc.style
- [x] Verify verbatim kernel move: git diff review of the six operator kernels + two pair-list builders + four integral builders + libcint bridge shows zero token changes inside loop/contraction bodies beyond signature, `*this`-alias line, indentation, and doc comments
- [x] Extend tests/architecture/no_factory_scrub.sh: 13-name tree-wide blacklist, cpu/scf-scoped `build_[a-z_0-9]*[[:space:]]*\(` wildcard with member-access/comment filters and the commented exemption whitelist (grid/minao/basis_set/shell_to_atom + linalg density.h scoping note); verify PASS on the completed tree and FAIL on an injected old-name probe
- [x] Run full check + test suite (`ATV_BUILD_HF=ON`, libcint 配置若可用则一并), recording test_integrals / test_rys_eri wall time against a pre-change baseline on the same node

## Testing strategy

- **快乐路径**：现有 `test_integrals`（AO 级 vs 块级 G parity，H2O/CH4/H2，1e-9）与 `test_rys_eri`（rys vs libcint 参考）改成员拼写后原断言全绿——这就是行为不变的主证据，无需新 fixture。
- **边界**：`ATV_HAS_LIBCINT` 关/开两种配置均编译通过（`LibcintContext` 构造函数在宏内）；`IntegralPack` 空 shells 分支（hamiltonian_builder.cpp:21 的 guard）行为不变；`compute_eri=true/false` 两路径不变。
- **门禁自证**：完成态跑 `no_factory_scrub.sh` 必须 PASS；向 `src/` 注入含 `build_g_rhf_direct(` 的临时文件必须 FAIL（随后删除探针）。
- **领域验证（science.required）**：数值逐位不变由 parity 单测 + 全套 SCF 单测（`ctest`）承载；性能量级不变由 test_integrals / test_rys_eri 前后计时对照承载（单遍 ERI 融合 + Schwarz 筛选逐字搬移的可观测面）。
- 单测均为秒级单函数测试，不新增任何 e2e（regressions/ 不动：无新 feature，纯拼写迁移）。

## Out of scope

- `linalg::build_density_matrix` 改名/收编（GPW 轮，含 gpw_solver 六调用点）；`linalg::build_orthogonalization_matrix` 同理。
- gpw/pw/semiempirical/md/neighbor/comm/cuda 树的一切 `build_*`（各自轮次；本轮门禁通配不覆盖，旧名黑名单覆盖全树）。
- cpu/scf 域内 `build_molecular_grid`、`build_minao_*`、`build_shell_list_from_table`、`build_shell_to_atom_from_integrals`（带注释豁免，随 IntegralPack 删除轮 / gradient 轮收编）。
- 任何 kernel 算法、并行策略、筛选阈值改动；`nuclear_repulsion_energy`（非 build_* 名）；`blocks_from_shells`；IntegralPack 本身的删除（后续债务轮）。

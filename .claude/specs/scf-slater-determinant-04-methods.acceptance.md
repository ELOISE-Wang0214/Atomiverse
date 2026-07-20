---
slug: scf-slater-determinant-04-methods
criteria:
  - id: ac-001
    summary: cpu/scf 头文件中 13 个 build_* 自由函数声明清零，成员方法就位
    type: code
    pass_when: |
      grep -rE 'build_(g_rhf_direct|g_rhf_direct_block|jk_uhf_direct|j_direct|k_direct_block|coulomb_exchange_fused|shell_pair_list|shell_block_pair_list|overlap|kinetic|nuclear_attraction|eri_dense|libcint_context)[[:space:]]*\(' include/atomiverse/cpu/scf src/cpu/scf 返回零命中（注释行除外）；internal/shell_pair.h 声明两构造函数 + g_rhf/jk_uhf/j 与 g_rhf/k 成员；internal/integrals.h 声明 coulomb_exchange_fused 公有成员与 overlap/kinetic/nuclear_attraction/eri_dense 私有成员；libcint_bridge.h 声明 LibcintContext(const Frame&)。
    status: pending
  - id: ac-002
    summary: operators.h 变 #error 毒丸，internal/two_electron.h 删除
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/operators.h 内容为与 scf/shell_pair.h 同款的 #error 毒丸 shim；include/atomiverse/cpu/scf/internal/two_electron.h 不存在于树中；全树无任何文件再 include 这两个路径（毒丸自身除外）。
    status: pending
  - id: ac-003
    summary: test_integrals AO-vs-block parity 以成员拼写通过
    type: runtime
    pass_when: |
      ./build/tests/test_integrals 全绿，其中 BlockG_H2O / BlockG_CH4 / BlockG_H2 使用 ShellPairList{pack} / ShellBlockPairList{pack} 构造与 .g_rhf(...) 成员拼写，1e-9 parity 断言不变。
    status: pending
  - id: ac-004
    summary: test_rys_eri 以构造函数拼写通过
    type: runtime
    pass_when: |
      ./build/tests/test_rys_eri 全绿，6 处 build_shell_block_pair_list 调用点全部替换为 ShellBlockPairList{...} 构造拼写。
    status: pending
  - id: ac-005
    summary: 全套单测在 ATV_BUILD_HF=ON 配置下全绿
    type: runtime
    pass_when: |
      cmake 配置含 -DATV_BUILD_HF=ON 时 `ctest --test-dir build --output-on-failure` 退出码 0（libcint 可用的环境下 ATV_HAS_LIBCINT 路径同样编译并通过）。
    status: pending
  - id: ac-006
    summary: kernel 逐字搬移——循环体零 token 改动
    type: code
    pass_when: |
      对 13 个被迁函数逐一做 git diff 审查：函数体内（不含签名行、`*this` 别名/`const int n = n_basis;` 引导行、`*this = std::move(out);` 收尾行、缩进与 doc 注释）零 token 改动；未引入任何表达式模板、算符代数或新的中间物化；OpenMP 归约、canonical 序、8 重对称散射、Schwarz 阈值逻辑逐字保留。
    status: pending
  - id: ac-007
    summary: no_factory_scrub.sh 在完成态全绿
    type: runtime
    pass_when: |
      `tests/architecture/no_factory_scrub.sh <repo-root>` 退出码 0 且输出 PASS；脚本含 13 个旧名黑名单条目、cpu/scf 域通配 pattern、成员/限定访问过滤，以及每项带注释说明的豁免白名单（grid/minao/basis_set/shell_to_atom + linalg density.h 作用域说明）。
    status: pending
  - id: ac-008
    summary: 门禁可证伪——旧名探针使脚本失败
    type: runtime
    pass_when: |
      在 src/ 下临时写入一行 `build_g_rhf_direct(`（非注释）的探针文件后运行 scrub 脚本退出码非零并指名该 pattern；删除探针后恢复 PASS。
    status: pending
  - id: ac-009
    summary: 性能量级不变——SCF 单测耗时对照
    type: performance
    pass_when: |
      同一节点、同一 Release 配置下，改动后 ./build/tests/test_integrals 与 ./build/tests/test_rys_eri 的 wall time 均 ≤ 改动前基线的 1.5 倍（基线在动手前测得并记录于 PR/commit 说明）。
    status: pending
  - id: ac-010
    summary: 新成员方法 doxygen 完整（含单位）
    type: docs
    pass_when: |
      每个新成员方法与构造函数带 @brief、@param（矩阵维度与 Hartree/Bohr 单位注明）、@return（如适用）；internal/one_electron.h 与 integrals.cpp 顶部示例代码更新为成员拼写，无残留 build_* 旧拼写示例。
    status: pending
---

# Acceptance criteria

- ac-001/ac-002 锁定 API 形态迁移完成且旧表面物理移除（毒丸 + 删除文件）。
- ac-003/ac-004/ac-005 是行为不变的运行时主证据（RED-before-GREEN 的 GREEN 面）。
- ac-006 是"逐字搬移"的绑定审查条款——引导行/收尾行是唯一允许的新增行。
- ac-007/ac-008 保证门禁既在完成态开闸全绿、又真的咬人（可证伪）。
- ac-009 承载"零性能回归"的量级判据；1.5× 阈值吸收节点噪声，量级漂移必被捕获。
- ac-010 对应 doc.style=doxygen 的强制项。

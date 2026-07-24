---
slug: scf-review-03-posthf-dedup
title: "Post-HF dedup: hoist EOM-CCSD intermediates + determinant-CI engine"
status: done
created: 2026-07-24
revised: 2026-07-24
domain: scf / post-hf
queue: scf-review
depends_on: [scf-review-01-linalg, scf-review-02-eri-hoists]
grilled: true
---

# Post-HF dedup: hoist EOM-CCSD intermediates + determinant-CI engine

## Summary

本 spec 消除 `cpu/scf` 后 HF（post-HF）代码中两处经架构评审确认的重复实现，且保证数值不变。第一处：`make_tau` 在 `ccsd.cpp`、`eom_ee_ccsd.cpp`、`eom_ea_ccsd.cpp`、`eom_ip_ccsd.cpp` 中共有 4 份逐字节相同的拷贝，`EomImds` + `build_imds` 在 EE 与 EA 之间 318 行逐字节相同、在 IP 中共享其中约 85 行开头；EE/CC2/SF 已正确共享一份 `build_imds`，而 EA/IP 各自分叉。将其提升到新的私有头 `internal/eom_ccsd_common.h`（配对 `.cpp`）。第二处：行列式 CI 机器（`choose`/`enumerate_strings`/`occ_list`/`jw_single_sign`/Slater–Condon 矩阵元/`ci_element`/稠密 CI 求解/`make_rdm12`）在 `fci.cpp` 与 `casscf.cpp` 中各实现一遍，而 FCI、CASCI、CASSCF、NEVPT2 四个方法本应由同一引擎驱动。将其提升到新的私有头 `internal/ci_engine.h`（配对 `.cpp`）。用户可见行为、公共 API 与所有已交付方法的能量数值均不改变。

## Domain basis

两处提升都涉及物理内核，正确性以解析值钉死，故列出定义与引用。

- EOM-CCSD 中间量 τ（Gauss–Stanton 表 III，PySCF `gintermediates.make_tau`）：

  τ_ij^ab = t_ij^ab + ½·fac·P(ij)P(ab) t_i^a t_j^b

  其中 P(pq) 为反对称化算符。`fac=1` 给出 τ，`fac=0.5` 给出 τ̃。此定义在 4 份拷贝中逐字节一致，提升后为唯一实现。

- 行列式基 Slater–Condon 规则（Szabo & Ostlund §2.3）给出 ⟨I|H|J⟩：
  - 对角（同一行列式）：Σ_i h_ii + Σ_{i<j}(⟨ii|jj⟩ − ⟨ij|ji⟩) + 核排斥。
  - 单激发 p→q：带 Jordan–Wigner 相位的 h_pq + Σ_j(⟨pq|jj⟩ − 交换项)。
  - 双激发：同自旋 (pr|qt)−(pt|qr)、异自旋仅直接库仑 (pa,ra|qb,sb)，相位为两次顺序单激发相位之积。

  规范实现取 `fci.cpp` 中经事后重构的 `single_element`/`same_spin_double_element`/`opposite_spin_double_element`（`fci.cpp:119-170`）；`casscf.cpp:97-179` 保留的是重构前分叉。由于两分叉的索引约定不同（`h(p,q)` vs `h(q,p)`，`(p,q,..)` vs `(q,p,..)`），仅在积分对称性成立时数值相等——因此本 spec 强制加入数值等价门与手算 2 电子矩阵元钉死单测。

参考文献：
- Stanton, Gauss, Watts, Bartlett, J. Chem. Phys. 94, 4334 (1991). DOI 10.1063/1.460284
- Gauss & Stanton, J. Chem. Phys. 103, 3561 (1995). DOI 10.1063/1.470372
- Szabo & Ostlund, *Modern Quantum Chemistry*, Dover (1996), §2.3, ISBN 0-486-69186-1

## Design

两处提升相互独立，落在同一 `cpu/scf` 域内，无新公共 API、无新顶层概念。

**提升一 — EOM-CCSD 公共中间量。** 新增私有头 `internal/eom_ccsd_common.h`（声明）+ `src/cpu/scf/internal/eom_ccsd_common.cpp`（定义），置于 `internal/ccsd_pack.h` 旁。内容为规范的 `make_tau` 与完整版 `EomImds`（含 `Wvvvv`/`Wvvvo`/`Wvovv`）。**Audit 修订（grill 2026-07-24）：`build_imds` 不得以自由函数形式存续**——它是 `no_factory_scrub` 在 dev 上现存的红灯违规之一（eom_ee:70 / eom_ea:66 / eom_ip:67 三处均被门禁点名）；提升时改为 **`EomImds` 的计算构造函数 `explicit EomImds(const CCSDSolution&)`**（house compute-in-ctor 先例：`CanonicalOrth(S,tau)`、`SpinResolution`），调用点变为 `const EomImds im(sol);`。此举同时偿清 3 处既有门禁违规（iron law paydown）。`make_tau` 未被门禁点名（纯数学自由函数，无 owner，属许可例外）——保持函数形态。`CCSDSolution`（`internal/ccsd_pack.h:20`）为已存在的共享输入类型，本头仅 `#include` 依赖它、绝不复制。改写点：`ccsd.cpp` 删除本地 `make_tau`；`eom_ee_ccsd.cpp`（含 EE/CC2/SF 三扇区，同一 TU）删除本地三符号；`eom_ea_ccsd.cpp` 删除逐字节相同的三符号；`eom_ip_ccsd.cpp` 删除其子集拷贝。

IP 子集决策：IP 只读取 `Foo/Fvv/Fov/Woooo/Wooov/Wovoo/Wovvo/Woovv`，不读 `Wvvvo/Wvvvv`。评审已确认 IP 的这些共享扇区与完整 `build_imds` 对应部分「逐字节（除注释）相同」，故 IP 改为调用完整 `build_imds`、忽略未用扇区，读取到的值逐位不变（数值保持）。代价：额外分配一份 `Wvvvv`（O(nv⁴)）与其构建开销——与 `CCSDSolution.vvvv`（IP 已持有的反对称 ⟨vv‖vv⟩ 块）同阶；对本路径面向的微小参考体系（n_mo≲14）可忽略。按 iron law 明确记名此代价，不隐藏；惰性/分扇区构建（给 `build_imds` 加 `bool need_vvvv`）是更干净的长期形态，但属超出「行为保持 dedup」的行为/性能改动，列入 Out of scope。实现时须逐点对比 IP 共享扇区与规范 `build_imds`——若发现任何差异即为潜在 bug，须上报（surface），不得静默吸收。

**提升二 — 行列式 CI 引擎。** 新增私有头 `internal/ci_engine.h` + `src/cpu/scf/internal/ci_engine.cpp`，同样置于 `internal/` 下。采用 `casscf.cpp` 的超集（`struct CiSolution`、`dense_ci_solve`、`make_rdm12`），并以 `fci.cpp` 事后重构的 Slater–Condon 元求值器为规范实现：`choose`/`enumerate_strings`/`occ_list`/`jw_single_sign`/`diagonal_element`/`single_element`/`same_spin_double_element`/`opposite_spin_double_element`/`ci_element`。改写点：`fci.cpp` 删除本地组合学与 `dense_ci_ground`，`fci_hamiltonian`/`fci_energy`/`casci_energy` 改用引擎的 `ci_element`/`dense_ci_solve`；`casscf.cpp` 删除本地 CI 分叉，改用引擎（`CiSolution`/`dense_ci_solve`/`make_rdm12`/`ci_element`）。`nevpt2.cpp` 经 `casscf_energy` 间接受益，无需直接改动。若采用规范求值器后 `casscf`/`casci` 数值相对旧分叉发生变化，即证明旧 `casscf` 分叉在某积分对称性不足处有误——须上报并以手算钉死单测判定正确值，不得静默择一。

命名/构造/错误处理沿用现有 post-HF 私有内核惯例：TINY 上限（n_mo≤14、N_det≤2000）throw、无回退（no fallback），与既有 `fci.cpp` 一致。不引入工厂函数、不引入公共 `*Config` 包、不新增公共符号——两个引擎头均在 `scf/internal/`，仅供 `cpu/scf` 内部与内核单测使用。

**Reuse decision**（依据评审所载 librarian 判定）：
- `reuse CCSDSolution`（`internal/ccsd_pack.h`）——EOM 的共享输入类型；`eom_ccsd_common.h` 依赖它，绝不复制。
- `reuse SlaterDeterminant`（`slater_determinant.h`）——单行列式容器；CI 引擎经既有 pack 入口消费它，**不扩展**该公共头。
- `new internal/eom_ccsd_common.h` —— 现无 EOM 公共中间量之家（`eom_*_pack.h` 皆为薄能量入口头）；EA/IP 分叉本应汇入 EE 已共享的一份。
- `new internal/ci_engine.h` —— 现无 CI 引擎之家（`fci_pack`/`casci_pack`/`casscf_pack`/`nevpt2_pack` 皆薄单函数入口头，非引擎家）；以 casscf 超集聚合。

## Files to create or modify

- `include/atomiverse/cpu/scf/internal/eom_ccsd_common.h` (new)
- `src/cpu/scf/internal/eom_ccsd_common.cpp` (new)
- `include/atomiverse/cpu/scf/internal/ci_engine.h` (new)
- `src/cpu/scf/internal/ci_engine.cpp` (new)
- `src/cpu/scf/ccsd.cpp`
- `src/cpu/scf/eom_ee_ccsd.cpp`
- `src/cpu/scf/eom_ea_ccsd.cpp`
- `src/cpu/scf/eom_ip_ccsd.cpp`
- `src/cpu/scf/fci.cpp`
- `src/cpu/scf/casscf.cpp`
- `tests/cpu/scf/test_eom_ccsd_common.cpp` (new)
- `tests/cpu/scf/test_ci_engine.cpp` (new)
- `tests/architecture/posthf_dedup_scrub.sh` (new)
- `regressions/scf/ci-engine-equiv-h2-sto3g/reg_ci-engine-equiv-h2-sto3g.cpp` (new)
- `regressions/scf/ci-engine-equiv-h2-sto3g/REFERENCE.md` (new)
- `CMakeLists.txt`
- `tests/CMakeLists.txt`
- `regressions/CMakeLists.txt`

## Tasks

- [x] Write failing unit tests for CiEngine Slater-Condon elements in tests/cpu/scf/test_ci_engine.cpp (hand-built 2-orbital/2-electron h,g → diagonal/single/double ci_element pinned to 1e-12), register in tests/CMakeLists.txt
- [x] Write failing unit test for make_tau in tests/cpu/scf/test_eom_ccsd_common.cpp (hand-computed τ = t2 + ½P(ij)P(ab)t1t1 to 1e-12), register in tests/CMakeLists.txt
- [x] Hoist make_tau + EomImds into include/atomiverse/cpu/scf/internal/eom_ccsd_common.h + src/cpu/scf/internal/eom_ccsd_common.cpp — build_imds becomes the computing constructor EomImds(const CCSDSolution&) (clears 3 pre-existing no_factory_scrub violations) — and rewire ccsd.cpp, eom_ee_ccsd.cpp, eom_ea_ccsd.cpp, eom_ip_ccsd.cpp (IP uses the full EomImds, ignores Wvvvv/Wvvvo; diff its shared sectors and surface any mismatch)
- [x] Hoist the CI engine (canonical fci Slater-Condon evaluators + casscf CiSolution/dense_ci_solve/make_rdm12) into include/atomiverse/cpu/scf/internal/ci_engine.h + src/cpu/scf/internal/ci_engine.cpp and rewire fci.cpp + casscf.cpp
- [x] Register src/cpu/scf/internal/eom_ccsd_common.cpp + src/cpu/scf/internal/ci_engine.cpp in CMakeLists.txt, and add grep-gate tests/architecture/posthf_dedup_scrub.sh registered in tests/CMakeLists.txt
- [x] Add regression example regressions/scf/ci-engine-equiv-h2-sto3g/reg_ci-engine-equiv-h2-sto3g.cpp + REFERENCE.md (public API RHF::fci / RHF::casci; hard-coded PySCF literal; no live oracle), register in regressions/CMakeLists.txt
- [x] Run full check + test suite

## Testing strategy

单元测试（仅 `tests/`，镜像 `src/cpu/scf/` → `tests/cpu/scf/`，每测单一函数）：

- `tests/cpu/scf/test_ci_engine.cpp`（`TestCiEngine`，`#include internal/ci_engine.h`，内核单测允许用 `scf/internal/`）：
  - happy：手构 2 轨道/2 电子的 `h`、`g` 张量，`diagonal_element` 命中解析值（硬编码字面量，1e-12）。
  - edge：`single_element` 单激发含 Jordan–Wigner 相位符号翻转的一例（1e-12）。
  - domain validation：`ci_element` 分派到异自旋双激发，命中手算 (pa,ra|qb,sb)（硬编码，1e-12）——此测正是「哪一分叉正确」的判定，若旧 casscf 分叉与之相悖则暴露为 bug。
- `tests/cpu/scf/test_eom_ccsd_common.cpp`（`TestEomCcsdCommon`，`#include internal/eom_ccsd_common.h`）：
  - happy + edge：`make_tau` 在 1-occ/1-vir 与 2×2 手算 τ 上逐元素命中（硬编码，1e-12），覆盖 P(ij)P(ab) 反对称四项。

架构 grep 门（`tests/architecture/posthf_dedup_scrub.sh`，镜像既有 `*_scrub.sh`，`add_test` 注册为 `test_posthf_dedup_scrub`）：断言 `make_tau`/`EomImds`/`build_imds` 的定义仅出现在 `internal/eom_ccsd_common.{h,cpp}`；`choose`/`enumerate_strings`/`occ_list`/`jw_single_sign`/四个 Slater–Condon 元求值器/`ci_element`/`dense_ci_solve`/`make_rdm12`/`CiSolution` 的定义仅出现在 `internal/ci_engine.{h,cpp}`；`src/cpu/scf/*.cpp` 内无本地再定义（排除注释行与引擎文件自身）。

回归例（`regressions/`，库外消费者，公共头，≤60s，硬编码字面量、无实时第三方）：`regressions/scf/ci-engine-equiv-h2-sto3g/`——用公共 API 计算 H₂/STO-3G 的 `RHF::fci` 与 `RHF::casci(2,2)`（二者张成同一满空间，应恒等），断言两者都等于 `REFERENCE.md` 中记录的 PySCF 2.14.0 满 CI 钉值 `-1.137283834488502 Ha`（工具/版本/命令/日期以注释记名）到 1e-6，且 `|E_fci − E_casci| < 1e-12`——此即统一 CI 引擎「合并无 bug 且两入口一致」的数值等价证明。

既有回归为不变门（本 spec 不改其字面量）：`eom-ee-ccsd`、`eom-ea-ccsd`、`eom-ip-ccsd`、`eom-cc2-h2`、`eom-sf-ccsd-o2`、`ccsd-h2-sto3g`、`ccsd-t-h2o-sto3g`、`cc2-h2-sto3g`、`ccsd-grad-h2o-sto3g`、`casci-h2-sto3g`、`casscf-h2-sto3g`、`nevpt2-h2-sto3g`。

## Out of scope

- EOM-IP 的惰性/分扇区 `build_imds`（避免分配未用的 `Wvvvo`/`Wvvvv`）——属性能改动，非行为保持 dedup，另开 spec。
- 统一 CI 引擎索引约定后对 `casscf`/`fci` 元求值器做进一步公式化简或性能优化。
- 触碰或扩展公共 `slater_determinant.h`（评审判定：单行列式容器，不得扩展）。
- `nevpt2.cpp` 的直接改动（仅经 `casscf_energy` 间接受益，本身无 CI 重复）。
- 链上 01/02/08 的职责：`cpu/linalg` 的 `eig_nonsymmetric`、`internal/ao2mo.h`、`internal/spin_orbital_eri.h` 与 CCSD 自旋轨道重写——均在本 spec 之外；本处提升的中间量须保持可被 08 的自旋轨道 oracle 调用。
- 已考虑但不采纳的替代：将两符号族做成 header-only `inline`；因 `build_imds` 达约 318 行、被 3 个 TU 包含会致头膨胀，改用与 `internal/ccsd_gradient.{h,cpp}` 一致的 decl+impl 分离。

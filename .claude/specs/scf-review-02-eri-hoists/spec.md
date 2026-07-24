---
slug: scf-review-02-eri-hoists
title: "scf-review-02: hoist ao2mo to internal header + unify spin-orbital ERI helpers"
status: done
created: 2026-07-24
domain: scf / post-hf
queue: scf-review
depends_on: [scf-review-01-linalg]
grilled: true
---

# scf-review-02: hoist ao2mo to internal header + unify spin-orbital ERI helpers

## Summary

本 spec 是审查链条第 2 环，纯行为保持（behavior-preserving）重构，**不改变任何数值**。它做三件事：(1) 把 `ao2mo_ovov` / `ao2mo_full` / `mo_hcore` 三个 AO→MO 变换声明从**公有头** `include/atomiverse/cpu/scf/mp2.h` 迁到**新的内部头** `internal/ao2mo.h`，切断 11 个后 HF 源文件为了拿变换而穿透包含公有 `mp2.h` 的分层泄漏；(2) 把在四个文件里成对重复的自旋轨道反对称化 ERI 助手（`ccsd.cpp`/`cisd.cpp` 的 `struct SO`+`so_eri`，`mp3.cpp`/`mp4.cpp` 的 `spat_`/`spin_`/`phys_`/`antisym_`）提升为**单一权威实现**放入新的 `internal/spin_orbital_eri.h`；(3) 把复制粘贴 ≥4 处的 MP 能量分母 aufbau 守卫收敛为一个共享内联助手 `mp_denominator`。用户可见 API 面（`RHF::mp2(frame)`、`MP2Result` 等）不变；交付判据是既有单元测试与既有后 HF 回归全部按引用文献值原样复现。

## Domain basis

被保持的物理恒等式（非新增物理，仅用于论证等价）：由化学家约定的 MO 双电子积分 `(pq|rs)` 得到物理学家约定的反对称化自旋轨道积分

\f[
  \langle pq \Vert rs \rangle = (pr\,|\,qs)\,\delta_{\sigma_p\sigma_r}\delta_{\sigma_q\sigma_s}
                              - (ps\,|\,qr)\,\delta_{\sigma_p\sigma_s}\delta_{\sigma_q\sigma_r}.
\f]

两种被删除的编码逐位计算的正是同一个 \f$\langle pq\Vert rs\rangle\f$（已对读源码核对：`so_eri` 与 `antisym_`/`phys_` 展开式逐项相同）。MP 二阶分母采用 Møller–Plesset aufbau 约定 \f$\Delta_{ij}^{ab}=\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b<0\f$，非负即退化/非 aufbau 参考。参考：Szabo & Ostlund, *Modern Quantum Chemistry*, Dover (1996), §2.4 与 §6.5。无新增方程、无新常数、无新单位。

## Design

**触及的实体（全部在 `cpu/scf/` 单一 domain 内）。** 三个新增/迁移的编译单位无新增 `.cpp`（`ao2mo` 定义原地保留在 `mp2.cpp`）：

- **`internal/ao2mo.h`（新）** —— 仅承载 `ao2mo_ovov` / `ao2mo_full` / `mo_hcore` 三个自由函数的声明（连同其现有 Doxygen `@brief`/`@param`/`@return` 原样搬运），命名空间维持 `ATV::cpu::scf`（与既有内部头 `internal/mp2_pack.h` 一致，后者也在 `internal/` 下用普通 `scf` 命名空间声明自由函数）。**定义留在 `mp2.cpp`**（我的取舍：这是行为保持重构，新建 `internal/ao2mo.cpp` 只是无收益的搬运；`mp2.cpp` 已在编译它们，改 include 即可）。选择维持 `scf` 而非 `scf::internal`，是为把 11 个调用点的改动收敛到「仅换 include 行」、不触碰任何调用表达式——分层修复的核心（不再出现在公有 `mp2.h` 包含路径上）由头文件位置本身完成，与命名空间无关。
- **`internal/spin_orbital_eri.h`（新）** —— header-inline 自由数学，命名空间 `ATV::cpu::scf`：
  - `struct SpinOrbital { int spat; int spin; static constexpr SpinOrbital from_packed(int so); };`（`from_packed` = `{so>>1, so&1}`，零成本解码位打包索引）。这是一个两 int 的**标签值类型**，不是参数袋（非 `*Config`）。
  - `inline double antisym_eri(const xt::xtensor<double,4>& g, const SpinOrbital&, ×4)` —— **唯一**的反对称化实现（等于原 `so_eri`/`antisym_` 逐项）。
  - `inline double antisym_eri(const xt::xtensor<double,4>& g, int p,int q,int r,int s)` —— 供位打包生产者（mp3/mp4）用的重载，叶端经 `from_packed` 转调上面唯一实现。
  - `inline double mp_denominator(double e_i,double e_j,double e_a,double e_b)` —— 返回 \f$\Delta<0\f$，`>=0` 抛 `std::invalid_argument`（单一权威守卫；见 Out of scope 对文案的说明）。
  - 这些是 CLAUDE.md 明确许可的「无自然属主的纯数学自由函数」例外（对张量做纯运算，风格对标 `cpu/linalg/eigensolve.h`），不违反 OOP-by-default 默认形状。
- **关键编码取舍（本 spec 必须决定并论证的开放点）——采用单一权威编码 `SpinOrbital{spat,spin}`，四个文件全部迁入。** 读两处调用模式后的依据：(a) 两种编码逐项算同一个 \f$\langle pq\Vert rs\rangle\f$，叶端做的都是同样的自旋比较 + 同样的 `g()` gather，struct 与位打包在热循环里**无性能差异**，故不选「更便宜」而选**改动最小 + 语义最清晰**；(b) `ccsd.cpp`/`cisd.cpp` 的 `occ/vir` 是**块自旋序**（先全部 α-occ 再全部 β-occ），`mp3.cpp`/`mp4.cpp` 的 `o_so/v_so` 是**交错序** `so=2·spat+spin`——这两种自旋轨道**枚举顺序是承重的**（决定振幅张量索引布局），把 ccsd/cisd 的循环变量强改成位打包整数会牵动全部下游张量索引，属于数值风险而非纯提升，**超出行为保持范围**；(c) 因此权威编码取 `SpinOrbital{spat,spin}`：ccsd/cisd 把本地 `struct SO` 重命名为 `SpinOrbital`（`.p→.spat`、`.s→.spin`）保持块自旋枚举不动；mp3/mp4 保留位打包循环索引、只在叶端 `antisym_eri(g,p,q,r,s)` 里 `from_packed` 解码保持交错枚举不动。**四文件全部落到同一份反对称化实现**，重复被完全消除，且每个文件的自旋轨道顺序原封不动 —— 真正数值等价。

**Reuse decision**（解析每个 librarian 候选）：
- `reuse ao2mo_ovov` / `reuse ao2mo_full` / `reuse mo_hcore` —— `mp2.cpp` 中的现有定义原样复用，只搬声明，绝不重实现。
- `generalize SO`+`so_eri`（ccsd/cisd）与 `spat_`/`spin_`/`phys_`/`antisym_`（mp3/mp4） —— 四份本地副本提升为唯一 `SpinOrbital`+`antisym_eri`（spawns 移植任务 5、6，而非并行「实现」任务）。
- `generalize` aufbau 守卫（≥4 份 `denom>=0.0` 抛出块） —— 提升为唯一 `mp_denominator`。
- `new — internal/ao2mo.h`：librarian 明确 `internal/mp2_pack.h` 是 **MP2 能量入口**（`mp2_energy`/`ump2_energy`/`romp2_energy`），**不得**把 ao2mo 折进去；通用后 HF 变换机件应有独立内部头。
- `new — internal/spin_orbital_eri.h`：无现存头贴近——`scf_types` 拥有的是 `SpinResolution`（不同概念）、公有 `scf_types.h` 是墓碑、`CCSDSolution` 存的是预物化块（不同表示）、`mp3_pack`/`mp4_pack` 是薄方法入口。命名与错误处理沿用最近模式（`internal/mp2_pack.h` 的普通 `scf` 命名空间自由函数 + `eigensolve.h` 的 header-inline 自由数学）。

## Files to create or modify

- `include/atomiverse/cpu/scf/internal/ao2mo.h` (new) — 三个 ao2mo/mo_hcore 声明 + Doxygen
- `include/atomiverse/cpu/scf/internal/spin_orbital_eri.h` (new) — `SpinOrbital`、`antisym_eri`（两重载）、`mp_denominator`
- `include/atomiverse/cpu/scf/mp2.h` — 删除三个 ao2mo/mo_hcore 声明，保留 `MP2Result`
- `src/cpu/scf/mp2.cpp` — 改 include 为 `internal/ao2mo.h`；aufbau 守卫改用 `mp_denominator`
- `src/cpu/scf/ccsd.cpp` — 删本地 `SO`/`so_eri`，改用 `SpinOrbital`/`antisym_eri`/`mp_denominator`；include `internal/ao2mo.h`
- `src/cpu/scf/cisd.cpp` — 删本地 `SO`/`so_eri`，改用 `SpinOrbital`/`antisym_eri`；include `internal/ao2mo.h`（其近零守卫见 Out of scope，不折）
- `src/cpu/scf/mp3.cpp` — 删本地 `spat_`/`spin_`/`phys_`/`antisym_`，改用位打包 `antisym_eri` + `mp_denominator`；include `internal/ao2mo.h`
- `src/cpu/scf/mp4.cpp` — 同 mp3.cpp
- `src/cpu/scf/casscf.cpp` — include `internal/ao2mo.h`（如未用 `MP2Result` 则去掉 `mp2.h`）
- `src/cpu/scf/fci.cpp` — include `internal/ao2mo.h`；若确实调用 MP2 入口则保留 `mp2.h`/`mp2_pack.h`
- `src/cpu/scf/g0w0.cpp` — include `internal/ao2mo.h`；按需去掉 `mp2.h`
- `src/cpu/scf/bse.cpp` — 同上
- `src/cpu/scf/adc.cpp` — 同上（其分母 clamp 语义不同，不折，见 Out of scope）
- `src/cpu/scf/gf2.cpp` — 同上
- `src/cpu/scf/internal/mp2_gradient.cpp` — include `internal/ao2mo.h`；按需去掉 `mp2.h`
- `tests/cpu/scf/test_spin_orbital_eri.cpp` (new) — 助手单元测试
- `tests/architecture/no_ao2mo_public_scrub.sh` (new) — grep 门禁
- `tests/CMakeLists.txt` — 注册上面两个新测试

## Tasks

- [x] Write failing unit tests for spin-orbital ERI helpers (tests/cpu/scf/test_spin_orbital_eri.cpp → TestSpinOrbitalEri): struct-vs-packed antisym_eri equivalence on a fixture ERI + mp_denominator aufbau throw
- [x] Create include/atomiverse/cpu/scf/internal/spin_orbital_eri.h with SpinOrbital, antisym_eri (struct + bit-packed overloads), mp_denominator
- [x] Create include/atomiverse/cpu/scf/internal/ao2mo.h with the ao2mo_ovov/ao2mo_full/mo_hcore declarations, and remove those three declarations from include/atomiverse/cpu/scf/mp2.h (keep MP2Result)
- [x] Repoint src/cpu/scf/mp2.cpp to include internal/ao2mo.h and adopt mp_denominator at its aufbau guards
- [x] Port src/cpu/scf/ccsd.cpp and src/cpu/scf/cisd.cpp to hosted SpinOrbital/antisym_eri and internal/ao2mo.h (delete local SO/so_eri)
- [x] Port src/cpu/scf/mp3.cpp and src/cpu/scf/mp4.cpp to hosted bit-packed antisym_eri + mp_denominator and internal/ao2mo.h (delete local spat_/spin_/phys_/antisym_)
- [x] Repoint the remaining post-HF consumers (casscf.cpp, fci.cpp, g0w0.cpp, bse.cpp, adc.cpp, gf2.cpp, internal/mp2_gradient.cpp) to internal/ao2mo.h and drop mp2.h where MP2Result/mp2_energy is unused
- [x] Add grep-gate scrub tests/architecture/no_ao2mo_public_scrub.sh and register both new tests in tests/CMakeLists.txt
- [x] Verify the existing post-HF regression families under regressions/scf/ reproduce their cited values unchanged (mp2/mp3/mp4/ccsd/cisd/adc/gf2/g0w0/bse/casscf/fci)
- [x] Run full check + test suite

## Testing strategy

- **单元测试**（`tests/cpu/scf/test_spin_orbital_eri.cpp` → `TestSpinOrbitalEri`，每测一函数）：
  - `antisym_eri` 的 struct 重载与位打包重载在一个小的硬编码 4-index 化学家 ERI 张量上**逐位相等**（覆盖同自旋、异自旋、交叉自旋各分支，含 `+(pr|qs)`/`-(ps|qr)` 两项都命中与只命中其一）。
  - `SpinOrbital::from_packed` 对若干 `so` 值给出 `{so>>1, so&1}`（硬编码期望）。
  - `mp_denominator` 对硬编码 aufbau 轨道能（如 `e_i=-0.5,e_j=-0.4,e_a=+0.3,e_b=+0.2`）返回负分母（硬编码 `-1.4`，`1e-12` 容差）；对非负分母（`e_a=e_b=-1.0`）抛 `std::invalid_argument`。
  - 无 e2e、无外部包、秒级。
- **架构门禁**（`tests/architecture/no_ao2mo_public_scrub.sh`，仿 `no_factory_scrub.sh`）：`grep -E 'ao2mo_(ovov|full)|mo_hcore' include/atomiverse/cpu/scf/mp2.h` 必须为空；`grep -rE 'struct SO\b|so_eri|inline (int|double) +(spat_|spin_|phys_|antisym_)\b' src/cpu/scf` 必须为空。
- **回归门禁（本 spec 的运行时判据；不新增回归——行为保持重构无新公有行为，绑定门禁是既有后 HF 家族按引用文献值原样复现）**：构建 `regressions` 目标后运行 `regressions/scf/` 下经由被触及代码路径的家族，逐个 PASS：`h2o_mp2_sto3g`、`ump2-o2-sto3g`、`romp2-o2-sto3g`、`frozen-core-mp2-h2o`、`scs-mp2-h2o-sto3g`、`sos-mp2-h2o-sto3g`、`h2o_ri_mp2`、`mp2-grad-h2o-sto3g`、`ump2-grad-o2-sto3g`、`mp3-h2-sto3g`、`mp4-h2-sto3g`、`ccsd-h2-sto3g`、`ccsd-grad-h2o-sto3g`、`ccsd-t-h2o-sto3g`、`cc2-h2-sto3g`、`eom-ip-ccsd`、`eom-ea-ccsd`、`eom-ee-ccsd`、`eom-cc2-h2`、`eom-sf-ccsd-o2`、`cisd-h2co-sto3g`、`adc2-h2co-sto3g`、`adc2x-h2co-sto3g`、`ip-adc2-h2o-sto3g`、`ea-adc2-h2-sto3g`、`gf2-h2-sto3g`、`g0w0-h2-sto3g`、`bse-h2-sto3g`、`casscf-h2-sto3g`、`casci-h2-sto3g`、`nevpt2-h2-sto3g`。每例文献值已由其 `REFERENCE.md` 硬编码（无 live oracle）。

## Out of scope

- **不改任何数值**：不动算法、不动振幅张量布局、不动自旋轨道枚举顺序（ccsd/cisd 块自旋序、mp3/mp4 交错序均保持）。
- **不折入不同语义的守卫**：`cisd.cpp` 用 `std::abs(denom)<1e-12` 抛 `std::runtime_error`、`adc.cpp`/`ucphf_pack.cpp` 用非零 clamp，与 `mp_denominator` 的 `>=0` 抛 `std::invalid_argument` **语义不同**，一律不折入（折入会改变行为）。`mp_denominator` 仅替换四处**逐字相同**的 `denom>=0.0→invalid_argument` 块（`mp2.cpp` 同自旋回路、`mp3.cpp`、`mp4.cpp`、`ccsd.cpp`）。这四处合并为**单一权威抛出文案**——这是**唯一**一处可观察改动，且仅落在退化/非 aufbau 的**错误路径**上（无任何回归或单元测试走该路径；实现时若发现有测试断言旧文案，Iron law 要求保留各自文案，改回按调用点传 tag）。
- **不迁命名空间**：ao2mo/自旋轨道助手维持 `ATV::cpu::scf`（对齐 `internal/mp2_pack.h` 先例），不改为 `scf::internal`，以把调用点改动收敛为「仅换 include」。若架构评审要求 `scf::internal`，另开 spec。
- **不动 `ao2mo` 定义位置**：留在 `mp2.cpp`，不新建 `internal/ao2mo.cpp`。
- **不触教程叙述**：`docs/tutorial/ch07-correlation.md` 仅把 `ao2mo_ovov` 作为实现细节提及，用户路径记的是 `rhf.mp2(frame)`，无需改。
- **不属于本环**：03 的 `internal/eom_ccsd_common.h`/`internal/ci_engine.h`、07 的后 HF 性能、08 的 CCSD 空间轨道重写与 `cc_diis` —— 均在本 spec 的 hoist 之上后续进行。

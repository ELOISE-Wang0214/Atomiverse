---
title: Tiny full-CI reference solver for closed-shell RHF
status: approved
created: 2026-07-04
builds_on: hf-ao2mo-mp2
science_required: true
acceptance: hf-fci-tiny.acceptance.md
---

# 微小体系全组态相互作用（Full CI）参考求解器

## Summary

为 Atomiverse 高斯基 SCF 域交付电子相关路线图的 Milestone 3：一个面向微小体系（H2/STO-3G、He/STO-3G 量级）的行列式基全组态相互作用（full CI）求解器 `fci_energy`，在给定有限基下给出**精确**的闭壳层相关能，用作 MP2 以及未来 CC/CI 方法的黄金参考。它是 `ATV::cpu::scf` 命名空间下的一个自由函数——与 `scf_gradient`、`mp2_energy` 同型的“SCF 结果的纯函数”，**不是** Driver、**不是** Component、不注册进任何 Plan，所有行列式 / Slater–Condon / CI 矩阵机制都封装在新建的 `fci.{h,cpp}` 内部（不新建 `cpu/ci/` 目录，以免违反后端支柱“一次性目录”规则）。本 spec 是 hf-ao2mo-mp2 的**硬链后继**（`builds_on: hf-ao2mo-mp2`）：FCI 的 Slater–Condon 规则跨越全部 MO 轨道对类别，需要一个完整的 (pq|rs) MO 张量与单电子 MO 矩阵 h_pq = Cᵀ·Hcore·C，而 hf-ao2mo-mp2 的 `ao2mo_ovov` 只产出 occ×vir×occ×vir 块——因此本 spec 在 `mp2.h`/`mp2.cpp`（AO→MO 变换族的自然归宿）中新增 `ao2mo_full` 与 `mo_hcore` 两个变换函数。因体系微小，构建 N_det×N_det 稠密 CI 矩阵后直接用 `xt::linalg::eigh` 取最低本征值即可。正确性以固定 PySCF 参考值为准：He/STO-3G 边界严格可推导（N_det=1 ⇒ FCI≡RHF ⇒ E_corr=0，无需 PySCF），H2/STO-3G 数字在实现时用 PySCF 再生并逐位钉死；单元测试层为毫秒级结构不变量；所有 PySCF 平价（含引擎自身 MO 积分喂入 `direct_spin1` 的交叉校验）置于验证 benchmark `benchmarks/qc/fci/validation/`（E2E/外部平价属 benchmark，非 tests/）。配套 `examples/14_fci.cpp` 与在 ch07 相关能章中新增一节 FCI 教学。

## Domain basis

记号：来自 ao2mo 的化学家式空间 MO 积分，实系数正则轨道。h[p,q]=单电子（动能+核吸引）矩阵，对称；(pq|rs)=∫∫φ_p(1)φ_q(1) r₁₂⁻¹ φ_r(2)φ_s(2)，实系数下 8 重对称。核排斥能 E_nuc 加到每个对角元一次（取自 `ints.V_nuc`）。

**1. 行列式表示与计数。** 每个行列式由两个占据位串（每自旋通道一个）表示：α_str, β_str ∈ [0, 2^n_mo)。第 p 位置 1 ⇔ 空间轨道 p 被该自旋占据；popcount(α_str)=n_alpha，popcount(β_str)=n_beta。规范自旋轨道**块序**（1α,2α,…,n_mo α, 1β,…,n_mo β）：整个 α 块先于整个 β 块，块内按空间指标升序——此约定固定所有符号。计数：N_det = C(n_mo, n_alpha)·C(n_mo, n_beta)。H2/STO-3G：n_mo=2、n_alpha=n_beta=1 ⇒ C(2,1)·C(2,1)=4 个行列式。He/STO-3G：n_mo=1、n_alpha=n_beta=1 ⇒ 1 个行列式 ⇒ FCI≡RHF 边界。

**2. Slater–Condon 规则**（化学家记号，自旋分辨；行列式/位串基，**非** CSF）。自旋轨道 i↔(空间 ī, 自旋 σ_i)。h[i,j]=h[ī,j̄]δ(σ_i,σ_j)；(ij|kl)_so=(ī j̄|k̄ l̄)δ(σ_i,σ_j)δ(σ_k,σ_l)。

- **Case A（I=J，对角）**：⟨I|H|I⟩ = E_nuc + Σ_i h[ī,ī] + Σ_{i<j}[(īī|j̄j̄) − (īj̄|j̄ī)δ(σ_i,σ_j)]。库仑项自旋无关；交换项仅同自旋对。
- **Case B（差一个自旋轨道 p→q，σ_p=σ_q 否则为 0）**：(−1)^Γ·{ h[p̄,q̄] + Σ_{j∈occ, j≠p}[(p̄q̄|j̄j̄) − (p̄j̄|j̄q̄)δ(σ_p,σ_j)] }。库仑遍历全部公共占据 j（两自旋），交换仅同自旋 j。
- **Case C（差两个，p,q→r,s）**：(−1)^Γ·{ (p̄r̄|q̄s̄)δ(σ_p,σ_r)δ(σ_q,σ_s) − (p̄s̄|q̄r̄)δ(σ_p,σ_s)δ(σ_q,σ_r) }。同自旋双激发两项都在；异自旋双激发只有直接项 (p̄r̄|q̄s̄) 存活。
- **Case D（差 ≥3）**：0。
- **相位 Γ（Jordan–Wigner 宇称，实现关键）**：单激发 p→q（同一自旋串，占据 b）：Γ=popcount(b & BETWEEN(p,q))，BETWEEN 掩码严格介于 min(p,q) 与 max(p,q) 之间的位；符号=(−1)^Γ。同自旋双激发=两个顺序单激发，第一步后更新 b 再计第二步，符号相乘。异自旋双激发（一个 α + 一个 β）：符号=（α 串单激发符号）×（β 串单激发符号），无交叉项（块序保证移动一个 β 电子绝不跨过任何 α 电子）。

**3. 参考协议**（PySCF 在 agent 环境中不可用——**禁止编造数字**；给出配方，实现时再生并钉死）。配方：`gto.M(atom, basis, cart=True, unit='Angstrom')`；`scf.RHF` `conv_tol=1e-12`；`fci.FCI(mf).kernel()[0]`；把几何+基组+`cart=True` 与每个钉死能量一并记录，精度 1e-8 Ha。He/STO-3G **严格可推导**：N_det=1 ⇒ E_FCI=E_RHF、E_corr=0 到机器精度（直接断言，无需 PySCF）。为把 ao2mo bug 与 FCI bug 隔离：另用 `fci.direct_spin1.kernel(h1e, h2e, norb, nelec, ecore=E_nuc)` 喂入**引擎自身**的 MO 积分，与均场驱动路径对比。几何复用 hf-koopmans/hf-ao2mo-mp2 表（H2 键长 0.7414 Å；He `0 0 0`）。

**固定参考值**（PySCF `cart=True`, SCF `conv_tol=1e-12`；E_corr = E_FCI − E_RHF，单位 Hartree）：

| 体系 | STO-3G |
|---|---|
| He | **0.0（严格为零——N_det=1，FCI≡RHF，机器精度断言，不经 PySCF）** |
| H2 | **⟨TODO：实现时用 PySCF `fci.FCI(mf).kernel()[0]` 再生并钉死 @1e-8 Ha；禁止此刻编造数字⟩** |

RHF 总能量锚点（复用 hf-koopmans 表，Ha）：He/STO-3G −2.8077839575；H2/STO-3G −1.1166843871。

**4. 待断言的不变量**：(i) Hermiticity——对角化前 H_CI 对称到 ~1e-12；(ii) 变分——E_FCI ≤ E_RHF + 1e-10；(iii) 无虚轨道恒等式——He/STO-3G E_FCI==E_RHF 到 1e-10（同时烟测 E_nuc + H_00 接线）；(iv) 相关能符号——E_corr = E_FCI − E_RHF ≤ 0。
🔴 **切勿断言 “E_FCI 恢复 ≥ |E_MP2|”——对这些二电子极小基目标它是可证伪的**。推导（二能级 H2/STO-3G，E_RHF=0，H=[[0,V],[V,Δ]]，V=(1̄2̄|1̄2̄)，gap Δ>0）：E_corr^FCI = ½(Δ−√(Δ²+4V²)) = −2V²/(√(Δ²+4V²)+Δ)；E_corr^MP2 = −V²/Δ（首阶）。因 √(Δ²+4V²)+Δ > 2Δ ⇒ |E_corr^FCI| < |E_corr^MP2|：MP2 过束缚，E_MP2 ≤ E_FCI < 0。**正确的交叉不变量**：断言 E_corr^MP2 ≤ E_corr^FCI < 0（MP2 更负），且/或首阶一致性 E_corr^FCI → E_corr^MP2 当 V/Δ→0。此交叉校验把两个 Milestone 绑在一起、同时护住二者。

**参考文献**：
- Slater, J. C. *Phys. Rev.* **34**, 1293 (1929). DOI 10.1103/PhysRev.34.1293
- Condon, E. U. *Phys. Rev.* **36**, 1121 (1930). DOI 10.1103/PhysRev.36.1121（Slater–Condon 规则）
- Knowles, P. J.; Handy, N. C. *Chem. Phys. Lett.* **111**, 315 (1984). DOI 10.1016/0009-2614(84)85513-X（行列式基 FCI 算法）
- Szabo, A.; Ostlund, N. S. *Modern Quantum Chemistry*, §2.3（Slater–Condon）、§4.3（FCI），H2/STO-3G 表（章节号钉前对物理书复核）。
- Jensen, F. *Introduction to Computational Chemistry*, 3rd ed. (2017), Ch. 4（仅框架，不引用其中任何断言数字；章节号钉前复核）。

## Design

**放置**：`cpu/scf` 域——新建 `include/atomiverse/cpu/scf/fci.h` + `src/cpu/scf/fci.cpp`（头/实现分离，镜像 `gradient.h`/`gradient.cpp`、`mp2.h`/`mp2.cpp`）。按仓库规则 “SCF 是内层求解器”，FCI 与 `scf_gradient`/`mp2_energy` 同类，是后 SCF 分析纯函数——**不是** Driver、**不是** Component、不注册进任何 Plan。所有行列式枚举、位串、Jordan–Wigner 相位、Slater–Condon 与 CI 矩阵机制都留在 `fci.{h,cpp}` 内部，不新增 `cpu/ci/` 目录（一次性多文件簇门槛未达，且会违反后端支柱规则）。

**新公共符号（FCI，在 `ATV::cpu::scf`）**：

- `struct FCIResult { double e_fci; double e_corr; };` — 单位 Hartree；`e_corr = e_fci − scf.energy`。
- `FCIResult fci_energy(const IntegralPack& ints, const SCFResult& scf)` — 从 `scf.C_alpha`、`scf.n_alpha` 取输入；**`n_mo = C_alpha.shape(1)`**（不是 `ints.n_basis`——线性相关剔除后 n_mo < n_ao，与 hf-ao2mo-mp2 spec line 50 同语义、对 hf-robust-02 的软兼容）。消费 `ints.Hcore`、`ints.ERI`（仅当 `build_integrals(..., compute_eri=true)` 时物化）、`ints.V_nuc`。流程：`h_pq = mo_hcore(ints.Hcore, C, n_mo)`；`g = ao2mo_full(ints.ERI, C, n_mo)`；枚举两自旋位串生成行列式列表；按 Slater–Condon Case A–D + Jordan–Wigner 相位显式回填稠密 `H_CI`（每对角元含 `ints.V_nuc`）；`auto [w, U] = xt::linalg::eigh(H_CI)`；`e_fci = w(0)`；`e_corr = e_fci − scf.energy`。

**变换族新增（在 `mp2.h`/`mp2.cpp`，与 `ao2mo_ovov` 并列——AO→MO 变换族的自然归宿，保持该族聚在一处；本 spec 因此硬链 hf-ao2mo-mp2）**：

- `xt::xtensor<double, 4> ao2mo_full(const xt::xtensor<double, 4>& ERI, const xt::xtensor<double, 2>& C, int n_mo)` — 完整 (pq|rs) MO 张量，形状 (n_mo, n_mo, n_mo, n_mo)。复用 `ao2mo_ovov` 的四步四分之一变换 BLAS 形状，但索引范围放宽到全 0..n_mo（每步 reshape 成矩阵后用 `xt::linalg`，BLAS-backed，禁止裸标量四重循环）。**这不违反 mp2 spec 的“最窄接口”原则**：FCI 内在需要全部 MO 轨道对类别，故完整张量**就是** FCI 的最窄充分接口。
- `xt::xtensor<double, 2> mo_hcore(const xt::xtensor<double, 2>& Hcore, const xt::xtensor<double, 2>& C, int n_mo)` — 单电子 MO 矩阵 h_pq = Cᵀ·Hcore·C；平凡两次 gemm：`xt::linalg::dot(xt::transpose(C), xt::linalg::dot(Hcore, C))`（此前仓库无生产者，已核实）。全程 `AccReal`=double。

**生命周期/所有权**：无状态自由函数；MO 积分与 H_CI 中间体由函数内局部持有，返回后释放，不缓存。Incore-only。

**守卫（无回退值规则——一律 `std::invalid_argument` 抛出，不返回哨兵；镜像 mp2）**：
- `ints.ERI` 未物化（`ERI.size() == 0`，调用方忘记 `build_integrals(..., /*compute_eri=*/true)`）→ 抛出，消息指明需 `compute_eri=true`；
- `!scf.converged` → 抛出；
- 开壳层输入（`scf.n_alpha != scf.n_beta`）→ 抛出（本 spec 仅闭壳层 RHF 参考）；
- **TINY 可解性守卫（承重）**：`n_mo > 14` 或 `N_det > 2000` → 抛出，消息说明这是校验/参考工具而非生产方法。理由：MO 积分 n_mo⁴ 与稠密 H_CI（大小 N_det×N_det、`eigh` 为 O(N_det³)）都会爆炸式增长；封顶保证毫秒到亚秒级、内存 <50 MB。He/STO-3G(1 det)、H2/STO-3G(4 dets) 远在封顶内。

**示例**：`examples/14_fci.cpp`（接续 13_mp2 之后；11/12/13 由 hf-koopmans/hf-delta-scf/hf-ao2mo-mp2 姊妹链提供，若其尚未落地则暂存编号空洞、不阻塞本 spec，与 mp2 spec 同一注记），打印 H2/STO-3G 的 E_RHF、E_FCI、E_corr 三行；在 `examples/CMakeLists.txt` 的 `if(ATV_BUILD_RHF)` 块注册。

**教程**：在 `docs/tutorial/ch07-correlation.md`（hf-ao2mo-mp2 已把该章由状态章改写为 MP2 实战章）**新增一节** “Exact correlation in a finite basis: full CI for tiny systems”（英文正文，保留章号）：行列式展开 → 相关能即激发行列式的混合 → FCI 为何指数级昂贵 → 为何它是参考而非生产方法（Jensen 风格）。钉死 He/STO-3G E_corr=0 边界与（再生后的）H2/STO-3G 数字作为教程回归常量，与测试同源（与 mp2 spec 把 ch07 数字绑到代码的做法一致）。

**Oracle**：`tests/cpu/scf/pyscf_reference_fci.py`（新，仿 `pyscf_reference_koopmans.py` 模式）：`scf.RHF`, `cart=True`, `conv_tol=1e-12`，输出 `fci.FCI(mf).kernel()[0]` 与 `direct_spin1` 路径的 E_FCI/E_corr `constexpr` 块——实现时跑一次再生、把 H2/STO-3G TODO 单元格钉死。

**测试与校验**：单元测试 `tests/cpu/scf/test_fci.cpp`（注册 `test_cpu_scf_fci`）只跑毫秒级**结构不变量**：He 边界（E_corr==0/finite/E_FCI==E_RHF）、Hermiticity、行列式计数、守卫抛出、变分与相关能符号——**不含任何外部引擎平价数字**。所有 PySCF 平价（H2/STO-3G 对 `fci.FCI`、`direct_spin1` 引擎-MO 交叉校验、MP2-vs-FCI 二能级对比、更大微小体系）属**验证 benchmark**：`benchmarks/qc/fci/validation/test_fci_vs_pyscf.cpp`，经 `atv_add_validation(val_fci_vs_pyscf ...)` 注册（`ATV_BUILD_BENCHMARKS`）。E2E/外部平价属 benchmark，非单元测试。

## Files to create or modify

- `include/atomiverse/cpu/scf/fci.h` (new) — `FCIResult`、`fci_energy` 声明 + Doxygen
- `src/cpu/scf/fci.cpp` (new) — 行列式枚举 / Slater–Condon / Jordan–Wigner 相位 / `xt::linalg::eigh` 实现
- `include/atomiverse/cpu/scf/mp2.h` — 新增 `ao2mo_full`、`mo_hcore` 声明 + Doxygen（与 `ao2mo_ovov` 并列）
- `src/cpu/scf/mp2.cpp` — 实现 `ao2mo_full`、`mo_hcore`
- `CMakeLists.txt` — 在 `ATV_BUILD_RHF` 源文件列表（约第 316 行 `src/cpu/scf/gradient.cpp` 旁）加入 `src/cpu/scf/fci.cpp`
- `tests/cpu/scf/test_fci.cpp` (new) — 结构单元不变量（无外部平价数字）
- `tests/CMakeLists.txt` — 注册 `test_cpu_scf_fci`
- `benchmarks/qc/fci/validation/test_fci_vs_pyscf.cpp` (new) — PySCF 平价（fci.FCI + direct_spin1 交叉校验 + MP2-vs-FCI）
- `benchmarks/qc/fci/validation/pyscf_reference_fci.py` (new) — PySCF oracle（FCI + direct_spin1）
- `benchmarks/CMakeLists.txt` — `if(ATV_BUILD_RHF)` 块内 `atv_add_validation(val_fci_vs_pyscf qc/fci/validation/test_fci_vs_pyscf.cpp)`
- `examples/14_fci.cpp` (new) — 打印 E_RHF / E_FCI / E_corr（H2/STO-3G）
- `examples/CMakeLists.txt` — `if(ATV_BUILD_RHF)` 块内注册 `14_fci`
- `docs/tutorial/ch07-correlation.md` — 新增 full-CI 教学一节

## Tasks

- [ ] Write failing structural unit tests for fci_energy in tests/cpu/scf/test_fci.cpp: He/STO-3G e_corr == 0.0 exactly and std::isfinite and E_FCI==E_RHF @1e-10, H_CI Hermiticity @1e-12, N_det == C(n_mo,n_alpha)*C(n_mo,n_beta), guards throw (ERI absent / unconverged / open-shell / over TINY cap), variational e_fci <= e_rhf+1e-10 and e_corr <= 0; register test_cpu_scf_fci in tests/CMakeLists.txt. NO external-engine parity numbers in the unit gate.
- [ ] Write the val_fci_vs_pyscf validation benchmark (benchmarks/qc/fci/validation/test_fci_vs_pyscf.cpp): H2/STO-3G e_fci vs regenerated PySCF fci.FCI @1e-8, the direct_spin1 engine-MO cross-check, the MP2-vs-FCI two-level cross-check (e_corr_mp2 <= e_corr_fci < 0), and any larger tiny systems; register via atv_add_validation in benchmarks/CMakeLists.txt (inside if(ATV_BUILD_RHF)). All PySCF parity lives here — not tests/.
- [ ] Write oracle script benchmarks/qc/fci/validation/pyscf_reference_fci.py (RHF cart=True conv_tol=1e-12; emits fci.FCI + direct_spin1 E_FCI/E_corr constexpr block; regenerate + pin the H2/STO-3G number at impl)
- [ ] Implement ao2mo_full + mo_hcore in include/atomiverse/cpu/scf/mp2.h + src/cpu/scf/mp2.cpp (full (pq|rs) MO tensor via four BLAS-shaped quarter transforms widened to 0..n_mo; h_pq = Cᵀ·Hcore·C via two gemms; AccReal=double, xt::linalg only)
- [ ] Implement FCIResult + fci_energy in include/atomiverse/cpu/scf/fci.h + src/cpu/scf/fci.cpp (determinant bitstrings, Slater-Condon Case A-D H_CI build, Jordan-Wigner phase, xt::linalg::eigh → eig(0); n_mo = C_alpha.shape(1), E_nuc on each diagonal, invalid_argument guards + TINY cap n_mo<=14 / N_det<=2000); add src/cpu/scf/fci.cpp to the ATV_BUILD_RHF source list in CMakeLists.txt
- [ ] Add docstring per doxygen style with units (Hartree) for every public symbol in fci.h (FCIResult + fields, fci_energy) and the two new mp2.h symbols (ao2mo_full, mo_hcore): @brief / @param / @return
- [ ] Implement examples/14_fci.cpp printing E_RHF, E_FCI, E_corr for H2/STO-3G; register in examples/CMakeLists.txt under ATV_BUILD_RHF
- [ ] Write the full-CI section in docs/tutorial/ch07-correlation.md ("Exact correlation in a finite basis: full CI for tiny systems"; determinant expansion, correlation as excited-determinant mixing, exponential cost, reference-not-production; embeds the He/STO-3G E_corr=0 boundary and the pinned H2/STO-3G number, Slater/Condon + Knowles-Handy + Szabo-Ostlund refs)
- [ ] Verify invariants + tutorial regression + gate: assert E_FCI <= E_RHF + 1e-10 (variational) and e_corr <= 0; assert ch07-quoted constants equal the test constants; measure default test_cpu_scf_fci stays ms-scale and default ctest keeps the ~16 s budget; run the benchmarks/ validation binaries once at acceptance
- [ ] Run full check + test suite

## Testing strategy

- **主路径（默认层）**：H2/STO-3G — `build_integrals(frame, BasisSet::STO3G, 0, /*compute_eri=*/true)` + `rhf_solve(ints, 1e-12, 100)` + `fci_energy`，E_FCI 对再生的 PySCF 固定值容差 1e-8 Ha；`e_corr == e_fci − scf.energy`。
- **边界情形（默认层）**：He/STO-3G 无虚轨道（N_det=1）→ E_corr 严格 == 0.0 且 `std::isfinite`、E_FCI==E_RHF @1e-10（无虚轨道恒等式，兼作 E_nuc + H_00 接线烟测）；ERI 未物化 / SCF 未收敛 / 开壳层输入 / 超 TINY 封顶 → `std::invalid_argument`。
- **单元不变量（默认层）**：H_CI 对角化前对称 ≤1e-12；N_det == C(n_mo,n_alpha)·C(n_mo,n_beta)；变分 E_FCI ≤ E_RHF + 1e-10；相关能符号 E_corr ≤ 0。
- **交叉校验（默认层）**：二能级 H2/STO-3G 上 E_corr^MP2 ≤ E_corr^FCI < 0（MP2 过束缚，把 hf-ao2mo-mp2 与本 Milestone 绑定；🔴 绝不断言 |E_FCI| ≥ |E_MP2|）。
- **领域验证（benchmark，`val_fci_vs_pyscf`）**：H2/STO-3G 对 `fci.FCI` @1e-8；`direct_spin1` 喂入引擎自身 MO 积分（h1e=mo_hcore、h2e=ao2mo_full、ecore=V_nuc）与均场驱动路径对比 @1e-8 Ha（隔离 ao2mo bug 与 FCI bug）；MP2-vs-FCI 二能级交叉校验；任何更大的微小体系。
- **教程回归**：测试断言 ch07 新节引用的 E_FCI/E_corr 数字（与测试固定常量同源），防止文档与代码漂移。
- **示例烟测**：`14_fci` 在 `ATV_BUILD_RHF=ON` 下编译并运行，输出含 E_RHF、E_FCI、E_corr，且 E_RHF + E_corr 在打印精度上等于 E_FCI。
- **门槛保护**：默认 `ctest` 总时长维持 ~16 s 预算；新增默认层用例为毫秒级（体系微小、封顶保证）。

## Out of scope

- 开壳层 FCI（UHF/ROHF 参考、n_alpha≠n_beta）——本 spec 仅闭壳层 RHF 参考，开壳层输入显式抛出。
- 生产级 FCI / 大体系：无 Davidson/Knowles–Handy sigma-vector 直接 CI、无字符串图/α-β 分块稀疏化——TINY 封顶显式拒绝超限体系；稠密 `eigh` 足够。曾考虑直接 CI，因本里程碑目标是参考正确性而非规模，明确排除。
- 截断 CI（CISD/CISDT 等）、活性空间（CASSCF/CASCI）、冻芯——本 spec 是全电子全 CI。
- FCI 解析梯度、密度矩阵、自然轨道、激发态（仅取最低本征值）。
- CUDA 后端 FCI（`cuda::scf` 无对应；CUDA SCFEvaluator 本身无梯度，平价路线未建立）。
- 任何 Driver/Component 包装（如 “FCIEvaluator”）——本 spec 明确采用 `scf_gradient`/`mp2_energy` 式纯函数；组件化留待真实需求出现。
- 重新钉死 hf-ao2mo-mp2 的 MP2 参考表——本 spec 只新增 `ao2mo_full`/`mo_hcore`，不改 `ao2mo_ovov`/`mp2_energy` 的既有接口或数值。

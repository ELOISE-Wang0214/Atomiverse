---
title: SCF Slater 行列式 · 03 SCFResult 消亡（方法持有波函数，消费面全迁）
slug: scf-slater-determinant-03-unbag
status: approved
chain: scf-slater-determinant
link: 03/04
depends: scf-slater-determinant-01-type, scf-slater-determinant-02-rehome
external-depends: scf-driver-refactor-06-criterion-unify
created: 2026-07-20
---

# SCF Slater 行列式 · 03 SCFResult 消亡

## Summary

`SCFResult` 袋消亡。前提已备：01-type 交付 `SlaterDeterminant`（双通道 C/P/本征值/占据 + `homo()`/`lumo()`），02-rehome 已把枚举、`SCFInitialGuess`、`GradientResult`、`IterationRecord` 疏散出 `scf_result.h`（只剩 `SCFResult` 一个 struct），外部链 scf-driver-refactor-06 已让 loop 节点拥有 `converged()`/`stop_reason()`（`cpu_driver.h:910,913`）。本篇把 `MeanFieldMethod`（CRTP，RHF/UHF/RKS/UKS 全家）改为持有 `SlaterDeterminant`（经现 `result_sp_` 的 shared_ptr 槽，`determinant()` 返回 `const&`），`run()` 从 `SCFResult&` 归 `void`（能量经 `frame.meta.get(kEnergy)` [eV] 与 `energy()` [Hartree]，收敛状态经 Driver 的 `converged()`/`stop_reason()`/新增 `iterations()`，波函数经 `determinant()`），`MeanFieldOutputs::result` 改型 `shared_ptr<const SlaterDeterminant>`；`last_result()` 与 `on_evaluate` 回调删除（后者树内零消费者，grep 实证）。内部消费者（gradient / integral_derivatives / mp2 / fci pack）与外部消费者（orbital_info_observer、两个 cube writer、13 个测试文件、1 个回归、8 个示例）全部迁移；`scf_result.h` 删除。嵌套 force-slot（AIMD：MD/Minimize/Evaluate 内 bind 具体方法）路径行为逐位保持。豁免点名不动：`scf::toy::SCFResult`（`solve.h`，ODR 注释）、NDDO/PW/GPW 自有 result 类型、python/（无 SCF 绑定）。

## Domain basis

本篇是重构，物理锚点是**行为保持**——五个 pinned 能量与本征值（`test_method_driver_parity.cpp:154-158`，[Hartree]）逐位不变；回归文献值 E_RHF/STO-3G(H₂O) = −74.961754 Ha ± 1e-3（LibreTexts Live Textbook of Physical Chemistry §28.2，PySCF 2.13 交叉核对，`regressions/scf/h2o_rhf_sto3g/REFERENCE.md`）。

涉及的两处领域内容：

- **自旋污染**：UHF 单行列式 ⟨S²⟩ = S(S+1) + N_β − Σ_{ij}|⟨φ_i^α|S|φ_j^β⟩|²（Amos & Hall, Proc. R. Soc. A 263, 483 (1961), DOI 10.1098/rspa.1961.0175；Szabo & Ostlund, *Modern Quantum Chemistry*, §3.8.4）。现实现 `spin_analysis.h:16` 与之一致（报告值为 ⟨S²⟩ − S(S+1) 偏差，纯自旋态为 0）。守卫不变量：H/STO-3G 单电子纯 doublet 偏差 ≈ 0（`test_method_driver_parity.cpp:211`）；ΔSCF doublet 阳离子 ⟨S²⟩ = 0.75 精确（`delta_scf_reference.h` 不变量 (3)）。
- **Koopmans / HOMO-LUMO**：−ε_HOMO 为冻结轨道电离能上界（Koopmans, Physica 1, 104 (1934), DOI 10.1016/S0031-8914(34)90011-2）；`OrbitalInfoObserver` 的 HOMO/LUMO/gap 输出 [eV] 数值不变，仅计算处从手工扫 `eigenvalues_alpha` 改为 `determinant().homo()/lumo()`。

单位：SCF 内部 Hartree；`frame.meta.get(kEnergy)` 恒为 eV（`units::kHartreeInEv`）；force-slot 写 eV/Å。

## Design

**结果通道四分（裁定汇总）：**

| 原 SCFResult 字段 | 去处 |
|---|---|
| `converged`, `iterations` | Driver：`converged()`/`stop_reason()`（06 已有）+ 新 `iterations()`（见下） |
| `energy` | `frame.meta.get(kEnergy)` [eV]（变更入口恒写）+ `energy()` 访问器 [Hartree] |
| `eigenvalues_*`, `P_*`, `C_*`, `n_alpha/n_beta`, `spin`, `spin_contamination` | `SlaterDeterminant`（`determinant()` / `MeanFieldOutputs`） |
| `residual`, `trace`, 能量分解五项, `integrated_n_*` | **消亡**——全树零读者（grep 实证）；同名 private 成员（`e_one_` 等）留在方法内部 |

**1. 行列式所有权（mean_field_method.h:59-63,167-178,326）。** `result_sp_: shared_ptr<SCFResult>` → `det_sp_: shared_ptr<SlaterDeterminant>`；`determinant()` 返回 `*det_sp_` 的 `const&`。保留 shared_ptr 槽（而非裸值成员）的理由是既有 handle 存活契约：`test_density_cube.cpp:90` 先取 `result_handle()` 再 `add(std::move(scf))`——handle 必须在方法被 move 进 driver 后仍有效。`MeanFieldOutputs::result` 改型 `shared_ptr<const SlaterDeterminant>`；`shells`/`ready` 字段与 aliasing 模式不动；`has_result()` 保留。

**2. `run()` 归 void。** 两个重载均 `void`：`run(Frame&, const SCFInitialGuess* = nullptr)`（`write_energy` 参数删除，**恒**写 `kEnergy` [eV]——Frame-out 原则，嵌套路径本已恒写）；`run(const Frame&, const SCFInitialGuess* = nullptr)`（保留，ΔSCF 需要 const 几何入口，`delta_scf_reference.h:241`；不触 frame）。`drive_scf_`/`drive_scf_core_`（`mean_field_method.cpp:213,290,541,551`）返回型随改 `void`，`DriveOpts::write_energy` 删除。

**3. `energy()` [Hartree] 与 `iterations()` 访问器。** 新增 `double energy() const`（读 `e_total_`）与 `int iterations() const`（epilogue 记录的 executed 计数，`run_nested_erased` 的返回值来源，`mean_field_method.cpp:207`）。理由：(a) const-Frame 入口不写 meta，Hartree 内部量（parity pins、mp2/fci 示例打印、ΔSCF）需一个非袋出口——标量访问器是方法自身 loop 状态（同 `last_grad_` 先例），不是平行世界模型；(b) `CPUDriver` 现无 `iterations()`（06 只交付了 converged/stop_reason），提升为 driver 通用访问器留给 driver 链，不在本篇顺手做。

**4. `on_evaluate` 删除。** `std::function<void(const SCFResult&)>` 回调（mean_field_method.h:189-196,353）全树**零 setter**（grep `\.on_evaluate\(` 无命中）——按最小面拍板删除，不改签名保留。

**5. spin_contamination 归属（拍板 + 理由）。** **存于 `SlaterDeterminant`（缓存字段），由 `drive_scf_core_` epilogue 经现有 `compute_spin_contamination`（`spin_analysis.h`，公式唯一所在，不动）填充。** 理由：⟨S²⟩ 偏差是行列式的内禀性质（谱表述只依赖双通道占据轨道），概念上归行列式；但在非正交 AO 基下求值需要重叠矩阵 S——S 是方法私有工作区（`IntegralPack`），消费者（`test_delta_scf.cpp:69,86,103`、parity doublet 检查）拿不到 S，按需计算的方法签名 `spin_contamination(S)` 对他们不可用。故填充时机在方法 epilogue（现 `mean_field_method.cpp:502-504` 原位），发布值随行列式走。restricted 恒 0（现语义不变）。

**6. 内部消费者签名（散装参数最小化）。** gradient 全链只读行列式数据（P/C/ε/n_α/n_β——`gradient.cpp:512-519` 实证）：`ScfGradient::compute(frame, ints, const SlaterDeterminant&)`（`gradient.h:61`）、`accumulate_xc_gradient`（`gradient.cpp:153`）、`accumulate_two_electron_gradient`（`integral_derivatives.h:89` / `.cpp:667`）同改，无需额外散装参数（method/spin/xc 本已单传）。MP2/FCI 额外要 `converged` 与 `energy`：`mp2_energy(ints, det, e_scf)`、`fci_hamiltonian(ints, det)`、`fci_energy(ints, det, e_scf)`（[Hartree]）；「SCF 未收敛即 throw」的 no-fallback 守卫（`mp2.cpp:120`、`fci.cpp:206`）上移到 `MeanFieldMethod::mp2()/fci()`（查 `this->converged()`），ERI/闭壳层/TINY 上限守卫原位保留。`MP2Result`/`FCIResult` 保留。

**7. 观察者/记录者。** `OrbitalInfoObserver`（`orbital_info_observer.cpp:43-60`）：手工 HOMO/LUMO/gap 扫描删除，改 `det.homo()`/`det.lumo()`（01 契约含占据有效性 throw），eV 转换与 TSV 行为不变。`DensityCubeWriter`（`density_cube_writer.cpp:44-53`，读 P_α/P_β/spin）与 `OrbitalCubeWriter`（`orbital_cube_writer.cpp:46-50`，读 C_α/C_β）经改型后的 `MeanFieldOutputs` 机械迁移，`#include scf_result.h` 删除。

**8. 消费面清册（librarian 穷举 + 本篇补漏）。** 请求清单之外验证新增的消费者：`tests/cpu/scf/{test_scf_solver,test_koopmans,test_spin_resolution,test_uhf_solver,test_uks_solver,test_mp2,test_fci}.cpp`、`tests/cpu/test_evaluate_driver.cpp:70,79`、`examples/{02_rhf.cpp:46, 02_rhf.cu:32（用的是 cpu::scf::RHF）, 03_uhf.cpp:28, 04_uks.cpp:27, 11_koopmans.cpp:39, 12_delta_scf.cpp:40,46}`、`include/atomiverse/cpu/driver/single_point.h:21`（注释）、`include/atomiverse/cpu/scf/methods.h:18`（注释）。豁免确认：`tests/mpi/test_serial_on_replica.cpp:147` 是 AM1（不动）；python 无 SCF 绑定。

**9. 嵌套 force-slot 不变。** `run_nested_erased`（`mean_field_method.cpp:197-208`）只改返回值来源（`result_sp_->iterations` → `iterations()`）；forces/energy 写 Frame 的路径逐位保持，`test_method_driver_wire` + `tests/cpu/driver/test_driver_nesting.cpp` 守卫。

## Files to create or modify

- `include/atomiverse/cpu/scf/mean_field_method.h` — det 槽、`determinant()`/`energy()`/`iterations()`、`MeanFieldOutputs` 改型、`run()` void、删 `last_result()`/`on_evaluate`
- `src/cpu/scf/mean_field_method.cpp` — `drive_scf_`/`drive_scf_core_` 构造 `SlaterDeterminant`、epilogue 填充 + spin_contamination、删 trace/分解字段喂入
- `include/atomiverse/cpu/scf/slater_determinant.h`（01-type 交付；实际路径以 01 落点为准）— 增 spin_contamination 缓存字段
- `include/atomiverse/cpu/scf/scf_result.h` (delete)
- `include/atomiverse/cpu/scf/internal/gradient.h` + `src/cpu/scf/gradient.cpp`
- `include/atomiverse/cpu/scf/internal/integral_derivatives.h` + `src/cpu/scf/integral_derivatives.cpp`
- `include/atomiverse/cpu/scf/internal/mp2_pack.h` + `src/cpu/scf/mp2.cpp`
- `include/atomiverse/cpu/scf/internal/fci_pack.h` + `src/cpu/scf/fci.cpp`
- `include/atomiverse/cpu/scf/orbital_info_observer.h` + `src/cpu/scf/orbital_info_observer.cpp`
- `src/cpu/recorder/density_cube_writer.cpp`、`src/cpu/recorder/orbital_cube_writer.cpp`
- `include/atomiverse/cpu/scf/methods.h`、`include/atomiverse/cpu/driver/single_point.h`（doc 注释）
- 测试：`tests/cpu/scf/{test_method_driver_parity,test_method_driver_wire,test_scf_solver,test_koopmans,test_spin_resolution,test_uhf_solver,test_uks_solver,test_mp2,test_fci,test_orbital_info_observer}.cpp`、`tests/cpu/scf/delta_scf_reference.h`、`tests/cpu/test_evaluate_driver.cpp`、`tests/cpu/recorder/test_density_cube.cpp`
- 回归：`regressions/scf/h2o_rhf_sto3g/reg_h2o_rhf_sto3g.cpp` + `REFERENCE.md`
- 示例：`examples/{02_rhf.cpp,02_rhf.cu,03_uhf.cpp,04_uks.cpp,11_koopmans.cpp,12_delta_scf.cpp,13_mp2.cpp,14_fci.cpp}`
- 文档（docs）：`docs/tutorial/{ch03-open-shell,ch04-basis-sets,ch06-hybrid-dft,ch07-correlation}.md`、`docs/developer/cpu-internals.md`、`docs/api/index.md`、`docs/roadmap/post-hf.md`

## Tasks

- [ ] Write failing tests for the unbagged surface: migrate tests/cpu/scf/test_method_driver_parity.cpp, tests/cpu/scf/test_method_driver_wire.cpp, tests/cpu/test_evaluate_driver.cpp to determinant()/energy()/iterations()/converged() (RED: surface absent; pinned Hartree literals unchanged)
- [ ] Implement SlaterDeterminant ownership in include/atomiverse/cpu/scf/mean_field_method.h + src/cpu/scf/mean_field_method.cpp: det_sp_ slot, determinant()/energy()/iterations(), typed MeanFieldOutputs, void run() both overloads (mutating entry always writes kEnergy [eV]), epilogue fills the determinant, delete last_result()/on_evaluate/DriveOpts::write_energy/trace-and-decomposition fills; doxygen with units
- [ ] Implement the spin_contamination cache on SlaterDeterminant (include/atomiverse/cpu/scf/slater_determinant.h, 01 landing) filled in drive_scf_core_ epilogue via spin_analysis.h compute_spin_contamination (formula home unchanged)
- [ ] Migrate internal packs to const SlaterDeterminant&: internal/gradient.h + gradient.cpp, internal/integral_derivatives.h + integral_derivatives.cpp, internal/mp2_pack.h + mp2.cpp (+e_scf param), internal/fci_pack.h + fci.cpp (+e_scf param); move the converged no-fallback guards into MeanFieldMethod::mp2()/fci()
- [ ] Migrate observers/recorders: orbital_info_observer.h/.cpp to det.homo()/lumo() (delete the manual eigenvalue scan), density_cube_writer.cpp and orbital_cube_writer.cpp to the typed handle (drop scf_result.h includes)
- [ ] Migrate remaining tests: test_scf_solver, test_koopmans, test_spin_resolution, test_uhf_solver, test_uks_solver, test_mp2, test_fci, test_orbital_info_observer, test_density_cube, delta_scf_reference.h (copy determinant scalars/matrices between warm/cold runs)
- [ ] Migrate regressions and examples: reg_h2o_rhf_sto3g.cpp + REFERENCE.md to frame.meta.get(kEnergy)/converged(); examples 02_rhf.cpp, 02_rhf.cu, 03_uhf.cpp, 04_uks.cpp, 11_koopmans.cpp, 12_delta_scf.cpp, 13_mp2.cpp, 14_fci.cpp to energy()/determinant()
- [ ] Delete include/atomiverse/cpu/scf/scf_result.h, reroute stale includes, update methods.h/single_point.h doc comments, and verify the SCFResult scrub (only scf::toy in solve.h/solve.cpp survives)
- [ ] Update the docs (tutorial ch03/ch04/ch06/ch07, developer/cpu-internals.md, api/index.md, roadmap/post-hf.md) off the SCFResult user path — docs
- [ ] Run full check + test suite

## Testing strategy

- **行为保持（单元，秒级）**：五个 pinned 能量/本征值（parity）逐位不变；`RunAndComponentLaunchAgree` 双入口同能量；`test_scf_solver` 双实例/重复 run 稳定性经 `energy()`/`determinant().P_alpha` 重述。
- **通道正确性**：变更入口后 `frame.meta.get(kEnergy)` [eV] 与 `energy()`×kHartreeInEv 一致（wire 测试既有断言）；`converged()`/`stop_reason()=="scf_converged"`（parity ConvergenceWiring）。
- **嵌套隔离/AIMD**：`test_method_driver_wire` + `test_driver_nesting` 绿——force-slot 绑定路径与内层收敛不停外层均不变。
- **handle 存活**：`test_density_cube` 先取 handle 再 move 进 driver，改型后仍绿。
- **Domain validation（science.required）**：(1) 自旋不变量——UHF-H、UKS-He 污染 ≈ 0，ΔSCF doublet ⟨S²⟩ 分析值（`test_delta_scf` 既有容差）；(2) 回归文献门 `reg_h2o_rhf_sto3g` −74.961754 Ha ± 1e-3；(3) MP2/FCI 在 He/STO-3G 零关联极限（parity + test_mp2/test_fci）。
- **债务终结（grep）**：`SCFResult` 在 include/src/tests/examples/regressions 仅剩 `scf::toy`（solve.h/solve.cpp）；`last_result|on_evaluate` 在 Gaussian-SCF 消费面为零；`scf_result.h` 不存在。
- **边界**：mp2()/fci() 对未收敛参考 throw（守卫上移后单测保留）；`run(const Frame&)` 不写 meta、不动坐标。
- 本篇无新 regressions 案例——纯重构，既有 h2o_rhf_sto3g 即行为守卫。

## Out of scope

- `build_*` 清理与 04-methods 的其余方法面工作（链下一环）。
- PW/GPW/NDDO/FCI/MP2 自有 Result 类型清理（`PWResult`/`GpwResult`/`NDDOResult`/`MP2Result`/`FCIResult` 保留，另立 spec）。
- `scf::toy::SCFResult`（solve.h，ODR 注释豁免）与 CUDA SCF 面。
- `CPUDriver` 通用 `iterations()` 访问器（driver 链的事，本篇只加方法本地计数）。
- python/ 绑定（无 SCF 绑定，无事可做）。

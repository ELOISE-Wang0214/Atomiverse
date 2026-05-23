---
title: CUDA MD categorical setter API parity with CPU
status: approved
created: 2026-05-23
parent: msu-pe-validation.md (Phase G.4)
chain: msu-pe-cuda-components-* (link 2/8)
---

# CUDA MD categorical setter API parity with CPU

## 综述 (Summary)

本规范是 `msu-pe-cuda-components-*` 链 8 个子链中的 2/8 链,实现 `msu-pe-validation.md` Phase G.4 的功能补齐:把 `cuda::md::MD`(`include/atomiverse/cuda/md/md.cuh`)从只暴露原始 `add()` / `on(phase).add()` 升级为与 `cpu::md::MD` 对齐的 LAMMPS-shaped 范畴化 setter API —— 11 个 setter 全部就位:`set_pair / set_bond / set_angle / set_dihedral / set_neighbor / set_integrator / set_thermostat / set_barostat / set_observer / set_record / set_thermo`,且配套 CPU 侧 `install_infra()` 在单 GPU 下成立的那一份自动注入逻辑(barostat 在场 → 自动注入 `cuda::VirialStress`,无 halo 交换 / 无 MPI)。原始 `add()` / `on(phase).add()` 路径保持不变,作为既有 escape hatch。这是 03 / 04 / 06 等下游链直接依赖的底层 —— 它们要靠 `set_integrator` / `set_thermostat` / `set_barostat` 把新加的 `NVELimit` / `FixLangevin` / `NoseHooverNPT` 接入 CUDA MD 流水线。

## Design

**作用域。** 修改 `include/atomiverse/cuda/md/md.cuh` 这一个 header,新增 11 个 setter + 一个 `install_infra()` override。不动 `cuda::CUDADriver` —— 它在 `cuda_driver.cuh` 的 doc 里明确"deliberately provides **no** categorical-setter API and no `install_infra()` auto-injection";本链让 `cuda::md::MD` 自己在驱动具体类一层补齐这套表面,不上提到 `CUDADriver`,以保持 CUDA 驱动 taxonomy 链的范围承诺不变。这是与 CPU 不同的位置选择(CPU 是放在 `CPUDriver`),但等价的外部用户体验。

**单数 vs 复数 slot,精确镜像 CPU(`cpu/cpu_driver.h` + `cpu/md/md.h`)。**

- 单数 + 显式 phase / order: `set_pair`(eval_phase = `md::Phase::force` @ 0)、`set_bond / set_angle / set_dihedral`(同 phase @ 0)、`set_neighbor`(同 phase @ -10,排在 evaluator 前)。idempotent 复设同 `name()`、不同 `name()` 抛 `std::invalid_argument`,与 CPU 字符串一致。
- 单数 + default_bindings: `set_integrator` —— 走组件自己的 `default_bindings()`,因为 `cuda::md::VelocityVerlet` / `NVELimit`(链 03)/ `NoseHooverNPT`(链 04)各自绑定到不同的子相位组合。
- 单数 cadence-only: `set_thermo(int every)` —— 绑定 `cuda::Thermo{every}` 到 `md::Phase::end_of_step @ +20`,与 CPU 排序一致。
- 复数 / 追加: `set_thermostat`、`set_barostat`、`set_observer`、`set_record` —— 每次调用直接 `this->add(...)`,保持插入顺序。

**自动注入(单 GPU 子集)。** 重写 `install_infra()`(override 自 `CUDADriver::build()` 调用前的钩子;`cuda::CUDADriver` 当前没有 `install_infra` 虚函数 —— 本链不修改 `CUDADriver`,而是在 `cuda::md::MD::build()` 中显式覆盖 `CUDADriver::build()`:先调用本地 `install_infra_()` 私有方法,再调 `CUDADriver::build()`)。规则只一条,镜像 CPU 的 stress 分支:扫描已注册的 pending 组件,若有任何组件 `role()==ComponentRole::Barostat` 且尚未注册 `name()=="virial_stress"` 的组件,自动 `add(cuda::VirialStress{})`;去重对照名字。CUDA `ComponentBase` 目前没有 `reads()` 虚函数(见 `cuda/component.cuh` 注释 + `cuda/md/barostat.cuh` 注释),所以以 `role()==Barostat` 作为 Stress 消费者的代理信号;`cuda::md::Barostat`(Berendsen)和 `NoseHooverNPT`(链 04)都设 `role()==Barostat`,语义吻合。**halo 自动注入 / MPI 不在 CUDA 单 GPU 范围**,显式不实现(`msu-pe-validation` 固定范围决策 1)。**neighbor-list 缺失检查**也复刻:若有 evaluator(`role()==Evaluator` 且 cutoff > 0)且未 `set_neighbor()`,`build()` 抛清晰错误。

**保持向后兼容。** 原 `add()` / `on(phase).add()` / 构造器里的 Zero* / KineticEnergy 绑定全部保留;新 setter 全部是 additive,不替换。`tests/cuda/md/test_md.cu` 现有 NVEEnergyConservation 等用例(直接走 `add(VelocityVerlet) + on(force).add(NeighborList).add(LJ)`)必须继续 pass,作为零回归证据。

**slot 状态变量。** 新增私有成员 `std::optional<std::string>` 五个:`pair_set_ / bond_set_ / angle_set_ / dihedral_set_ / neighbor_set_ / integrator_set_`(integrator 也是单数),外加 `std::optional<int> thermo_every_`。setter helper(`set_singular_at_` 与 `set_singular_default_`)在 header 里以模板内联实现,与 CPU 同名同语义,不重用 `cpu::CPUDriver` 的实现(避免跨后端模板纠缠)。

## Files to create or modify

- `include/atomiverse/cuda/md/md.cuh` — 新增 11 个 setter、`install_infra_()` 与 `build()` override、slot 成员、setter helper 模板;构造器里现有 Zero* / KineticEnergy 绑定保留。
- `tests/cuda/md/test_md.cu` — 追加范畴化 setter parity 测试(见 Testing strategy);现有用例不动。
- `tests/cuda/md/test_cuda_md_components.cpp` (new) — 链上游已规划的统一 CUDA MD 组件测试文件,本链先建文件并在其中放置 `set_thermostat` / `set_barostat` / `set_observer` / `set_record` / `set_thermo` / install_infra(`VirialStress` 自动注入)的 parity 测试;链 03 / 04 之后会向这个文件继续追加 `NVELimit` / `FixLangevin` / `NoseHooverNPT` 用例。
- `CMakeLists.txt` — 把 `tests/cuda/md/test_cuda_md_components.cpp` 注册到 `gtest_discover_tests` 的 CUDA test 列表,保证 CI 跑得到。

## Tasks

- [ ] Write failing tests for CUDA MD categorical setter parity and VirialStress auto-injection (tests/cuda/md/test_md.cu, tests/cuda/md/test_cuda_md_components.cpp)
- [ ] Implement set_pair / set_bond / set_angle / set_dihedral / set_neighbor with singular-slot semantics in include/atomiverse/cuda/md/md.cuh
- [ ] Implement set_integrator / set_thermostat / set_barostat / set_observer / set_record / set_thermo in include/atomiverse/cuda/md/md.cuh
- [ ] Implement install_infra_() VirialStress auto-injection (role()==Barostat trigger) and build() override in include/atomiverse/cuda/md/md.cuh
- [ ] Implement missing-NeighborList build-time error in install_infra_() mirroring CPU
- [ ] Register tests/cuda/md/test_cuda_md_components.cpp in CMakeLists.txt
- [ ] Verify zero regression against existing tests/cuda/md/test_md.cu raw add() / on().add() paths
- [ ] Run full check + test suite

## Testing strategy

- **Happy path** — `tests/cuda/md/test_cuda_md_components.cpp::CudaMdCategoricalApi.SmokeAllSetters` 构造 `cuda::md::MD sim; sim.set_pair(cuda::LJ{1.0, 1.0, 3.0}); sim.set_neighbor(cuda::NeighborList{3.0}); sim.set_integrator(cuda::md::VelocityVerlet{0.001}); sim.set_thermo(5);` 然后 `sim.build()` + `sim.run(frame, 10)` 在 64 原子 LJ 立方系统上完成,无 NaN、能量为有限实数 —— 证明范畴化 API 自洽。`CudaMdCategoricalApi.SetterParityWithCpu` 用静态断言 + 名字字符串列表确认 `cuda::md::MD` 暴露的 11 个 setter 名与 `cpu::md::MD` 完全一致(通过 SFINAE detection idiom 检查 `decltype(sim.set_xxx(...))` 良构性)。
- **Edge cases** — `CudaMdCategoricalApi.SingularRejectsMismatch` 同一 setter 两次不同 component 名调用,第二次抛 `std::invalid_argument`,异常 message 包含两个 name;`CudaMdCategoricalApi.SingularIdempotentOnSame` 同一 setter 两次同 name,不抛、只绑定一次(用 pending 计数验证)。`CudaMdCategoricalApi.PluralAppends` `set_observer` / `set_record` / `set_thermostat` 每次都追加,插入顺序保留。`CudaMdCategoricalApi.MissingNeighborThrows` set_pair without set_neighbor → `build()` 抛运行时错误,message 命名 evaluator。`CudaMdCategoricalApi.RawAddStillWorks` 重跑 `tests/cuda/md/test_md.cu::MD.NVEEnergyConservation` 路径作为零回归底线。
- **Domain validation** — `CudaMdCategoricalApi.VirialStressAutoInjectedOnBarostat` 设 `sim.set_pair(cuda::LJ{...}); sim.set_neighbor(...); sim.set_integrator(cuda::md::VelocityVerlet{0.001}); sim.set_barostat(cuda::md::Barostat{0.0, 100.0, 0.001}); sim.build();` 之后,内部 pending 列表中应出现 `name()=="virial_stress"` 的组件;`CudaMdCategoricalApi.VirialStressDedupedWhenExplicit` 用户先 `sim.add(cuda::VirialStress{})` 再 `sim.set_barostat(...)`,build 后 pending 中 `virial_stress` 只有一份。`CudaMdCategoricalApi.NoAutoInjectWithoutBarostat` 无 barostat 时不注入 `VirialStress`。

## Out of scope

- 多 GPU / CUDA+MPI 自动注入(halo forward / reverse)—— 由 `msu-pe-validation.md` 固定范围决策 1 显式排除;本链不做。考虑过把 halo 占位接口先放进去,被否决:CUDA 单 GPU 路径根本没有 ghost atom 概念,占位接口只会污染 API。
- 在 `cuda::ComponentBase` 上加 `reads()` 虚函数 —— `cuda/md/barostat.cuh` 注释把这件事"deferred to a later spec",本链遵守该 deferral,自动注入用 `role()==Barostat` 代理。考虑过顺手补 `reads()`,被否决:那会触及所有 CUDA 组件子类,远超 G.4 范围。
- 修改 `cuda::CUDADriver` 让它持有共享 setter API —— 与 `cuda/cuda_driver.cuh` 中"deliberately provides **no** categorical-setter API"的范围承诺冲突。考虑过把 setter 上提以与 CPU 完全镜像,被否决:跨链范围漂移,留给后续 cuda-driver-taxonomy refinement。
- `polyethylene_tension.cpp` 的 CUDA 路径 —— 那是 G.9 / 链 06 的工作,本链只提供 setter 表面。

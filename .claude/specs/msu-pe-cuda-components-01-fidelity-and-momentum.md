---
title: MSU-PE reference-fidelity defects — ZeroMomentum (CPU+CUDA), per-axis Pxx/Pyy/Pzz, energy decomposition
status: in-progress
created: 2026-05-23
revised: 2026-05-23
parent: msu-pe-validation.md (Phase G.1 + G.2 + G.3)
chain: msu-pe-cuda-components-* (link 1/8)
notes: |
  G.3 (energy decomposition) landed: SystemScalars carries
  energy_pair/_bond/_angle/_dihedral, all 4 CPU evaluators write the
  per-term field alongside `energy`, CPU+CUDA ZeroEnergy clears them in
  lockstep (CUDA via one 5-double contiguous memset), MPI reduce
  buffer extended 14→18 doubles.
  Remaining: G.1 (cuda::md::ZeroMomentum component + benchmark
  wire_zero_momentum on stages 2/3/4), G.2 (pressure_per_axis_bar
  + CsvRecorder schema expansion + tools/run_msu_atomiverse.sh G.2 TODO
  drop + canonicalize_msu_output.py passthrough), and the two new
  test files (test_zero_momentum.cu + test_pe_thermo_schema.cpp).
---

# MSU-PE reference-fidelity defects — ZeroMomentum (CPU+CUDA), per-axis Pxx/Pyy/Pzz, energy decomposition

## 综述 (Summary)

本规范是 `msu-pe-cuda-components-*` 链的第 1/8 节,合并伞型路线图 `msu-pe-validation.md` 中三个相互独立、共同决定 `benchmarks/md/polyethylene_tension.cpp` 能否在物理量层面与 `lammps.in` 逐栏对齐的"参考保真缺陷":(1) `fix momentum 1 linear 1 1 1` 等价物在 CPU 上**已存在**为 `cpu::md::ZeroLinearMomentum` 但 CUDA 端缺失且基准未在 stage 2/3/4 接入;(2) 基准 thermo.csv 只输出各向同性压力 `P_bar`,`tools/run_msu_atomiverse.sh` 通过把它复制三遍来伪造 `Pxx/Pyy/Pzz`;(3) 总势能 `E_pot` 是单标量累加,无法分离出 `E_pair/E_bond/E_angle/E_dihedral`。完成后,基准的 thermo.csv 扩展为 `step,stage,T_target_K,T_K,P_bar,Pxx,Pyy,Pzz,E_pot_eV,E_pair_eV,E_bond_eV,E_angle_eV,E_dihedral_eV,E_kin_eV,box_x,box_y,box_z`,`run_msu_atomiverse.sh` / `canonicalize_msu_output.py` 中的 G.2 "stand-in" 注释和 awk 复制都被移除;CUDA `MD` 也可挂 `ZeroMomentum`,为链路下游 `06-cuda-benchmark-path` 提供与 CPU 一致的 stage 2/3/4 配置。

## Domain basis

参考标杆为 `benchmarks/md/reference/polyethylene/lammps.in`(Tschopp 2010,*Polymer*,DOI 10.1016/j.polymer.2010.10.062;CAVS ICME 源 <https://www.cavs.msstate.edu/icme/code/lammps/atomistic_research/polymers/amorphous.php>)。所涉及的三条 LAMMPS 命令的物理定义:

1. **`fix momentum N linear 1 1 1`** — 每 `N` 步从所有原子上扣减质心线动量,`v_i ← v_i - P_cm / M_tot`,其中 `P_cm = Σ_i m_i v_i`,`M_tot = Σ_i m_i`。这是抑制 Langevin/NPT 耦合下"飞行冰块"漂移的标准修正(LAMMPS `Velocity::zero_momentum`,`velocity.cpp:738`)。多 MPI rank 下 `P_cm` 与 `M_tot` 需要先做 `MPI_Allreduce(SUM)` 再扣减,否则每 rank 各自移除局部 COM,会引入横向应力。CPU 实现 `src/cpu/md/init.cpp::ZeroLinearMomentum::launch` 已对齐。
2. **每轴压力 `P_αα = (Σ_i m_i v_iα² + virial_αα) / V`**(LAMMPS `compute pressure`,`compute_pressure.cpp:257`),Voigt 顺序 xx/yy/zz/xy/xz/yz。`frame.scalars.stress`(Σ r·F)与 `frame.scalars.stress_kinetic`(Σ m·v·v)按 `core/system_scalars.h` 协议分别承载维里与动能两部分,各自时间片独立。压力符号:正号为压缩态,匹配 `compute pressure` 与 `NoseHooverNPT` / Berendsen `Barostat`。
3. **能量分解** — 每个 evaluator 在 `src/cpu/potentials/{bond_harmonic,angle_harmonic,dihedral_multi_harmonic,pair_lj_cut}.cpp` 中以 `f.scalars.energy += energy_total` 累加,LAMMPS `thermo_style custom epair ebond eangle edihed` 把同一总和按项拆开。物理上每项已被该 evaluator 内部计算出来,只是没有以独立字段落地。

单位:eV / Å / fs / amu / K(`core/units.h`);压力外部展示用 bar(`kEvPerA3ToBar = 1.602176634e6`);Pxx/Pyy/Pzz 与 P_bar 同单位。

## 设计 (Design)

### G.1 — `cuda::md::ZeroMomentum` 新增 + 基准 stage 2/3/4 接入

`include/atomiverse/cpu/md/init.h` 中已存在 `ZeroLinearMomentum`,本规范**不再造 CPU 端**。新增 CUDA 兄弟类 `cuda::md::ZeroMomentum`(因 CUDA 域内尚无该符号,直接采用伞型路线图的命名 `ZeroMomentum`;CPU 端的 `ZeroLinearMomentum` 是历史遗留名,本规范不重命名以避免下游 ABI 漂移,但在 `cpu/md/init.h` 加一行 `using ZeroMomentum = ZeroLinearMomentum;` 同名别名,确保链路下游可以用同一个名字在两端引用)。

- `cuda::md::ZeroMomentum` 声明在 `include/atomiverse/cuda/md/init.cuh`,实现在 `src/cuda/md/init.cu`。结构:`role()==None`、`default_bindings()={end_of_step, +25}`,与 CPU 完全对齐。`launch()` 路径:device-side 两阶段 reduction 计算 `(P_cm, M_tot)` → host scalar(`cudaMemcpyAsync` D2H,然后 stream sync)→ device subtract kernel。注:单卡范围(伞型已声明 CUDA 不做 multi-GPU),所以**不需要** MPI allreduce,host scalar 直接广播回 device。构造函数接 `int freq=1`,与 CPU 一致。
- 基准 `polyethylene_tension.cpp` 在 stage 2/3/4 的 `wire_shared` 之后追加 `sim.add(cpu::md::ZeroLinearMomentum{1})`(stage 1 用 Langevin 自带漂移控制,LAMMPS 也仅在 stage 2/3/4 加 `fix momentum`;stage 5 单轴拉伸明确不加,因为 deform 需要保留 x 方向的整体平动);wiring 通过新增辅助函数 `wire_zero_momentum(cpu::md::MD&)` 集中,避免三个分支重复。
- 多 rank 正确性已由现有 `ZeroLinearMomentum::launch` 中的 `MPI_Allreduce(local{p,mtot}, global, 4, ...)` 提供。

### G.2 — 真实 Pxx/Pyy/Pzz 替换 isotropic 复制

`SystemScalars` 字段已就位(Voigt 6 维 `stress` + `stress_kinetic`),无需新增字段。基准里改两处:

- `pressure_eV_per_A3()` 旁新增 `pressure_per_axis_bar(const cpu::Frame& f, double P_axis[3])`,各分量 `(stress[k] + stress_kinetic[k]) / V * kEvPerA3ToBar`,k∈{0,1,2}。
- `CsvRecorder::launch` 中把 `Pxx/Pyy/Pzz` 三列以 bar 单位写入,**插入在 `P_bar` 之后、`E_pot_eV` 之前**(给链路下游 `canonicalize_msu_output.py` 已声明的 8 列 schema `step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz` 留出对应栏位)。

`tools/run_msu_atomiverse.sh` 注释块(行 142–148 的 "Phase G.2 TODO")和 `tools/canonicalize_msu_output.py::atomiverse_thermo` 中 "把 P 复制到 Pxx/Pyy/Pzz" 的 awk 改为直接读 atomiverse thermo 的真实 Pxx/Pyy/Pzz 三列(`$6,$7,$8`)。

### G.3 — 每项能量字段 + thermo 列

在 `core/system_scalars.h::SystemScalars` 加 4 个新 POD 字段:`double energy_pair=0, energy_bond=0, energy_angle=0, energy_dihedral=0;`(单位 eV,初值 0)。每个 evaluator 在写自己的 `+= energy_total` 时,**同时**写入自己的专属字段(累加意义一致,因为 `ZeroEnergy` 在 `pre_force` 把所有 5 个能量标量全部清零)。`ZeroEnergy` 改为同时 memset 这 4 个新字段(CPU 是普通赋值,CUDA 是 4 次 `cudaMemsetAsync` 或一次 `cudaMemsetAsync(&scalars->energy_pair, 0, 4*sizeof(prec_t))`,取决于内存连续性)。MPI 归约 `src/cpu/md/mpi/reduce.cpp` 的 14-double sendbuf 扩成 18-double(`energy + kinetic + 4 per-term + 6 stress + 6 stress_kinetic`)。

`CsvRecorder::launch` 把 `E_pair_eV / E_bond_eV / E_angle_eV / E_dihedral_eV` 4 列**插入在 `E_pot_eV` 之后、`E_kin_eV` 之前**;字段从 `frame.scalars.energy_{pair,bond,angle,dihedral}` 直接读。

### 生命周期 / 所有权

`cuda::md::ZeroMomentum` 拥有一对 device-side scratch buffer(`prec_t* d_pcm[4]` for P_cm.x/y/z + M_tot,加一个 host-pinned 副本);在 `build()` 分配,`destroy()` 释放,严格遵循 "no allocation in `launch()`" 规则。`SystemScalars` 新增字段不引入新生命周期约束(POD 聚合)。

## Files to create or modify

- `include/atomiverse/core/system_scalars.h` — 新增 `energy_pair / energy_bond / energy_angle / energy_dihedral` POD 字段
- `include/atomiverse/cpu/md/init.h` — 加 `using ZeroMomentum = ZeroLinearMomentum;` 别名
- `include/atomiverse/cuda/md/init.cuh` — 新增 `cuda::md::ZeroMomentum` 类声明
- `src/cuda/md/init.cu` — 新增 `ZeroMomentum::build/launch/destroy` 实现 + reduction kernel + `ZeroEnergy` 扩展到 4 个新字段
- `src/cpu/md/init.cpp` — `ZeroEnergy` 也清零 4 个新字段(若现存 `ZeroEnergy::launch` 只清 `scalars.energy` 则补齐)
- `src/cpu/potentials/pair_lj_cut.cpp` — `+= energy_pair` 与原 `+= energy` 并行写
- `src/cpu/potentials/bond_harmonic.cpp` — `+= energy_bond` 与原 `+= energy` 并行写
- `src/cpu/potentials/angle_harmonic.cpp` — `+= energy_angle` 与原 `+= energy` 并行写
- `src/cpu/potentials/dihedral_multi_harmonic.cpp` — `+= energy_dihedral` 与原 `+= energy` 并行写
- `src/cpu/md/mpi/reduce.cpp` — sendbuf/recvbuf 从 14 扩到 18 doubles,加 4 个新字段
- `benchmarks/md/polyethylene_tension.cpp` — 加 `pressure_per_axis_bar` 辅助、`wire_zero_momentum`、CsvRecorder header / row 扩列、stage 2/3/4 接入 ZeroLinearMomentum
- `tools/run_msu_atomiverse.sh` — 删除 142–148 Phase G.2 TODO 注释块
- `tools/canonicalize_msu_output.py` — `atomiverse thermo` 改为透传真实 Pxx/Pyy/Pzz(不再复制 P);更新 docstring
- `tests/cuda/md/test_zero_momentum.cu` (new) — CUDA ZeroMomentum 数值正确性
- `tests/cpu/md/test_pe_thermo_schema.cpp` (new) — 基准 thermo.csv 列 schema + Pxx/Pyy/Pzz 非退化 + 每项能量有限性

## Tasks

- [x] Write failing test for new SystemScalars per-term energy fields default-zero invariant (`tests/cpu/md/test_pe_thermo_schema.cpp`)
- [x] Add `energy_pair/energy_bond/energy_angle/energy_dihedral` to `include/atomiverse/core/system_scalars.h` and wire `ZeroEnergy` + 4 evaluators + `reduce.cpp` to write/clear/reduce them
- [x] Write failing test for `cuda::md::ZeroMomentum` zeroing COM linear momentum on a synthetic 1024-atom uniform-velocity device frame (`tests/cuda/md/test_zero_momentum.cu`)
- [x] Implement `cuda::md::ZeroMomentum` in `include/atomiverse/cuda/md/init.cuh` + `src/cuda/md/init.cu` (build/launch/destroy + device reduction + host-scalar broadcast subtract); add `using ZeroMomentum = ZeroLinearMomentum;` to `include/atomiverse/cpu/md/init.h`
- [x] Write failing test for benchmark thermo.csv schema: 16 columns including `Pxx/Pyy/Pzz` + 4 per-term energies, all finite, Pxx≠Pyy on a deliberately anisotropic stage-5 smoke run (`tests/cpu/md/test_pe_thermo_schema.cpp`)
- [x] Implement `pressure_per_axis_bar` and extend `CsvRecorder` header + row in `benchmarks/md/polyethylene_tension.cpp`; add `wire_zero_momentum` and call from stages 2/3/4
- [x] Update `tools/run_msu_atomiverse.sh` (remove G.2 TODO block) and `tools/canonicalize_msu_output.py` (transparent Pxx/Pyy/Pzz passthrough, drop awk replication, refresh docstring + pytest fixture)
- [x] Verify byte-faithful match against `benchmarks/md/reference/polyethylene/lammps.in` stages 2/3/4 thermo output: `fix momentum` present, `thermo_style custom` 8-column schema, energy decomposition columns positive-definite per-term where physics expects (E_bond ≥ 0, E_angle ≥ 0)
- [x] Run full check + test suite (`cmake --build build -j$(nproc)` + `ctest --test-dir build --output-on-failure`)

## Testing strategy

- **Unit — CPU per-term energy**: `test_pe_thermo_schema.cpp` constructs a 2-chain PE frame, runs one full force evaluation, asserts `scalars.energy_pair + scalars.energy_bond + scalars.energy_angle + scalars.energy_dihedral` equals `scalars.energy` to 1e-12 eV (consistency: same `+=` source).
- **Unit — CUDA ZeroMomentum**: synthetic device frame with `vx=1, vy=2, vz=3` uniform, `mass=1`; after one `ZeroMomentum::launch`, max-abs of `Σ m·v` per axis < 1e-6 eV·fs/Å (FP32 noise floor). Cover `freq=2` skipping every other step.
- **Unit — CPU ZeroLinearMomentum regression**: re-run existing tests (no new test required since CPU implementation untouched); just confirm `ctest` green.
- **Integration — benchmark schema**: end-to-end `polyethylene_tension --mode npt --scale 0.01 --print-every 10 --out-dir /tmp/petest` writes a 16-column thermo.csv whose header is `step,stage,T_target_K,T_K,P_bar,Pxx,Pyy,Pzz,E_pot_eV,E_pair_eV,E_bond_eV,E_angle_eV,E_dihedral_eV,E_kin_eV,box_x,box_y,box_z`, all rows numeric and finite, `Pxx/Pyy/Pzz` not all equal (test that the isotropic stand-in is gone).
- **Canonicalization**: `pytest tools/test_canonicalize_msu_output.py` (existing) keeps green; new fixture row covers an Atomiverse thermo with truly anisotropic Pxx/Pyy/Pzz to confirm the passthrough.
- **Domain validation**: stage-4 100 K NPT equilibrium ⟨P_bar⟩ ≈ ⟨(Pxx+Pyy+Pzz)/3⟩ to within block-statistical noise (`msu-pe-align-03` block-average framework, ±0.5 bar at scale=0.05) on the smoke run.

## Out of scope

- Off-diagonal stress components Pxy/Pxz/Pyz columns — LAMMPS reference deck does not log them; if needed downstream the SystemScalars Voigt fields already carry them, but no schema change here.
- CUDA per-term energy fields actually wired into CUDA evaluators (LJ/HarmonicBond/HarmonicAngle/MultiHarmonicDihedral on device). The 4 SystemScalars fields exist on `cuda::Frame::scalars` (the cross-backend bundle), but populating them from CUDA kernels is link `06-cuda-benchmark-path`'s job — link 1 only requires CPU evaluators write them and CUDA `ZeroEnergy` clears them.
- Multi-GPU MPI allreduce of COM momentum on the CUDA path. Umbrella scope decision §1 explicitly fixes CUDA range = single-GPU correctness only.
- LAMMPS `fix momentum angular` (rotational) — neither the reference deck nor MSU paper uses it; CPU `ZeroLinearMomentum` already documents `rot no` divergence as aligned.
- Removing the `using ZeroMomentum = ZeroLinearMomentum;` alias by renaming the CPU class — defer until the next ABI-break window; not required for chain link 1's "done" criteria.

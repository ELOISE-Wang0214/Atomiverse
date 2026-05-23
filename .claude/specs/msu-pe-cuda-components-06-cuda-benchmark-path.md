---
title: MSU PE CUDA components — link 06 — CUDA path in polyethylene_tension benchmark
status: approved
created: 2026-05-23
parent: msu-pe-validation.md (Phase G.9)
chain: msu-pe-cuda-components-* (link 6/8)
---

# MSU PE CUDA components — link 06 — CUDA path in polyethylene_tension benchmark

## Summary

把 `benchmarks/md/polyethylene_tension.cpp` 从 CPU-only 扩展为 **双后端基准**:在 `#ifdef ATV_BUILD_CUDA` 守卫下加入一条与 CPU 路径完全对位的 CUDA 分支,五个 stage 各自构造一个 `cuda::md::MD` 驱动,沿用前序链节 (01 `ZeroMomentum`、02 categorical API、03 `NVELimit`+`FixLangevin`、04 `NoseHooverNPT` iso+aniso、05 LJ 1-2/1-3/1-4 exclusions) 落地的 CUDA 组件,读相同的 `msu_pe.data`、产出相同的 thermo / stress-strain CSV schema (含 G.2 `Pxx/Pyy/Pzz` 与 G.3 `E_pair/E_bond/E_angle/E_dihedral` 列)。新 CLI 旗标 `--backend {cpu|cuda}` 在运行期选路,`tools/run_msu_atomiverse.sh` 同步增加 `--backend cuda` 传递。本链节是 Phase G.9 的实现;下游链节 08 (multi-engine validation) 在 CUDA 参考 CSV 落地之前不能跑。

## Domain basis

物理协议与 CPU 路径**逐字节相同** —— 不变量来自 `benchmarks/md/reference/polyethylene/lammps_msu_paper.in` (Tschopp 2010;M.A. Tschopp et al., *Polymer* 2010, doi:10.1016/j.polymer.2010.10.030;CAVS ICME deck)。五阶段时步 / 步数 / 系综 / 力场系数与 `msu-pe-validation.md` Phase G domain-basis 表一致:stage 1 `nve/limit 0.05` + Langevin (`T=500 K`, seed `904297`, dt 1.0 fs, 10 000 步);stage 2 iso NPT (`T=500 K`, `P=0`, dt 0.5 fs, 50 000 步) + `fix momentum`;stage 3 ramped iso NPT (500→100 K, dt 0.5 fs, 50 000 步) + `fix momentum`;stage 4 iso NPT (100 K, dt 0.5 fs, 50 000 步) + `fix momentum`;stage 5 aniso y/z NPT + `fix deform x erate 1e-5` (100 K, dt 1.0 fs, 171 800 步)。

单位:eV / Å / fs / amu / K (`core/units.h`);`1 eV/Å³ = 1.602176634×10⁶ bar = 160.21766 GPa` (常数已存在于 `polyethylene_tension.cpp` 顶部,CUDA 路径直接复用)。CUDA 是 **single-precision** 算力 (`ATV_CUDA_REAL=float`, `ATV_CUDA_WORK=double` 混合精度,见 `.claude/NOTES.md` 2026-05-22 条目),本链节**不做** Atomiverse-CUDA ↔ LAMMPS 的逐步或比特级对比 —— 那是 Phase I / 链节 08 的统计 block-average 比较 (`msu-pe-align-03`, Flyvbjerg & Petersen, *J. Chem. Phys.* 91 (1989) 461, doi:10.1063/1.457480) 的职责。本链节的物理验收只看 **CUDA 路径自洽**:stage-4 ⟨T⟩ ≈ 100 K (NoseHooverNPT 温控生效)、stage-4 ⟨P⟩ ≈ 0 bar (NoseHooverNPT 压控生效)、全协议无 NaN。

## Design

**总体策略:CPU 路径不动一行。** 现有的 `make_pe_box`、`CsvRecorder`、`TrajectoryRecorder`、`SharedComponents`、`wire_shared`、五个 `if (stage_active(...))` 块全部保留 —— 它们是 `--backend cpu` (默认) 路径,对外行为零变化。

**新增一个 `--backend` CLI 选项** (枚举 `Backend::cpu | Backend::cuda`,默认 `cpu`)。`parse_args` 在已有的 `take_s`/`take_d`/`take_i` lambda 之后加一段 `if (std::strcmp(s, "--backend") == 0)` 分支,与 `--mode` 同款解析。`--backend cuda` 在 **未编 `ATV_BUILD_CUDA`** 的二进制上立即报错退出 (返回码 2),不静默回退到 CPU —— 用户的选择必须被尊重。

**CUDA 分支以 `#ifdef ATV_BUILD_CUDA` 块整体出现在 `main()` 内**,与 CPU 五个 stage 块互斥 (按 `args.backend` 分支)。块内结构:

1. **共享 CUDA 组件 (lvalue, outlive every per-stage driver):**
   - `cuda::Frame d_frame;` — device-resident frame。
   - `cpu::to_cuda(frame, d_frame);` — 一次 H2D,把 `cpu::load_lammps_data(args.data_path)` 得到的 `cpu::Frame` 拷到设备。
   - `cuda::NeighborList nl{kLjCutoff};`
   - `cuda::LJ lj{lj_eps_eV(), kLjSigma, kLjCutoff};` 接 `set_special_bonds(0.0, 0.0, 0.0)` (链节 05 落地的 API)。
   - `cuda::HarmonicBond hb`、`cuda::HarmonicAngle ha`、`cuda::MultiHarmonicDihedral hd` —— 与 CPU 同力场常数。
   - **CSV recorder:** 一个新的 `CudaCsvRecorder`,类内通过 `AsyncBridge` 在 `defer_stream` 上 D2H 一份 `SystemScalars` (含 `stress`, `stress_kinetic`, `energy_pair/bond/angle/dihedral`, `kinetic`, `box`),host 回调里写入与 CPU 路径**字段顺序与列名完全相同**的 `thermo.csv` / `stress_strain.csv`。`begin_stage(stage_id, T_target, record_strain)` 与 CPU 同签名,保持 CSV `step` 列跨 stage 单调。
   - 重用 CPU 路径的同名 `std::ofstream thermo_csv, stress_csv;` —— 两个后端共用一个文件 (`--backend` 二选一,互不竞争)。

2. **`wire_shared_cuda(cuda::md::MD& sim, CudaSharedComponents& sc)` 帮助函数**,与 CPU `wire_shared` 同形:通过 G.4 的 categorical API `sim.set_neighbor / set_pair / set_bond / set_angle / set_dihedral` 注册评估器,把 `CudaCsvRecorder` 通过 `sim.set_record(...)` 或 `sim.on(md::Phase::end_of_step).add(rec, 25)` 绑到 `end_of_step` slot (与 CPU 一致的 `+25` order)。`VirialStress` 由 `cuda::md::MD::build()` 在 `set_integrator(NoseHooverNPT{...})` 触发 `reads()==Stress` 时自动注入 (categorical-API 的 auto-inject 契约,链节 02 行为承诺)。

3. **五个 stage 块** —— 与 CPU 块一一对应、签名相同:
   - Stage 1:`cuda::md::NVELimit nve{dt_stage1, 0.05}`;`cuda::md::InitVelocity init_v{T_hot, seed}`;`cuda::md::FixLangevin langevin{T_hot, T_hot, langevin_damp, dt_stage1, seed+1}`;`sim.set_integrator(nve); sim.add(init_v); sim.set_thermostat(langevin);` 不带 `ZeroMomentum`。
   - Stage 2:`cuda::md::NoseHooverNPT npt{{.T_target=T_hot, .tau_T, .P_target=0.0, .tau_P, .dt=dt_stage2, .drag=2.0}}`;`cuda::md::ZeroMomentum zm{};`;`sim.set_integrator(npt); sim.add(zm);`。
   - Stage 3:与 stage 2 相同 NPT 但 `.T_stop=T_cold`,附加 `npt.set_run_span(0, steps-1);`;带 `ZeroMomentum`。
   - Stage 4:iso NPT @ `T_cold`;带 `ZeroMomentum`。
   - Stage 5:aniso NPT (`.couple_axis = {false, true, true}`) + `cuda::md::Deform deform{cuda::md::Deform::Mode::erate, strain_rate, dt_stage5, 0}`;**不带** `ZeroMomentum` (与 LAMMPS deck 一致)。

4. **生命周期:** `cuda::Frame` 与 CPU 路径里的 `cpu::Frame` 角色一致 —— 它是 stage 间结转的设备状态;五个 per-stage `cuda::md::MD` 各自 `build() → run(d_frame, steps) → destroy()`,但 `d_frame` 持续存在。run 全部结束后 `cuda::to_cpu(d_frame, frame_final);` 是 **可选** 的 (CSV 已经在 `CudaCsvRecorder` 里写完了),本链节不强制做最终 D2H 拉回。

5. **`tools/run_msu_atomiverse.sh` 扩展:** 在脚本顶部解析 `--backend cuda` (在 `RANKS` 之后;用 case-shift 而非 getopts,与脚本现有风格一致),把它转发给 `polyethylene_tension --backend cuda`。`--backend cuda` 时 `LAUNCH=()` (单进程 GH200),禁用 `mpirun`/`srun` 多进程 (RANKS 强制 1)。`RUN_DIR` 后缀从 `atomiverse_${RANKS}r_${TS}` 改为 `atomiverse_${BACKEND}_${RANKS}r_${TS}` 以避免与 CPU 运行混叠。Canonical artifact 文件名同步加 backend 标签:`thermo_atomiverse_cuda_${RANKS}r.csv` 等 —— 下游 `canonicalize_msu_output.py` 的 `engine` 标签暂传 `atomiverse` (本链节不改 canonicalizer),但文件名前缀区分。

**显式不做的事:** 不重写 CPU 路径;不引入 ZeroMomentum 的 CPU/CUDA 实现 (链节 01 的事);不实现 categorical setter API (链节 02 的事);不实现 NVELimit/FixLangevin/NoseHooverNPT/LJ-exclusions 的 CUDA kernel (链节 03/04/05 的事);不做 CUDA↔CPU / CUDA↔LAMMPS 的统计比较 (链节 08 / Phase I);不动 `canonicalize_msu_output.py`。

## Files to create or modify

- `/cfs/klemming/home/j/jichenli/work/Atomiverse/benchmarks/md/polyethylene_tension.cpp` — 增 `Backend` 枚举、`--backend` 解析、`#ifdef ATV_BUILD_CUDA` CUDA 路径 (含 `CudaCsvRecorder`、`CudaSharedComponents`、`wire_shared_cuda`、五个 stage 块)。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/benchmarks/CMakeLists.txt` — 在 `atv_add_cpp_md_bench(polyethylene_tension)` 之后,加 `if(ATV_BUILD_CUDA) target_compile_definitions(polyethylene_tension PRIVATE ATV_BUILD_CUDA) target_link_libraries(polyethylene_tension PRIVATE atomiverse_cuda) endif()` (验证 `atomiverse_cuda` target 名 — 若实际名为 `atomiverse-cuda` 或 `atv_cuda` 则用之;可由实现者在编译时落地)。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tools/run_msu_atomiverse.sh` — 加 `--backend {cpu|cuda}` 解析、CUDA 时禁用 `mpirun`/`srun`、`RUN_DIR` 与 `THERMO_REF`/`STRESS_REF` 文件名带 backend 标签。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tests/benchmarks/test_polyethylene_tension_cuda_smoke.sh` (new) — gtest-free PASS/FAIL shell 测试,按 `benchmarks/CMakeLists.txt` 中 "gtest-free PASS/FAIL harness" 模式 (注释 73-74 行):编译 `ATV_BUILD_CUDA=ON` 后跑 `polyethylene_tension --backend cuda --mode full --scale 0.01 --data benchmarks/md/reference/polyethylene/msu_pe.data --out-dir <tmp>`,校验 `thermo.csv` 存在、列名匹配 canonical schema、所有数值有限 (无 NaN/Inf)。CTest 在 `ATV_BUILD_CUDA=ON` 时注册。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tests/CMakeLists.txt` — 在 `if(ATV_BUILD_CUDA)` 块下注册 `test_polyethylene_tension_cuda_smoke` (`add_test(NAME ... COMMAND bash ...)`),仅在 GPU 可见时才跑 (用 `set_tests_properties(... PROPERTIES LABELS "cuda;gpu_required")`)。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/.claude/specs/INDEX.md` — 加一行 `msu-pe-cuda-components-06-cuda-benchmark-path`。

## Tasks

- [ ] Write failing smoke test `tests/benchmarks/test_polyethylene_tension_cuda_smoke.sh` (new) that compiles + runs `polyethylene_tension --backend cuda --mode full --scale 0.01` and asserts CSV schema match + finite values
- [ ] Implement `--backend {cpu|cuda}` CLI parsing in `benchmarks/md/polyethylene_tension.cpp` (`Backend` enum + parse-arg branch + help-text line) with CPU default preserving byte-identical behavior
- [ ] Implement `CudaCsvRecorder` and `wire_shared_cuda` in `benchmarks/md/polyethylene_tension.cpp` under `#ifdef ATV_BUILD_CUDA`, emitting the canonical schema columns (`step,stage,T_target,T,P_bar,Pxx,Pyy,Pzz,E_pair,E_bond,E_angle,E_dihedral,E_pot,E_kin,box_x,box_y,box_z`)
- [ ] Implement the five `cuda::md::MD` per-stage blocks in `benchmarks/md/polyethylene_tension.cpp` mirroring the CPU stages 1–5 via the G.4 categorical API + G.5–G.8 CUDA components + G.1 `cuda::md::ZeroMomentum` on stages 2/3/4
- [ ] Wire `ATV_BUILD_CUDA` compile-definition and CUDA-library link into `benchmarks/CMakeLists.txt` for the `polyethylene_tension` target
- [ ] Extend `tools/run_msu_atomiverse.sh` with `--backend {cpu|cuda}` (default cpu), force single-process when cuda, and tag `RUN_DIR` + canonical-artifact filenames with the backend
- [ ] Register `test_polyethylene_tension_cuda_smoke` in `tests/CMakeLists.txt` under `if(ATV_BUILD_CUDA)` with `LABELS "cuda;gpu_required"`
- [ ] Add `msu-pe-cuda-components-06-cuda-benchmark-path` entry to `.claude/specs/INDEX.md`
- [ ] Verify against `--backend cuda --mode full --scale 0.01` smoke run on a GH200: thermo CSV has all 5 stages, finite values, ⟨T⟩ stage-4 within ±10 K of 100 K, ⟨P⟩ stage-4 within ±200 bar of 0
- [ ] Add Doxygen-style banner to the new CUDA section per `doc.style: doxygen`, document units (`T [K]`, `P [bar]`, `dt [fs]`, `strain_rate [1/fs]`, `epsilon [eV]`, `sigma [Å]`)
- [ ] Run full check + test suite (`cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON -DATV_BUILD_CUDA=ON -DATV_BUILD_BENCHMARKS=ON && cmake --build build -j && ctest --test-dir build --output-on-failure -L cuda`)

## Testing strategy

- **Happy path (CPU regression):** `polyethylene_tension --mode full --scale 0.01 --data msu_pe.data` (no `--backend`, default cpu) and `... --backend cpu ...` both produce **byte-identical** `thermo.csv` and `stress_strain.csv` — proves the `--backend` plumbing does not perturb the existing CPU path.
- **Happy path (CUDA smoke):** `--backend cuda --mode full --scale 0.01` on GH200 completes; thermo CSV has the 17-column canonical schema (post G.2 + G.3), one row per `print_every` steps, `step` column monotonic across stages, all numeric cells finite (`isfinite()`); stage-4 ⟨T⟩ within ±10 K of 100 K (NoseHooverNPT temperature control), stage-4 ⟨P⟩ within ±200 bar of 0 bar (NoseHooverNPT pressure control); stage-5 `strain` column monotonically increases.
- **Edge case (build):** `cmake -DATV_BUILD_CUDA=OFF` build succeeds and `--backend cuda` exits with code 2 and a clear error message (`"cuda backend not compiled in; rebuild with -DATV_BUILD_CUDA=ON"`), no silent fall-through to cpu.
- **Edge case (mode subset):** `--backend cuda --mode stage:1 --scale 0.05` runs stage 1 only (NVELimit + FixLangevin) without NaN — proves stochastic stage 1 alone is sound on CUDA before the full pipeline.
- **Edge case (script):** `tools/run_msu_atomiverse.sh --backend cuda --help` is rejected only if cluster lacks CUDA; otherwise the script generates a `RUN_DIR` matching `atomiverse_cuda_1r_*` and writes `thermo_atomiverse_cuda_1r.csv` distinct from any CPU run dir.
- **Domain validation:** per Domain basis — single-precision CUDA does not match double-precision CPU bit-for-bit; this spec validates only **CUDA self-consistency** (finite values + NoseHooverNPT temperature/pressure control reaching steady-state setpoints in stage 4). Cross-engine statistical comparison is link 08 / Phase I.
- **No new unit test under `tests/cuda/md/`:** components (NVELimit, FixLangevin, NoseHooverNPT, LJ exclusions, ZeroMomentum) are tested in their own chain links (01, 03, 04, 05); this link tests only the *integration*, which is the smoke run.

## Out of scope

- ZeroMomentum, NVELimit, FixLangevin, NoseHooverNPT, LJ exclusions, categorical-API setter implementations (chain links 01–05).
- Async correctness for the CUDA recorder hot-path (Phase H.2 / H.3 — a `defer_stream` snapshot pattern is implied by this spec's `CudaCsvRecorder` design but the *correctness* gate (recorder frequency must not stall the loop) is enforced in Phase H, not here).
- Statistical CUDA↔CPU or CUDA↔LAMMPS-GPU comparison (Phase I / chain link 08).
- `canonicalize_msu_output.py` extension for the CUDA engine label — currently the script passes `engine=atomiverse`; if downstream Phase I needs `atomiverse_cuda` as a distinct engine label, the canonicalizer change rides with link 08, not here.
- PDC `gpugh` partition sbatch wrapper around `run_msu_atomiverse.sh --backend cuda` — operator-action artifact for Phase F.
- LAMMPS `lmp_cuda` GPU-mode reference run.
- Multi-GPU / CUDA+MPI: explicitly out per `msu-pe-validation.md` fixed-scope decision #1.
- Removing the `# Phase G.2 TODO` comment in `run_msu_atomiverse.sh` about Pxx/Pyy/Pzz being faked — G.2 itself lands in another chain; this link only adds the `--backend` plumbing.

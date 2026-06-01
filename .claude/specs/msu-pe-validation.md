---
title: MSU amorphous-PE tension -- CPU MD correctness and performance validation vs LAMMPS
status: in-progress
created: 2026-05-22
last-updated: 2026-06-01 (DESCOPE to CPU MD only -- user dropped the CUDA/GPU backend, the async-CUDA work, and the Python bindings from this validation. Removed Phase I (CUDA single-GPU), the CUDA arms of Phases G/H/D/E/F, and the four-config 2x2 matrix; folded the surviving CPU work into a single-backend roadmap. Acceptance contract re-cut from 21 dual-backend criteria to 9 CPU-only criteria. Re-verified the CPU state by build+run: build-main (x86) green; test_msu_reference_inputs 11/13 pass, 8-rank test_pe_msu_pipeline PASSED, polyethylene_tension --mode full runs all 5 stages on the real 10k-atom msu_pe.data with no NaN. Prior dual-backend spec is retained in git history.)
supersedes: msu-pe-validation (in place)
---

# MSU amorphous-PE validation (CPU MD)

## 综述 (Summary)

本规范验证 Atomiverse 的 **CPU/MPI MD 后端** 在 MSU 非晶聚乙烯单轴拉伸基准上对 LAMMPS 的正确性与性能。最终诉求:用与 `lammps_msu_paper.in` 逐 fix / 逐 ensemble / 逐步数 **物理上完全一致** 的 5 阶段协议,从同一个 `msu_pe.data` 出发跑完整轨迹,产出能量(含分解)、密度、应力-应变曲线和每阶段墙钟时间,并与 LAMMPS-CPU 参考在 **统计分块平均** 意义上一致。凡缺失的功能先实现;非热路径的 recorder/thermo 改为异步,不阻塞积分线程。

本规范保持 **伞型路线图 (umbrella roadmap)** 形态:持有阶段状态快照、依赖图、各阶段退出门。Phase 0/A/B/C(引擎阻塞、单/多 rank 正确性、体系构建)、功能补齐 G(CPU 部分)、异步正确性 H.1/H.4、性能 E.8、复现工具 F(CPU 配置)全部 ✅ landed。剩余 **operator-action**:全尺寸 AV-CPU 参考运行(Phase A 收尾)+ Phase D 的 CPU↔LAMMPS-CPU 统计校验。

> **2026-06-01 DESCOPE.** CUDA/GPU 后端、async-CUDA 工作、Python 绑定均移出本规范范围,详见 `Out of scope`。CUDA 侧的工具/代码已落地的部分仍在 git 历史与 `msu-pe-cuda-components-*` 链中,本规范不再跟踪。

## 已确认的范围决定 (fixed scope decisions)

1. **后端 = 仅 CPU/MPI。** 验证只针对 CPU 后端在单/多 rank 下的正确性与扩展性。CUDA / 多 GPU / CUDA+MPI **不在范围内**。
2. **"每阶段速度" = 每阶段墙钟时间。** 给基准加每阶段墙钟计时,用于性能描述。
3. **对比 = AV-CPU ↔ LAMMPS-CPU。** 单条参考轨道:Atomiverse-CPU ↔ `lmp` CPU 模式。协议是随机性的(stage 1 是带种子的 Langevin),**所有跨引擎数值对比一律用集合均值的统计分块平均**(Flyvbjerg & Petersen 1989,遵循 `msu-pe-align-03`);不做对 LAMMPS 的逐步或比特级对比。
4. **集群 = PDC Dardel `main` (CPU) 分区, SLURM。** 完整 ~360k 步协议是操作员动作;运行脚本通过 `sbatch`/`srun` 提交,需要 `-A <project>` 账户标志。

## 状态快照 (Status snapshot)

| Phase | Sub-item | Status |
|---|---|---|
| 0 / A / B / C | 引擎阻塞、单/多 rank 正确性、MSU 体系构建(多 rank 1e-10 能量匹配) | ✅ landed |
| G (CPU) | 保真缺口:`fix momentum`(`cpu::md::ZeroLinearMomentum`)、真实 per-axis `Pxx/Pyy/Pzz`、能量分解 `E_pair/E_bond/E_angle/E_dihedral` | ✅ landed |
| H.1 | CPU `Thermo` 改为 `AsyncComponentBase`,reduction+I/O 走 `AnalysisExecutor` | ✅ landed |
| H.4 (CPU) | 异步正确性退出门:Thermo-on 墙钟在 Thermo-off 的 5% 以内 | ✅ landed |
| A (收尾) | **全尺寸 AV-CPU 5 阶段参考运行**(`--scale 1.0`,10k 原子,~332k 步) | 🟡 operator-action — 机制已验证(`--scale 0.001` smoke 跑完 5 段无 NaN),只差全尺寸跑一遍出平衡数据 |
| D (CPU) | AV-CPU ↔ LAMMPS-CPU 统计分块平均对比 | 🟡 工具 ✅ landed(`tools/compare_msu.py`);差 operator 的 LAMMPS-CPU 参考运行 |
| E.2 / E.3 | CPU 强/弱扩展性 1–32 ranks | ⏳ external — 需 >8 rank 的集群分配 |
| E.8 | 每阶段墙钟计时 | ✅ landed(`polyethylene_tension.cpp`,commit `4c3e109`) |
| F (CPU) | 复现工具:`run_msu_atomiverse.sh`、`run_msu_lammps.sh --mode cpu`、`run_msu_pdc.sh`、`compare_msu.py` | ✅ tooling landed;operator 复现运行仍是外部动作 |

**Effective state.** CPU 代码侧全部 landed 且已 build+run 验证(见下)。剩余两步都是 operator-action 且**不需要 GPU**:(1) 全尺寸 AV-CPU 跑出平衡参考数据(ac-007);(2) LAMMPS-CPU 参考运行 + `compare_msu.py` 统计判定(ac-008)。

**2026-06-01 build+run 验证(x86 `build-main`):**
- `tests/cpu/md/test_msu_reference_inputs` → 11 passed / 2 skipped(2 skip 仅因尚无真机 run 产出的参考 CSV)。
- `tests/mpi/test_pe_msu_pipeline`(8-rank 全栈)→ PASSED(215 ms)。
- `polyethylene_tension --mode full --scale 0.001` 在真实 10k 原子 `msu_pe.data` 上跑完全部 5 段,exit 0,无 NaN,thermo 含完整 `Pxx/Pyy/Pzz` + 能量分解列。

---

## Domain basis

The protocol physics is the immutable upstream deck `benchmarks/md/reference/polyethylene/lammps_msu_paper.in` (Tschopp 2010; M.A. Tschopp et al., "Molecular dynamics simulations of deformation in amorphous polyethylene," *Polymer* 2010; CAVS ICME source <https://www.cavs.msstate.edu/icme/code/lammps/atomistic_research/polymers/amorphous.php>). System: 10 chains x 1000 united-atom CH2 = 10 000 atoms, cubic box 80.0586 A, mass 14.02 amu per UA bead (`benchmarks/md/reference/polyethylene/SOURCES.md`).

Force field -- Dreiding-UA, verbatim from the deck: harmonic bond *k* = 350 kcal/mol/A^2, *r*_0 = 1.53 A; harmonic angle *k* = 60 kcal/mol/rad^2, theta_0 = 109.5 deg; multi/harmonic dihedral *A*_1..*A*_5 = 1.73, -4.49, 0.776, 6.99, 0.0 kcal/mol; LJ eps = 0.112 kcal/mol, sigma = 4.01 A, *r*_cut = 10.5 A. Units: eV, A, fs, amu, K (`core/units.h`); LAMMPS `units real` -> Atomiverse eV via `1 kcal/mol = 1/23.060547830619026 eV`.

Five-stage protocol (every fix, ensemble, timestep and step count is byte-faithful ground truth):

| Stage | Ensemble / fix | T | P | Steps | dt |
|---|---|---|---|---|---|
| 1 | `nve/limit 0.05` + `langevin 500 500 10.0 904297` | 500 K | -- | 10 000 | 1.0 fs |
| 2 | `npt ... iso` + `momentum 1 linear 1 1 1` | 500 K | 0 bar | 50 000 | 0.5 fs |
| 3 | `npt temp 500 100 ... iso` + `momentum` (cooling ramp) | 500->100 K | 0 bar | 50 000 | 0.5 fs |
| 4 | `npt ... iso` + `momentum` | 100 K | 0 bar | 50 000 | 0.5 fs |
| 5 | `npt ... y ... z` (aniso) + `deform 1 x erate 1e-5` | 100 K | 0 bar (y,z) | 171 800 | 1.0 fs |

Stage 1 is load-bearing in two ways: (a) `msu_pe.data` is a built-but-**unrelaxed** configuration with non-bonded overlaps (initial `E_pair ~ 3.4e8 kcal/mol`) -- a displacement-capped `nve/limit` integrator is mandatory, a plain Velocity-Verlet integrator blows the system up; (b) `fix langevin` is **seeded** (`904297`), so stage 1 is a **stochastic** thermostat. Per `msu-pe-align-03` (block averaging -- Flyvbjerg & Petersen 1989, *J. Chem. Phys.* 91, 461, doi:10.1063/1.457480), two engines running the same stochastic-thermostat protocol with different RNG streams do **not** share a trajectory -- they sample the same ensemble. **Every cross-engine value comparison in this spec is therefore a statistical block-average of ensemble means; no per-step or bit-level comparison against LAMMPS.**

Stress for the stress-strain curve: engineering strain eps = box_x/box_x0 - 1; sigma_xx in GPa (`1 eV/A^3 = 160.21766 GPa`). Mechanical descriptors: elastic modulus *E* = dsigma/deps over eps < 0.02, yield sigma_y = max sigma_xx over eps in [0,0.1], hardening slope dsigma/deps over eps in [0.5,1.0].

### Strict reference-fidelity requirement

The benchmark protocol must be **byte-faithful in physics** to `lammps_msu_paper.in`: every fix, ensemble, integrator, timestep, step count and force-field coefficient. The three fidelity defects called out historically are all fixed (Phase G):

- **`fix momentum`** -- the deck applies `fix momentum 1 linear 1 1 1` in stages 2/3/4; implemented as `cpu::md::ZeroLinearMomentum` (alias `ZeroMomentum`), bound at `end_of_step`, multi-rank-correct via `MPI_Allreduce`.
- **Per-axis pressure** -- thermo emits real diagonal-virial `Pxx/Pyy/Pzz`; `tools/run_msu_atomiverse.sh` no longer copies isotropic P into the per-axis columns.
- **Energy decomposition** -- thermo emits `E_pair/E_bond/E_angle/E_dihedral`, reduced across ranks with the rest of the scalars.

The per-stage timesteps and step counts in `polyethylene_tension.cpp` already match the reference (verified -- stage 5 is dt = 1.0 fs, 171 800 steps).

---

## Dependency graph

```
Phase 0/A/B/C ✅ (engine blockers, FF correctness, system construction)
        |
Phase G  ✅ (fix momentum, per-axis P, energy decomposition)
        |
Phase H.1 ✅ (CPU Thermo async)  ----  H.4 ✅ (async-correctness gate)
        |
Phase A收尾  🟡 operator: full-scale AV-CPU run (ac-007)
        |
Phase D-CPU  🟡 operator: LAMMPS-CPU run + compare_msu.py (ac-008)
        |
Phase E.8 ✅ (per-stage timing)   F ✅ (CPU reproducibility tooling)
```

---

## Phase G -- Functional-gap implementation (CPU) ✅ landed

- **ZeroMomentum (CPU).** `cpu::md::ZeroLinearMomentum` (alias `ZeroMomentum`) in `include/atomiverse/cpu/md/init.h` -- bound at `end_of_step @ +25`, multi-rank-correct via 4-double `MPI_Allreduce` of `{Sum m, Sum m*v_d}`. Benchmark stages 2/3/4 wired.
- **Per-axis `Pxx/Pyy/Pzz`.** `polyethylene_tension.cpp` thermo emits real diagonal-virial columns; `tools/run_msu_atomiverse.sh` no longer fakes them.
- **Energy decomposition.** thermo CSV emits `E_pair/E_bond/E_angle/E_dihedral`; reduced across ranks in the same `MPI_Allreduce` as the rest of the scalars.

Exit gate (a)-(c) satisfied.

---

## Phase H.1 / H.4 -- Async correctness (CPU) ✅ landed

H.1: CPU `cpu::recorder::Thermo` subclasses `AsyncComponentBase` and routes reduction + I/O through `AnalysisExecutor`, off the integration thread. H.4 exit-gate micro-benchmark: CPU steps/s with Thermo on (high frequency) within 5% of Thermo off.

---

## Phase A收尾 / D -- CPU reference run + LAMMPS comparison (operator-action)

Goal: demonstrate Atomiverse-CPU and LAMMPS-CPU, starting from the same `msu_pe.data` and running the same protocol, produce statistically matching observables.

- **ac-007 (full-scale AV-CPU run).** `polyethylene_tension --mode full --scale 1.0` on the real 10k-atom `msu_pe.data` (PDC `main`) completes all 5 stages with no NaN; post-stage-4 density within 5% of ~0.94 g/cm3. The mechanism is verified (`--scale 0.001` smoke runs all 5 stages clean); only the full-length equilibration run remains.
- **ac-008 (CPU<->LAMMPS-CPU comparison).** `tools/run_msu_lammps.sh --mode cpu` produces the LAMMPS reference; `tools/compare_msu.py` does the statistical block-average parity (Flyvbjerg & Petersen 1989) over ⟨T⟩, ⟨E_pot⟩ + decomposition, density (stage 4), and the three stress-strain descriptors (E, sigma_y, hardening slope), tolerances ~3% on means / +-5% on descriptors.

**Exit gate**: ac-007 completes with correct density; ac-008 reports PASS and exits 0.

---

## Phase E -- Performance benchmarks (CPU)

CPU/MPI scaling (E.1-E.7: single-rank breakdown, strong/weak scaling, head-to-head vs LAMMPS-CPU, profile decomposition, hybrid MPI x OpenMP, multi-node) stays as in the prior spec; E.2/E.3 at ranks > 8 are cluster-gated (⏳ external).

### E.8 Per-stage wall-clock timing ✅ landed

Per-stage wall-clock timer in `benchmarks/md/polyethylene_tension.cpp` brackets each of the five stages and emits a per-stage timing CSV (commit `4c3e109`).

---

## Phase F -- Reproducibility artifacts (CPU, PDC Dardel)

- **F.1 Inputs under version control** -- `benchmarks/md/reference/polyethylene/` (`force_field.yaml`, `msu_pe.data`, `lammps_msu_paper.in`, reference CSVs, `SOURCES.md`, `README.md`).
- **F.2 Scripts** -- `tools/run_msu_lammps.sh --mode cpu` (lmp CPU mode); `tools/run_msu_atomiverse.sh` (CPU); `tools/run_msu_pdc.sh` (PDC SLURM submission, `-A` account flag required, `main` partition); `tools/compare_msu.py` + `tools/plot_msu_compare.py` (block-average validator + overlay).
- **F.3 README** -- `benchmarks/md/reference/polyethylene/README.md` documents the CPU-config reproduction (lmp CPU build/usage, PDC `main` `sbatch` with `-A`, expected-results table). Reference: <https://support.pdc.kth.se/doc/run_jobs/job_scheduling/>.

**Exit gate**: a fresh clone + the documented operator commands reproduce the AV-CPU and LAMMPS-CPU runs on PDC Dardel and `tools/compare_msu.py` exits 0.

---

## Tasks

CPU implementation tasks are ✅ landed; remaining work is operator-action SLURM runs on PDC Dardel `main`.

- [x] Implement `cpu::md::ZeroLinearMomentum` (fix momentum), wire benchmark stages 2/3/4 (Phase G)
- [x] Implement real per-axis Pxx/Pyy/Pzz + E_pair/E_bond/E_angle/E_dihedral CSV columns (Phase G)
- [x] Implement CPU Thermo on `AsyncComponentBase` (Phase H.1)
- [x] Async-correctness gate: Thermo-on within 5% of Thermo-off (Phase H.4)
- [x] Per-stage wall-clock timer (Phase E.8)
- [x] CPU comparison tooling -- `tools/compare_msu.py` + `tools/plot_msu_compare.py`; `run_msu_*.sh` CPU paths (Phase F)
- [ ] **Operator-action.** Submit the full-scale AV-CPU 5-stage run (`--scale 1.0`) on PDC `main`; confirm all stages complete, no NaN, density within 5% of ~0.94 g/cm3 (ac-007)
- [ ] **Operator-action.** Submit the LAMMPS-CPU reference run on PDC `main`, then run `tools/compare_msu.py` and confirm the block-average invariants PASS (ac-008)

## Testing strategy

- **Happy path** -- `tests/cpu/md/test_zero_momentum.cpp` drives COM momentum to zero; `tests/cpu/md/test_msu_reference_inputs` checks FF constants + thermo schema; the 8-rank `tests/mpi/test_pe_msu_pipeline` runs the full bonded + LJ special-bonds + halo + migration + NVT + barostat + deform stack for 1000 steps without NaN.
- **Edge cases** -- `NVELimit` on the overlapped `msu_pe.data` (initial E_pair ~ 3.4e8 kcal/mol) does not blow up; `tools/compare_msu.py` with a missing engine CSV exits non-zero with a clear error.
- **Domain validation** -- stochastic-protocol parity (Flyvbjerg & Petersen 1989 block averaging): the AV-CPU <-> LAMMPS-CPU comparison uses block-averaged ensemble means (⟨T⟩, ⟨E_pot⟩, energy decomposition, density) over the equilibrated stages, never per-step. Stress-strain descriptors (E, sigma_y, hardening slope) agree to +-5%. Reference fidelity: an assertion confirms the benchmark applies a `fix momentum` equivalent in stages 2/3/4 and emits real (non-faked) `Pxx/Pyy/Pzz`.
- **Async correctness** -- the H.4 micro-benchmark measures steps/s with Thermo off vs at high frequency; the gate requires high-frequency steps/s within 5% of Thermo-off.

## Out of scope

- **CUDA / GPU backend (descoped 2026-06-01).** The entire CUDA arm -- Phase I (single-GPU 5-stage correctness), the CUDA components/categorical API, async-CUDA hot-path desync (former H.2/H.3), the GPU-mode LAMMPS comparison, and the four-config 2x2 matrix -- is removed from this spec per the user's CPU-only descope. Whatever CUDA code/tooling already landed lives in git history and the `msu-pe-cuda-components-*` chain; it is not tracked here.
- **Python bindings (descoped 2026-06-01).** Not part of the CPU MD validation.
- **Multi-GPU and CUDA+MPI scaling** -- out of scope (follows from the CPU-only descope).
- **Per-step or bit-level comparison against LAMMPS** -- physically impossible for this stochastic protocol; all cross-engine comparison is statistical block-average (Flyvbjerg & Petersen 1989).
- **Matching MSU paper numbers directly** -- LAMMPS under identical inputs is the reproducible reference; MSU paper-number matching is a downstream addendum.
- **Changing the protocol physics** -- `lammps_msu_paper.in` is immutable ground truth; the benchmark is corrected to match it, never the reverse.

## Rationale

- **Why descope to CPU?** The user explicitly narrowed the goal to the CPU MD backend and dropped the GPU backend and Python bindings. Carrying the dual-backend roadmap added 12 GPU-only criteria and a four-config matrix that no longer reflect intent; collapsing to a single CPU backend makes "done" reachable and the contract honest.
- **Why call out `fix momentum`, per-axis P, and energy decomposition as fidelity defects?** The benchmark must strictly match `lammps_msu_paper.in`: omitting `fix momentum` is a physics deviation, faking `Pxx/Pyy/Pzz` is wrong for the stress-strain curve, and the energy decomposition is what the user compares. These were correctness defects, now fixed.
- **Why statistical block-average for the LAMMPS comparison?** The MSU protocol is stochastic (seeded stage-1 Langevin); two engines never share a trajectory -- they sample the same ensemble. The only physically meaningful cross-engine statement is agreement of block-averaged ensemble means (Flyvbjerg & Petersen 1989).
- **Why PDC Dardel `main`?** The full-length reference runs are long operator actions; PDC Dardel's `main` CPU partition runs both the AV-CPU and LAMMPS-CPU references, and the wall-time cap is why the runs are operator actions, not CI steps.

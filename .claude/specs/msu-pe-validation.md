---
title: MSU amorphous-PE tension — dual-backend (CPU + CUDA) correctness and performance validation vs LAMMPS
status: in-progress
created: 2026-05-22
last-updated: 2026-05-22 (SUPERSEDE/REFINE of msu-pe-validation — extends the umbrella to dual-backend: adds single-GPU CUDA correctness validation against LAMMPS GPU mode, a functional-gap phase that implements the missing CUDA components, an async-correctness phase, and a 2×2 four-config comparison. CPU Phase 0/A/B/C remain ✅ as landed. Phase D sub-chain `msu-pe-phase-d-*` and `msu-pe-align-03` are referenced, not re-specified.)
supersedes: msu-pe-validation (in place)
---

# MSU amorphous-PE validation (dual-backend)

## 综述 (Summary)

本规范是 `msu-pe-validation` 的就地 SUPERSEDE/REFINE：把原本只针对 CPU/MPI 引擎的 MSU 非晶聚乙烯单轴拉伸验证，扩展为 **CPU 与 CUDA 双后端** 对 LAMMPS 的验证。用户的最终诉求是：MSU 非晶 PE 拉伸基准要同时验证 Atomiverse 的 CPU 后端和 CUDA 后端,并在四种引擎配置上对比能量、密度、每阶段墙钟时间和应力-应变曲线;凡是缺失的功能必须**先实现**;所有非热路径的 recorder / thermo / observer 必须改为异步,尤其在 CUDA 上 GPU 绝不能因为 analyzer/recorder 而停等。

本规范保持 **伞型路线图 (umbrella roadmap)** 形态:它持有阶段状态快照、依赖图、各阶段的退出门 (exit gate),并把详细的可实现切片继续下放给各阶段子链 (Phase D → `msu-pe-phase-d-*`)。它不被拆成链 —— 用户明确要求更新这一个文件。CPU 侧的 Phase 0 / A / B / C 已全部完成,其 ✅ 状态原样保留;新工作以新阶段的形式编织进来:功能补齐阶段 (G)、异步正确性阶段 (H)、CUDA 单卡正确性验证阶段 (I),并把 Phase D / E / F 升级为双后端、四配置。用户可见的成果:一次基准运行同时给出 CPU 与 CUDA 对 LAMMPS 的统计学一致性结论,以及四引擎的逐阶段性能对比,且 GPU 永不为 recorder 停等。

## 已确认的范围决定 (fixed scope decisions)

1. **CUDA 范围 = 单卡正确性验证。** CUDA 在单张 GPU 上跑完整 5 阶段协议;结果与 LAMMPS GPU 模式、以及 Atomiverse CPU 运行对比。多 GPU / CUDA+MPI 扩展性 **不在范围内** —— 现有的 CPU/MPI Phase E (强/弱扩展性, E.1–E.7) 保持 CPU-only 不变。
2. **"每阶段速度" = 每阶段墙钟时间。** 给基准加每阶段墙钟计时,在四种引擎配置间对比。这是计时/性能对比,不是原子速度的物理观测量。
3. **对比是同类相比的 2×2。** 用户构建了单个 LAMMPS 二进制 `lmp_cuda`,可跑 CPU 模式和 GPU 模式两种模式。两条参考轨道:Atomiverse-CPU ↔ `lmp_cuda` CPU 模式;Atomiverse-CUDA ↔ `lmp_cuda` GPU 模式;外加一条 Atomiverse-CPU ↔ Atomiverse-CUDA 交叉校验。由于协议是随机性的 (stage 1 是带种子的 Langevin) 且 CUDA 是单精度,**所有跨引擎数值对比一律用集合均值的统计分块平均 (statistical block-average of ensemble means)** —— 遵循 `msu-pe-align-03`;不做对 LAMMPS 的逐步或比特级对比。
4. **集群 = PDC Dardel, SLURM。** 参考运行是操作员动作 (完整 ~360k 步协议约 24 小时级墙钟,PDC 作业上限 24h)。LAMMPS CPU + Atomiverse CPU 跑 `main` 分区;LAMMPS GPU 模式 + Atomiverse CUDA 跑 `gpugh` 分区 (Nvidia GH200 节点)。运行脚本通过 `sbatch`/`srun` 提交,需要 `-A <project>` 账户标志。

## 状态快照 (Status snapshot)

| Phase | Sub-item | Status |
|---|---|---|
| 0.1 | Ghost-aware CellList | ✅ landed |
| 0.2 | AtomMigration (atom data via MPI_Alltoallv) | ✅ landed |
| 0.3 | Tag-based bonded topology + cross-rank bonds + LJ exclusions | ✅ landed |
| 0.4 | ReduceEnergy 8-double allreduce | ✅ landed |
| 0.5 | NoseHoover global-KE allreduce | ✅ landed |
| 0.6 | VirialStress + Berendsen Barostat (multi-rank-correct) | ✅ landed |
| 0.7 | Deform uniaxial affine box update | ✅ landed |
| A.2 | C20H42 single-chain force-field check | ✅ landed |
| A.3 | LJ special-bonds exclusions + multi-chain box | ✅ landed |
| A.4 | Single-chain NVE conservation 1000 steps | ✅ landed |
| A.5 | Single-chain NoseHoover NVT stability | ✅ landed |
| B.2 | 8-rank PE bonded NVE conservation | ✅ landed |
| B.3 | 8-rank NoseHoover NVT, real chain crossing rank boundary | ✅ landed |
| End-to-end | 8-rank pipeline (bonded + LJ + halo + NVT + NPT + tension) | ✅ landed |
| C | MSU 10k-atom system construction (data file) | ✅ downloaded 2026-05-19 — `benchmarks/md/reference/polyethylene/msu_pe.data` |
| C.input | MSU paper-scale LAMMPS input deck | ✅ downloaded 2026-05-19 — `benchmarks/md/reference/polyethylene/lammps_msu_paper.in` |
| G.1 | Reference-fidelity defect: `fix momentum` (ZeroMomentum, CPU+CUDA) | 🔴 not started — hard fidelity defect |
| G.2 | Reference-fidelity defect: per-axis `Pxx/Pyy/Pzz` virial output | 🔴 not started — hard fidelity defect (run script currently FAKES it) |
| G.3 | Reference-fidelity defect: energy decomposition `E_pair/E_bond/E_angle/E_dihedral` | 🔴 not started — hard fidelity defect |
| G.4 | CUDA categorical setter API parity (`set_pair/set_bond/...`) | 🔴 not started — blocks CUDA benchmark wiring |
| G.5 | CUDA `NVELimit` integrator | 🔴 not started — blocks CUDA stage 1 |
| G.6 | CUDA `FixLangevin` thermostat | 🔴 not started — blocks CUDA stage 1 |
| G.7 | CUDA `NoseHooverNPT` (iso + aniso y/z) | 🔴 not started — blocks CUDA stages 2–5 |
| G.8 | CUDA `LJ` 1-2/1-3/1-4 special-bonds exclusions | 🔴 not started — silently corrupts bonded polymer on CUDA |
| G.9 | CUDA path in `polyethylene_tension.cpp` | 🔴 not started — depends on G.4–G.8 |
| H.1 | CPU `Thermo` → `AsyncComponentBase` | 🔴 not started |
| H.2 | CUDA driver hot-path desync (no per-phase-slot `defer_stream` sync) | 🔴 not started |
| H.3 | CUDA `VirialStress` deferrable when observed-only | 🔴 not started |
| H.4 | Async-correctness exit gate (step time independent of recorder freq) | 🔴 not started |
| I | CUDA single-GPU 5-stage correctness validation | 🔴 not started — blocked by G + H |
| D | Dual-backend LAMMPS comparison (2×2 + CPU↔CUDA) | 🟡 CPU arm sub-chained — `msu-pe-phase-d-*` + `msu-pe-align-03`; CUDA arm blocked by G/H/I |
| E.2/E.3 | CPU strong / weak scaling at 1–32 ranks | ⏳ external — needs cluster access for ranks > 8 |
| E.4 | CPU vs LAMMPS-CPU head-to-head | ⏳ depends on D |
| E.8 | Per-stage wall-clock timing across AV-CPU/AV-CUDA/LAMMPS-CPU/LAMMPS-GPU | 🔴 not started |
| F | Reproducibility artifacts — `lmp_cuda`, PDC `gpugh` sbatch, multi-engine validator | 🟡 CPU-only artifacts partially landed; CUDA/PDC artifacts not started |

What can be done **right now** without external/cluster tooling: all of Phase G (component implementation + unit tests), Phase H (async refactor + the async-correctness micro-benchmark gate), and the multi-engine extension of `validate_polyethylene.py`. What needs **external setup**: the operator cluster runs on PDC Dardel (`main` + `gpugh`), the LAMMPS `lmp_cuda` reference runs, and any quantitative LAMMPS comparison.

---

## Domain basis

The protocol physics is the immutable upstream deck `benchmarks/md/reference/polyethylene/lammps_msu_paper.in` (Tschopp 2010; M.A. Tschopp, M.F. Horstemeyer, M.K. Hossain et al., "Molecular dynamics simulations of deformation in amorphous polyethylene," *Polymer* 2010; CAVS ICME source <https://www.cavs.msstate.edu/icme/code/lammps/atomistic_research/polymers/amorphous.php>). System: 10 chains × 1000 united-atom CH₂ = 10 000 atoms, cubic box 80.0586 Å, mass 14.02 amu per UA bead (`benchmarks/md/reference/polyethylene/SOURCES.md`).

Force field — Dreiding-UA, verbatim from the deck: harmonic bond *k* = 350 kcal/mol/Å², *r*₀ = 1.53 Å; harmonic angle *k* = 60 kcal/mol/rad², θ₀ = 109.5°; multi/harmonic dihedral *A*₁..*A*₅ = 1.73, −4.49, 0.776, 6.99, 0.0 kcal/mol; LJ ε = 0.112 kcal/mol, σ = 4.01 Å, *r*_cut = 10.5 Å. Units: eV, Å, fs, amu, K (`core/units.h`); LAMMPS `units real` → Atomiverse eV via `1 kcal/mol = 1/23.060547830619026 eV`.

Five-stage protocol (every fix, ensemble, timestep and step count is byte-faithful ground truth):

| Stage | Ensemble / fix | T | P | Steps | dt |
|---|---|---|---|---|---|
| 1 | `nve/limit 0.05` + `langevin 500 500 10.0 904297` | 500 K | — | 10 000 | 1.0 fs |
| 2 | `npt ... iso` + `momentum 1 linear 1 1 1` | 500 K | 0 bar | 50 000 | 0.5 fs |
| 3 | `npt temp 500 100 ... iso` + `momentum` (cooling ramp) | 500→100 K | 0 bar | 50 000 | 0.5 fs |
| 4 | `npt ... iso` + `momentum` | 100 K | 0 bar | 50 000 | 0.5 fs |
| 5 | `npt ... y ... z` (aniso) + `deform 1 x erate 1e-5` | 100 K | 0 bar (y,z) | 171 800 | 1.0 fs |

Stage 1 is load-bearing in two ways: (a) `msu_pe.data` is a built-but-**unrelaxed** configuration with non-bonded overlaps (initial `E_pair ≈ 3.4×10⁸ kcal/mol`) — a displacement-capped `nve/limit` integrator is mandatory, a plain Velocity-Verlet integrator blows the system up; (b) `fix langevin` is **seeded** (`904297`), so stage 1 is a **stochastic** thermostat. Per `msu-pe-align-03` (block averaging — Flyvbjerg & Petersen 1989, *J. Chem. Phys.* 91, 461, doi:10.1063/1.457480), two engines running the same stochastic-thermostat protocol with different RNG streams and (for CUDA) different precision do **not** share a trajectory — they sample the same ensemble. **Every cross-engine value comparison in this spec is therefore a statistical block-average of ensemble means; no per-step or bit-level comparison against LAMMPS.** CUDA being single-precision makes this doubly binding.

Stress for the stress-strain curve: engineering strain ε = box_x/box_x0 − 1; σ_xx in GPa (`1 eV/Å³ = 160.21766 GPa`). Mechanical descriptors (per `msu-pe-phase-d-04`): elastic modulus *E* = dσ/dε over ε < 0.02, yield σ_y = max σ_xx over ε ∈ [0,0.1], hardening slope dσ/dε over ε ∈ [0.5,1.0].

### Strict reference-fidelity requirement

The benchmark protocol must be **byte-faithful in physics** to `lammps_msu_paper.in`: every fix, ensemble, integrator, timestep, step count and force-field coefficient. Three current fidelity defects are called out and fixed in Phase G:

- **`fix momentum` missing.** The deck applies `fix momentum 1 linear 1 1 1` (zero center-of-mass linear momentum every step) in stages 2/3/4. Atomiverse has no equivalent component on either backend. This is a physics deviation, not cosmetic.
- **Per-axis pressure faked.** The deck dumps `pxx pyy pzz` in stages 4 and 5; the stress-strain curve needs real σ_xx. The benchmark thermo CSV emits only isotropic `P_bar`, and `tools/run_msu_atomiverse.sh` currently *copies* isotropic P into the `Pxx/Pyy/Pzz` columns. That stand-in is wrong and must be replaced by a real per-axis virial path.
- **Energy decomposition missing.** The deck dumps `epair ebond eangle edihed` in stages 4 and 5; the user asks to compare energy. The benchmark emits only total `E_pot/E_kin` and must additionally emit `E_pair/E_bond/E_angle/E_dihedral`.

The per-stage timesteps and step counts in `polyethylene_tension.cpp` already match the reference (verified — stage 5 is dt = 1.0 fs, 171 800 steps).

---

## Dependency graph

```
Phase 0 ✅ ─→ Phase A ✅ ─→ Phase B ✅ ─→ Phase C ✅
                                            │
                ┌───────────────────────────┤
                ▼                           ▼
   Phase G (functional gaps)      Phase D-CPU arm
   CUDA components + fidelity     (msu-pe-phase-d-* chain
   defects — HARD-BLOCKS          + msu-pe-align-03) — CPU
   CUDA validation                quantitative comparison
                │
                ▼
   Phase H (async correctness)
   CPU Thermo async + CUDA
   driver hot-path desync
                │
                ▼
   Phase I (CUDA single-GPU 5-stage correctness)
                │
                ▼
   Phase D-CUDA arm  →  2×2 + CPU↔CUDA cross-check
                │
                ▼
   Phase E (per-stage timing + CPU scaling unchanged)
                │
                ▼
   Phase F (reproducibility artifacts — lmp_cuda, PDC gpugh)
```

Phase 0/A/B/C are landed. Phase G **hard-blocks** the entire CUDA validation arm — no CUDA stage runs without the missing CUDA components, and the fidelity defects (G.1–G.3) corrupt both backends until fixed. Phase H is independent of correctness and can proceed in parallel with G. Phase I requires both G and H. The Phase D-CUDA arm requires I. Phase D-CPU arm is already delegated to the `msu-pe-phase-d-*` sub-chain and `msu-pe-align-03` and is not re-specified here.

---

## Phase 0 — Engine blockers (✅ landed, unchanged)

Items 0.1–0.7 from `parallel-md-engine.md` — all ✅. Exit criterion (32-rank LJ + harmonic-bond run on 10⁵ atoms matches 1-rank to 1e-10 relative energy through 10 000 steps) met. Retained verbatim; see prior spec history.

## Phase A — Single-rank force-field correctness (✅ landed, unchanged)

A.2–A.5 all green on one CPU rank. Dreiding-UA force field migrated 2026-05-19 (`pe-ff-migrate-dreiding-ua`). Retained verbatim.

## Phase B — Multi-rank correctness (✅ landed, unchanged)

B.2/B.3 and the end-to-end 8-rank pipeline green. Retained verbatim.

## Phase C — MSU amorphous-PE system construction (✅ landed, unchanged)

`msu_pe.data` and `lammps_msu_paper.in` downloaded from the CAVS ICME wiki (`SOURCES.md`). System is 10 chains × 1000 UA = 10 000 atoms. Retained.

---

## Phase G — Functional-gap implementation (CUDA components + reference-fidelity defects)

Goal: implement every piece of missing functionality required before the CUDA backend can run the protocol at all, and fix the three reference-fidelity defects on both backends. **This phase hard-blocks the CUDA validation arm.** Each sub-item is a discrete component or output-path change; detailed implementable slicing of the larger CUDA components is delegated to a `msu-pe-cuda-components-*` chain at `/mol:impl` time — this umbrella states the gates.

### G.1 `ZeroMomentum` component (CPU + CUDA) — fidelity defect

Reference `fix momentum 1 linear 1 1 1` zeroes the center-of-mass linear momentum every step in stages 2/3/4. New `cpu::md::ZeroMomentum` and `cuda::md::ZeroMomentum` components (role `None` structural / Integrator-adjacent, no `Capability` write beyond velocities). Bound at the end-of-step phase slot; multi-rank-correct on CPU (COM momentum reduced across ranks before subtraction). Without it the CPU and CUDA runs both drift from the reference.

### G.2 Per-axis `Pxx/Pyy/Pzz` virial output — fidelity defect

The benchmark must emit a real per-axis pressure tensor, not isotropic P copied three ways. `VirialStress` (CPU `include/atomiverse/cpu/md/`, CUDA `include/atomiverse/cuda/observer/stress.cuh`) already computes the full virial tensor; the gap is the thermo/CSV path. The benchmark thermo writer must surface `Pxx/Pyy/Pzz` columns from the diagonal of the virial tensor, and `tools/run_msu_atomiverse.sh` must stop faking those columns.

### G.3 Energy decomposition output — fidelity defect

The benchmark must emit `E_pair/E_bond/E_angle/E_dihedral` columns (the reference deck dumps `epair ebond eangle edihed`). Each evaluator component already computes its own energy contribution; the gap is a per-term accumulation surfaced into the thermo CSV.

### G.4 CUDA categorical setter API parity

`cuda::md::MD` (`include/atomiverse/cuda/md/md.cuh`) currently exposes only raw `add()` / `on(phase).add()`. CPU `cpu::md::MD` has `set_pair/set_bond/set_angle/set_dihedral/set_neighbor/set_integrator/set_thermostat/set_barostat/set_observer/set_record/set_thermo`. CUDA must reach parity with the same categorical surface (and the same auto-injection behavior — halo/VirialStress auto-inject — where applicable to single-GPU). This is needed so the CUDA path in the benchmark mirrors the CPU path.

### G.5 CUDA `NVELimit` integrator

New `cuda::md::NVELimit` (`include/atomiverse/cuda/md/nve_limit.cuh` (new) + `src/cuda/md/nve_limit.cu` (new)). Displacement-capped integrator (`nve/limit 0.05`) — load-bearing for stage 1 on the unrelaxed `msu_pe.data`. CPU equivalent at `include/atomiverse/cpu/md/nve_limit.h`.

### G.6 CUDA `FixLangevin` thermostat

New `cuda::md::FixLangevin` (`include/atomiverse/cuda/md/fix_langevin.cuh` (new) + `src/cuda/md/fix_langevin.cu` (new)). Stage-1 Langevin thermostat: per-atom RNG state, friction drag, random kick, temperature ramp. Thermostat role, also Force-capable. CPU equivalent at `include/atomiverse/cpu/md/fix_langevin.h`.

### G.7 CUDA `NoseHooverNPT` integrator

New `cuda::md::NoseHooverNPT` (`include/atomiverse/cuda/md/nose_hoover_npt.cuh` (new) + `src/cuda/md/nose_hoover_npt.cu` (new)). Integrated NPT (integrator + thermostat + barostat) supporting **both** isotropic coupling (stages 2–4) and **anisotropic y/z-only** coupling (stage 5), plus the 500→100 K temperature ramp. CUDA today has only separate `cuda::md::NoseHoover` + `cuda::md::Barostat`; CPU equivalent at `include/atomiverse/cpu/md/nose_hoover_npt.h`.

### G.8 CUDA `LJ` special-bonds exclusions

CUDA `LJ` / `pair_lj_cut` (`include/atomiverse/cuda/potentials/pair_lj_cut.cuh`, `include/atomiverse/cuda/potentials/lj_kernel.cuh`, `src/cuda/potentials/pair_lj_cut.cu`) has no 1-2/1-3/1-4 special-bonds exclusion support — the kernel processes every neighbor pair unconditionally, silently corrupting any bonded polymer on CUDA. Add a `set_special_bonds` API and exclusion masking in the kernel, matching CPU `pair_lj_cut.h`.

### G.9 CUDA path in `polyethylene_tension.cpp`

`benchmarks/md/polyethylene_tension.cpp` is CPU-only — five per-stage `cpu::md::MD` instances, no CUDA branch. Add an `#ifdef`-guarded CUDA path that builds the equivalent five-stage `cuda::md::MD` pipeline via the G.4 categorical API, reading the same `msu_pe.data`, writing the same CSV schema (incl. G.2 per-axis P and G.3 energy decomposition). Selectable at runtime (engine flag).

**Exit gate for Phase G**: (a) `ZeroMomentum` (CPU+CUDA), `cuda::md::NVELimit`, `cuda::md::FixLangevin`, `cuda::md::NoseHooverNPT` (iso + aniso), CUDA `LJ` special-bonds exclusions all have green unit tests in `tests/cuda/md/` / `tests/cpu/md/`; (b) CUDA categorical setter API compiles and a smoke CUDA MD run via `set_*` completes; (c) the benchmark emits `Pxx/Pyy/Pzz` and `E_pair/E_bond/E_angle/E_dihedral` real columns on both backends and `run_msu_atomiverse.sh` no longer fakes per-axis P; (d) `polyethylene_tension.cpp` builds and runs a short (≤1000-step) protocol on CUDA without NaN.

---

## Phase H — Async correctness (no recorder ever stalls the loop)

Goal: make every non-hot-path recorder / thermo / observer genuinely asynchronous so the integration loop — and on CUDA the GPU compute stream — never stalls on an analyzer/recorder.

### H.1 CPU `Thermo` → `AsyncComponentBase`

CPU `Thermo` (`include/atomiverse/cpu/recorder/thermo.h`, `src/cpu/recorder/thermo.cpp`) runs synchronously inside `launch()` on the integration thread — it is not an async observer. The CPU async machinery already exists: `AsyncComponentBase` (`include/atomiverse/cpu/async_component_base.h`) + `AnalysisExecutor` (`include/atomiverse/cpu/analysis_executor.h`, triple-buffered snapshot ring + worker pool) — CPU `RDF` already uses it. Move `Thermo` onto `AsyncComponentBase` so its reduction + I/O run off the integration thread.

### H.2 CUDA driver hot-path desync

CUDA `AsyncBridge` (`include/atomiverse/cuda/async_bridge.cuh`) is sound — a sidecar `defer_stream` + CUDA event lets recorders' D2H + host callback run off the compute stream. But the CUDA driver `src/cuda/driver.cu` `Backend::step()` calls `cudaStreamSynchronize(stream)` **and** `cudaStreamSynchronize(defer_stream)` after **every phase slot, every step**, which defeats the deferral: the host blocks on the recorder's D2H+callback before the next phase slot, so recorders do stall the loop and the GPU ends up waiting on the analyzer/recorder. The fix: deferred work must lag behind the compute stream across steps — snapshot device buffers into recorder-owned buffers so compute can overwrite the live frame; remove the per-phase-slot `defer_stream` sync from the hot path; a bounded backlog drains at run teardown. The `stream` (compute) sync stays where genuinely required for correctness; only the recorder/observer deferral path is desynced.

### H.3 CUDA `VirialStress` deferrable when observed-only

CUDA `VirialStress` (`include/atomiverse/cuda/observer/stress.cuh`, `src/cuda/observer/stress.cu`) launches its kernel on the main compute stream — correct when a barostat consumes the stress that same step (genuinely hot-path), but when stress is only being observed/recorded it must be deferrable onto the `defer_stream` like other observers. Make `VirialStress` aware of whether a same-step consumer (barostat) reads its output: hot-path when consumed, deferred when observed-only.

### H.4 Async-correctness exit gate

A micro-benchmark in `benchmarks/md/` measures steps/s on the CUDA backend with recorders/observers off versus at high recording frequency. The gate: step wall-time must be independent of recorder/observer frequency within a small tolerance (target: high-frequency steps/s ≥ 95% of recorders-off steps/s) — proving the GPU never stalls on a recorder. A CPU counterpart confirms the same for the now-async `Thermo`.

**Exit gate for Phase H**: H.1, H.2, H.3 implemented with unit/integration tests green; H.4 micro-benchmark shows CUDA steps/s with recorders at high frequency within 5% of recorders-off, and CPU `Thermo`-on within 5% of `Thermo`-off.

---

## Phase I — CUDA single-GPU 5-stage correctness validation

Goal: mirror the CPU Phase A–C structure for the CUDA backend on a single GPU — prove the CUDA force field is correct, the full 5-stage protocol completes, and CUDA agrees with both the Atomiverse CPU run and LAMMPS GPU mode. Blocked by G (components) and H (async). All cross-engine comparisons are statistical block-averages per `msu-pe-align-03`.

### I.1 CUDA force-field correctness

New `tests/cuda/md/test_pe_cuda_ff.cpp` (new): evaluate `E_pair/E_bond/E_angle/E_dihedral` for the C₂₀H₄₂ perturbed-geometry cross-check (the same geometry as `tests/cpu/md/test_pe_energy_xcheck.cpp`, with the closed-form perturbation in `SOURCES.md`) on CUDA. Because CUDA is single-precision the tolerance is statistical/loose, not the CPU 1e-6 eV bound: each energy term agrees with the CPU double-precision value to a relative tolerance consistent with float accumulation over 10 000 atoms (target ≤ 1e-3 relative). Confirms G.5–G.8 components produce physically correct forces, including special-bonds exclusions.

### I.2 CUDA full 5-stage protocol completes

CUDA path of `polyethylene_tension.cpp` runs the full 5-stage protocol on `msu_pe.data` on a single GH200 GPU to completion without NaN, writing the canonical CSV schema (thermo with per-axis P + energy decomposition, stress-strain). Operator action on the PDC `gpugh` partition. Density after stage 4 within 5% of the MSU reported amorphous-PE density at 100 K (~0.94 g/cm³).

### I.3 CUDA ↔ CPU cross-check

`tools/validate_polyethylene.py` compares the Atomiverse-CUDA and Atomiverse-CPU runs on block-averaged ensemble means: ⟨T⟩, ⟨E_pot⟩ and the energy decomposition over stage 4, density, and the three stress-strain descriptors. Tolerance is the statistical bar widened for single vs double precision (target ~2–3% relative on means, ±5% on stress-strain descriptors). This isolates precision/RNG effects from genuine engine bugs before the LAMMPS comparison.

### I.4 CUDA ↔ LAMMPS-GPU statistical comparison

`tools/validate_polyethylene.py` compares the Atomiverse-CUDA run against the `lmp_cuda` GPU-mode reference run on block-averaged ensemble means — the GPU-mode arm of the 2×2. Same statistical methodology and tolerances as the CPU↔LAMMPS arm (`msu-pe-align-03` / `msu-pe-phase-d-*`).

**Exit gate for Phase I**: I.1 green; I.2 completes the 5-stage CUDA run with density within 5% of reference; I.3 CUDA↔CPU statistical comparison passes; I.4 CUDA↔LAMMPS-GPU statistical comparison passes.

---

## Phase D — Dual-backend LAMMPS comparison (2×2 + CPU↔CUDA cross-check)

Goal: demonstrate Atomiverse and LAMMPS, starting from the same `msu_pe.data` and running the same protocol, produce statistically matching observables on **both** backends. Phase D is now dual-backend and four-config.

### D-CPU arm — already sub-chained (referenced, not re-specified)

The CPU arm of Phase D is delegated to the `msu-pe-phase-d-*` chain and `msu-pe-align-03`:
- `msu-pe-phase-d-01-reference-runs` — produces the LAMMPS + Atomiverse-CPU reference CSVs.
- `msu-pe-align-03-statistical-thermo` (supersedes `msu-pe-phase-d-02`) — block-average statistical thermo parity, `validate_polyethylene.py --check thermo`.
- `msu-pe-phase-d-03-equilibrium` — equilibrium-state observables (ρ, ⟨E_pot⟩, RDF peak, ⟨R²⟩), `--check equilibrium`.
- `msu-pe-phase-d-04-stress-strain` — stress-strain descriptors + Phase D exit gate, `--check stress-strain` / `--check all`.

This umbrella does **not** redo those sub-specs. The CUDA arm below stays consistent with them: same statistical block-average methodology, same `validate_polyethylene.py` surface, same exit-gate shape.

### D.5 The 2×2 comparison matrix

Four engine configs, three comparison tracks — all statistical block-average of ensemble means:

| Track | Engine A | Engine B | Methodology |
|---|---|---|---|
| CPU reference | Atomiverse-CPU | `lmp_cuda` CPU-mode | block-average — `msu-pe-phase-d-*` chain |
| GPU reference | Atomiverse-CUDA | `lmp_cuda` GPU-mode | block-average — Phase I.4 |
| CPU↔CUDA cross-check | Atomiverse-CPU | Atomiverse-CUDA | block-average — Phase I.3 |

Compared observables in every track: energy (incl. the `E_pair/E_bond/E_angle/E_dihedral` decomposition from G.3), density ρ, equilibrium observables (RDF peak, ⟨R²⟩), and the stress-strain curve descriptors (E, σ_y, hardening slope). `tools/validate_polyethylene.py` becomes multi-engine: it must accept the four engine CSVs and run all three tracks, keeping the `msu-pe-align-03` statistical methodology and the `--check {thermo,equilibrium,stress-strain,all}` surface, extended with energy-decomposition columns and a track selector.

**Exit gate for Phase D**: the CPU arm passes per the `msu-pe-phase-d-*` chain exit gates; the GPU-reference and CPU↔CUDA tracks pass per Phase I.3/I.4; `validate_polyethylene.py` runs all three tracks and exits 0.

---

## Phase E — Performance benchmarks (CPU scaling unchanged + per-stage timing)

The CPU/MPI scaling work (E.1–E.7: single-rank breakdown, strong scaling, weak scaling, head-to-head vs LAMMPS-CPU, profile decomposition, hybrid MPI×OpenMP, multi-node) stays **CPU-only and unchanged** — see prior spec body. Multi-GPU / CUDA+MPI scaling is explicitly out of scope. One new sub-item is added:

### E.8 Per-stage wall-clock timing across four engine configs

Add per-stage wall-clock timing to `benchmarks/md/polyethylene_tension.cpp` (a timer bracketing each of the five stages) and emit a per-stage timing CSV. Compare per-stage wall time across the four configs — Atomiverse-CPU, Atomiverse-CUDA, `lmp_cuda` CPU-mode, `lmp_cuda` GPU-mode — via `tools/validate_polyethylene.py` (a new `--check timing` track). This is a timing/performance comparison, not a physics observable; it is descriptive (report the four-way per-stage breakdown) rather than a hard pass/fail tolerance gate, except that Atomiverse-CUDA per-stage time must be finite and the run must complete.

**Exit gate for Phase E**: E.1–E.5 pass their CPU efficiency targets (unchanged); E.8 emits a four-config per-stage timing table and `validate_polyethylene.py --check timing` produces it without error. E.6/E.7 remain stretch.

---

## Phase F — Reproducibility artifacts (dual-backend, PDC Dardel)

Everything an external reviewer needs to reproduce the dual-backend results on PDC Dardel.

### F.1 Input files under version control

Unchanged set under `benchmarks/md/reference/polyethylene/` (`force_field.yaml`, `msu_pe.data`, `lammps_msu_paper.in`, the reference CSVs, `SOURCES.md`, `README.md`), plus the GPU-mode and CUDA reference CSVs (engine-suffixed per the `SOURCES.md` Phase D.1 filename scheme, extended with a GPU/CUDA suffix).

### F.2 Scripts — `lmp_cuda` + PDC `gpugh`

- `tools/run_msu_lammps.sh` — extended to launch `lmp_cuda` in **both** CPU mode and GPU mode (not the old CPU/MPI-only path); PDC `sbatch` submission, `-A <project>` account flag, `main` partition for CPU mode and `gpugh` partition for GPU mode.
- `tools/run_msu_atomiverse.sh` — extended with a CUDA / GPU-mode launch path on the `gpugh` partition; stops faking per-axis P (G.2).
- `tools/run_msu_pdc.sh` (new) — PDC Dardel SLURM submission wrapper: `sbatch`/`srun`, `-A` account flag (concrete values: `naiss2025-1-25` for the `main` CPU partition, `naiss2025-1-25-gpu` for the `gpugh` GPU partition), partition selection (`main` vs `gpugh`), 24h wall-time cap awareness. The flag stays required with no default in the script CLI; the concrete account belongs in docs/reproduction commands, not hardcoded.
- `tools/validate_polyethylene.py` — the multi-engine 2×2 validator (D.5 / E.8): three comparison tracks, energy-decomposition columns, per-stage timing, statistical block-average methodology retained from `msu-pe-align-03`.
- `tools/plot_scaling.py`, `tools/plot_profile.py` — unchanged.

### F.3 README

`benchmarks/md/reference/polyethylene/README.md` updated: the `lmp_cuda` build/usage (CPU mode + GPU mode), the PDC Dardel `gpugh` / `main` partition `sbatch` procedure with the `-A` account flag, the four-config 2×2 reproduction commands, and the dual-backend expected-results table. Reference: <https://support.pdc.kth.se/doc/run_jobs/job_scheduling/>.

**Exit gate for Phase F**: a fresh clone plus the documented operator commands reproduce all four engine runs on PDC Dardel and `tools/validate_polyethylene.py` runs the three 2×2 tracks + the timing track green.

---

## Files to create or modify

### New files

- `include/atomiverse/cpu/md/zero_momentum.h` (new) — Phase G.1
- `src/cpu/md/zero_momentum.cpp` (new) — Phase G.1
- `include/atomiverse/cuda/md/zero_momentum.cuh` (new) — Phase G.1
- `src/cuda/md/zero_momentum.cu` (new) — Phase G.1
- `include/atomiverse/cuda/md/nve_limit.cuh` (new) — Phase G.5
- `src/cuda/md/nve_limit.cu` (new) — Phase G.5
- `include/atomiverse/cuda/md/fix_langevin.cuh` (new) — Phase G.6
- `src/cuda/md/fix_langevin.cu` (new) — Phase G.6
- `include/atomiverse/cuda/md/nose_hoover_npt.cuh` (new) — Phase G.7
- `src/cuda/md/nose_hoover_npt.cu` (new) — Phase G.7
- `tests/cpu/md/test_zero_momentum.cpp` (new) — Phase G.1
- `tests/cuda/md/test_cuda_md_components.cpp` (new) — Phase G.5/G.6/G.7/G.8 unit tests (ZeroMomentum-CUDA, NVELimit, FixLangevin, NoseHooverNPT, LJ exclusions)
- `tests/cuda/md/test_pe_cuda_ff.cpp` (new) — Phase I.1 CUDA force-field correctness
- `tools/run_msu_pdc.sh` (new) — Phase F.2 PDC Dardel SLURM wrapper

### Modified files

- `include/atomiverse/cuda/md/md.cuh` — Phase G.4 categorical setter API
- `src/cuda/driver.cu` — Phase G.4 (setter wiring) + Phase H.2 (hot-path `defer_stream` desync)
- `include/atomiverse/cuda/potentials/pair_lj_cut.cuh` — Phase G.8 `set_special_bonds`
- `include/atomiverse/cuda/potentials/lj_kernel.cuh` — Phase G.8 exclusion masking in kernel
- `src/cuda/potentials/pair_lj_cut.cu` — Phase G.8 exclusion launch wiring
- `benchmarks/md/polyethylene_tension.cpp` — Phase G.2/G.3 (per-axis P + energy-decomposition CSV columns, per-stage timer for E.8), Phase G.9 (CUDA path)
- `include/atomiverse/cpu/recorder/thermo.h` — Phase H.1 `AsyncComponentBase`
- `src/cpu/recorder/thermo.cpp` — Phase H.1 async migration
- `include/atomiverse/cuda/recorder/thermo.cuh` — Phase G.2/G.3 per-axis P + energy-decomposition columns on the CUDA thermo recorder
- `include/atomiverse/cuda/observer/stress.cuh` — Phase H.3 deferrable-when-observed-only
- `src/cuda/observer/stress.cu` — Phase H.3 defer-stream path
- `tools/validate_polyethylene.py` — Phase D.5 multi-engine 2×2 tracks + energy decomposition + Phase E.8 `--check timing`
- `tools/run_msu_atomiverse.sh` — Phase F.2 CUDA/gpugh launch path; stop faking per-axis P
- `tools/run_msu_lammps.sh` — Phase F.2 `lmp_cuda` CPU-mode + GPU-mode launch
- `tools/test_validate_polyethylene.py` — Phase D.5/E.8 test coverage for the new tracks
- `benchmarks/md/reference/polyethylene/README.md` — Phase F.3 dual-backend + PDC reproduction
- `benchmarks/md/reference/polyethylene/SOURCES.md` — Phase F.1 GPU/CUDA reference-CSV filename scheme
- `CMakeLists.txt` — register new CUDA components, new tests, the async-correctness micro-benchmark
- `.claude/specs/INDEX.md` — update this spec's umbrella entry to dual-backend
- `.claude/NOTES.md` — record the dual-backend MSU validation umbrella + CUDA-component implementation

## Tasks

> This umbrella spec delegates detailed implementable slicing of the larger CUDA components to a `msu-pe-cuda-components-*` chain at `/mol:impl` time. The tasks below are the umbrella-level tracker; each maps to a phase exit gate.

- [ ] Write failing tests for ZeroMomentum (CPU+CUDA), CUDA NVELimit/FixLangevin/NoseHooverNPT, CUDA LJ special-bonds exclusions (tests/cpu/md/test_zero_momentum.cpp, tests/cuda/md/test_cuda_md_components.cpp)
- [ ] Implement ZeroMomentum, cuda::md::NVELimit, cuda::md::FixLangevin, cuda::md::NoseHooverNPT, and CUDA LJ exclusions per Phase G.1/G.5/G.6/G.7/G.8
- [ ] Implement CUDA categorical setter API in include/atomiverse/cuda/md/md.cuh and src/cuda/driver.cu per Phase G.4
- [ ] Implement per-axis Pxx/Pyy/Pzz and E_pair/E_bond/E_angle/E_dihedral CSV columns and the CUDA path in benchmarks/md/polyethylene_tension.cpp per Phase G.2/G.3/G.9
- [ ] Implement CPU Thermo on AsyncComponentBase and CUDA driver hot-path defer_stream desync per Phase H.1/H.2/H.3
- [ ] Verify the async-correctness gate: CUDA steps/s at high recorder frequency within 5% of recorders-off (Phase H.4 micro-benchmark)
- [ ] Verify CUDA force-field correctness against CPU double-precision values in tests/cuda/md/test_pe_cuda_ff.cpp per Phase I.1
- [ ] Implement the multi-engine 2×2 tracks and --check timing in tools/validate_polyethylene.py per Phase D.5/E.8
- [ ] Implement lmp_cuda CPU/GPU-mode launch and PDC gpugh sbatch in tools/run_msu_lammps.sh, tools/run_msu_atomiverse.sh, tools/run_msu_pdc.sh per Phase F.2
- [ ] Add doxygen docstrings with units for all new CUDA components and update README.md per Phase F.3
- [ ] Verify the CUDA↔CPU and CUDA↔LAMMPS-GPU statistical block-average comparisons pass per Phase I.3/I.4 (operator cluster run)
- [ ] Run full check + test suite

## Testing strategy

- **Happy path** — `tests/cuda/md/test_cuda_md_components.cpp` and `tests/cpu/md/test_zero_momentum.cpp` exercise each new component on a small synthetic system: `ZeroMomentum` drives COM momentum to zero; `NVELimit` caps displacement at the configured bound; `FixLangevin` drives temperature to the setpoint with a fixed seed; `NoseHooverNPT` relaxes a box to target pressure in both iso and aniso y/z modes; CUDA `LJ` with special-bonds excludes 1-2/1-3/1-4 pairs. The CUDA benchmark path runs a short (≤1000-step) protocol without NaN.
- **Edge cases** — `NVELimit` on the overlapped `msu_pe.data` (initial E_pair ≈ 3.4×10⁸ kcal/mol) does not blow up; CUDA `LJ` with an empty special-bonds list behaves identically to the un-masked kernel; `validate_polyethylene.py` with a missing engine CSV exits non-zero with a clear error naming the absent track; the `--check timing` track on a CSV missing a stage column fails gracefully, not with a crash.
- **Domain validation** — CUDA force field: `test_pe_cuda_ff.cpp` confirms each of `E_pair/E_bond/E_angle/E_dihedral` agrees with the CPU double-precision reference to ≤1e-3 relative on the C₂₀H₄₂ perturbed geometry. Stochastic-protocol parity: per `msu-pe-align-03`, the CUDA↔CPU and CUDA↔LAMMPS-GPU comparisons use block-averaged ensemble means (⟨T⟩, ⟨E_pot⟩, energy decomposition, density) over the equilibrated stages, never per-step; CUDA single precision widens the means tolerance to ~2–3% relative. Stress-strain descriptors (E, σ_y, hardening slope) agree to ±5%. Reference fidelity: an assertion confirms the benchmark applies a `fix momentum` equivalent in stages 2/3/4 and emits real (non-faked) `Pxx/Pyy/Pzz`.
- **Async correctness** — the H.4 micro-benchmark measures steps/s with recorders/observers off vs at high frequency; the gate requires high-frequency steps/s ≥ 95% of recorders-off steps/s on CUDA, proving the GPU never stalls on a recorder; the CPU counterpart confirms the same for the now-async `Thermo`.

## Out of scope

- Multi-GPU and CUDA+MPI scaling — explicitly out of scope per fixed scope decision 1; CUDA validation is single-GPU correctness only. The CPU/MPI Phase E (E.1–E.7 strong/weak scaling) stays CPU-only and unchanged. (Considered adding a CUDA scaling phase; rejected — the user fixed CUDA scope to single-GPU correctness.)
- Re-specifying the CPU Phase D sub-chain — `msu-pe-phase-d-01..04` and `msu-pe-align-03` are referenced and kept consistent with, not rewritten here. (Considered folding them back into the umbrella; rejected — the umbrella delegates detailed slicing to per-phase chains by design.)
- Per-step or bit-level comparison against LAMMPS — physically impossible for this stochastic, single-precision-on-CUDA protocol; all cross-engine comparison is statistical block-average per `msu-pe-align-03`.
- Matching MSU paper numbers directly — LAMMPS under identical inputs is the reproducible reference; MSU paper-number matching remains a downstream addendum, unchanged from the prior spec rationale.
- Changing the protocol physics — `lammps_msu_paper.in` is immutable ground truth; the benchmark is corrected to match it, never the reverse.

## Rationale

- **Why supersede in place instead of splitting?** The user explicitly asked to update this one umbrella file. The umbrella already delegates detailed implementable slicing to per-phase chains (`msu-pe-phase-d-*`); the new CUDA-component work likewise delegates to a `msu-pe-cuda-components-*` chain at impl time. Keeping one umbrella preserves the single roadmap view.
- **Why a dedicated functional-gap phase (G) before any CUDA validation?** The user said missing functionality must be implemented first. The CUDA backend literally cannot run the protocol — it has no `NVELimit`, `FixLangevin`, `NoseHooverNPT`, no special-bonds exclusions (which silently corrupt the polymer), and no categorical setter API. Phase G hard-blocks the CUDA arm; making it explicit prevents starting validation on a non-functional backend.
- **Why call out `fix momentum`, per-axis P, and energy decomposition as fidelity defects?** The benchmark does not yet strictly match `lammps_msu_paper.in`: it omits `fix momentum` (a physics deviation), fakes `Pxx/Pyy/Pzz` by copying isotropic P (wrong for the stress-strain curve), and omits the energy decomposition the user wants to compare. These are correctness defects, not enhancements.
- **Why statistical block-average for every CUDA comparison?** `msu-pe-align-03` established that the MSU protocol is stochastic (seeded stage-1 Langevin). CUDA additionally runs single precision. Two engines therefore never share a trajectory — they sample the same ensemble. The only physically meaningful cross-engine statement is agreement of block-averaged ensemble means.
- **Why a dedicated async-correctness phase (H)?** The user's hard requirement is that the GPU never waits on a recorder. The CUDA `AsyncBridge` exists but is defeated by per-phase-slot `cudaStreamSynchronize(defer_stream)` in the driver hot path. An explicit phase with a measurable exit gate (step time independent of recorder frequency) turns "make it async" into a verifiable contract.
- **Why PDC Dardel as the cluster?** Fixed scope decision 4: the reference runs are ~24h-class operator actions; PDC Dardel's `main` (CPU) and `gpugh` (GH200 GPU) partitions map cleanly onto the 2×2 engine matrix, and the 24h job cap is why the runs are operator actions, not CI steps.

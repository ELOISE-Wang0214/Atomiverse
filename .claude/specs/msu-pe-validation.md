---
title: MSU amorphous-PE tension — dual-backend (CPU + CUDA) correctness and performance validation vs LAMMPS
status: in-progress
created: 2026-05-22
last-updated: 2026-05-24 (cleanup pass 2 — corrected the status snapshot to reflect actual reality: G.1–G.9 + H.1–H.3 are ✅ landed via the `msu-pe-cuda-components-01..07` chain; E.8 per-stage timing landed; Phase F tooling all present. Collapsed Phase G + H section bodies into ✅ summaries pointing at chain specs. Remaining work is operator-action SLURM runs. Earlier 2026-05-24 cleanup pass removed Phase 0/A/B/C placeholder sections. 2026-05-22 SUPERSEDE/REFINE: extended the umbrella to dual-backend with single-GPU CUDA correctness validation against LAMMPS GPU mode, a functional-gap phase G, an async-correctness phase H, and a 2×2 four-config comparison; Phase D sub-chain `msu-pe-phase-d-*` and `msu-pe-align-03` referenced, not re-specified.)
supersedes: msu-pe-validation (in place)
---

# MSU amorphous-PE validation (dual-backend)

## 综述 (Summary)

本规范是 `msu-pe-validation` 的就地 SUPERSEDE/REFINE：把原本只针对 CPU/MPI 引擎的 MSU 非晶聚乙烯单轴拉伸验证，扩展为 **CPU 与 CUDA 双后端** 对 LAMMPS 的验证。用户的最终诉求是：MSU 非晶 PE 拉伸基准要同时验证 Atomiverse 的 CPU 后端和 CUDA 后端,并在四种引擎配置上对比能量、密度、每阶段墙钟时间和应力-应变曲线;凡是缺失的功能必须**先实现**;所有非热路径的 recorder / thermo / observer 必须改为异步,尤其在 CUDA 上 GPU 绝不能因为 analyzer/recorder 而停等。

本规范保持 **伞型路线图 (umbrella roadmap)** 形态：它持有阶段状态快照、依赖图、各阶段的退出门 (exit gate)。CPU 侧的 Phase 0 / A / B / C 与功能补齐阶段 G (G.1–G.9)、异步正确性阶段 H 的 H.1–H.3、性能阶段 E.8、复现工具阶段 F 全部 ✅ landed。剩余 **operator-action**：Phase D 的 CPU↔LAMMPS 重跑校验 (recent fidelity 修复后)、Phase I 的 CUDA 单卡 5 阶段验证 + 与 LAMMPS-GPU 对比 (delegated to `msu-pe-cuda-components-08`)。用户可见的成果：一次基准运行同时给出 CPU 与 CUDA 对 LAMMPS 的统计学一致性结论，以及四引擎的逐阶段性能对比，且 GPU 永不为 recorder 停等。

## 已确认的范围决定 (fixed scope decisions)

1. **CUDA 范围 = 单卡正确性验证。** CUDA 在单张 GPU 上跑完整 5 阶段协议;结果与 LAMMPS GPU 模式、以及 Atomiverse CPU 运行对比。多 GPU / CUDA+MPI 扩展性 **不在范围内** —— 现有的 CPU/MPI Phase E (强/弱扩展性, E.1–E.7) 保持 CPU-only 不变。
2. **"每阶段速度" = 每阶段墙钟时间。** 给基准加每阶段墙钟计时,在四种引擎配置间对比。这是计时/性能对比,不是原子速度的物理观测量。
3. **对比是同类相比的 2×2。** 用户构建了单个 LAMMPS 二进制 `lmp_cuda`,可跑 CPU 模式和 GPU 模式两种模式。两条参考轨道:Atomiverse-CPU ↔ `lmp_cuda` CPU 模式;Atomiverse-CUDA ↔ `lmp_cuda` GPU 模式;外加一条 Atomiverse-CPU ↔ Atomiverse-CUDA 交叉校验。由于协议是随机性的 (stage 1 是带种子的 Langevin) 且 CUDA 是单精度,**所有跨引擎数值对比一律用集合均值的统计分块平均 (statistical block-average of ensemble means)** —— 遵循 `msu-pe-align-03`;不做对 LAMMPS 的逐步或比特级对比。
4. **集群 = PDC Dardel, SLURM。** 参考运行是操作员动作 (完整 ~360k 步协议约 24 小时级墙钟,PDC 作业上限 24h)。LAMMPS CPU + Atomiverse CPU 跑 `main` 分区;LAMMPS GPU 模式 + Atomiverse CUDA 跑 `gpugh` 分区 (Nvidia GH200 节点)。运行脚本通过 `sbatch`/`srun` 提交,需要 `-A <project>` 账户标志。

## 状态快照 (Status snapshot)

| Phase | Sub-item | Status |
|---|---|---|
| 0 / A / B / C | Engine blockers, single-/multi-rank correctness, MSU system construction | ✅ landed |
| G.1–G.9 | All functional gaps — fidelity defects (`fix momentum`, per-axis `Pxx/Pyy/Pzz`, energy decomposition) + CUDA components (categorical API, `NVELimit`, `FixLangevin`, `NoseHooverNPT` iso+aniso, LJ special-bonds) + CUDA benchmark path (`--backend {cpu\|cuda}`) | ✅ landed — chain `msu-pe-cuda-components-01..06`. CPU `ZeroLinearMomentum` (alias `ZeroMomentum`) in `include/atomiverse/cpu/md/init.h`; CUDA `cuda::md::ZeroMomentum` in `init.cuh`; benchmark stages 2/3/4 wired; thermo emits real per-axis P + 4-term energy decomposition. |
| H.1 / H.2 / H.3 | CPU `Thermo` async + CUDA driver hot-path desync + `VirialStress` deferrable when observed-only | ✅ landed — chain `msu-pe-cuda-components-07` (`cuda::AsyncRecordBase` + zero-`cudaStreamSynchronize` hot-path guard). |
| H.4 | Async-correctness exit gate (step time independent of recorder freq) | 🟡 9/10 ACs verified — 25×-frequency wall-clock < 5% gate passes; ac-007 `compute-sanitizer memcheck` deferred (tool absent from `~/atv-build-env`). |
| I | CUDA single-GPU 5-stage correctness validation | 🟡 tooling code-complete (chain `msu-pe-cuda-components-08`) — ac-007..010 pending operator-action `gpugh` SLURM runs. |
| D | Dual-backend LAMMPS comparison (2×2 + CPU↔CUDA) | 🟡 CPU arm: `msu-pe-phase-d-*` + `msu-pe-align-03` code-complete; D.2 / D.3 / D.4 re-verify pending the operator reference run. CUDA arm = Phase I. |
| E.2 / E.3 | CPU strong / weak scaling at 1–32 ranks | ⏳ external — needs cluster access for ranks > 8. |
| E.4 | CPU vs LAMMPS-CPU head-to-head | ⏳ depends on D operator runs. |
| E.8 | Per-stage wall-clock timing across AV-CPU/AV-CUDA/LAMMPS-CPU/LAMMPS-GPU | ✅ landed — per-stage timer in `polyethylene_tension.cpp` (commit `4c3e109`). |
| F | Reproducibility artifacts — `lmp_cuda`, PDC `gpugh` sbatch, multi-engine validator | ✅ tooling landed — `tools/run_msu_atomiverse.sh --backend {cpu\|cuda}`, `tools/run_msu_lammps.sh --mode {cpu\|gpu}`, `tools/run_msu_pdc.sh`, `tools/compare_msu.py`, `tools/plot_msu_compare.py`. The operator reproduction runs themselves remain external. |

**Effective state.** All code-side functional gaps (G) and async-correctness work (H.1–H.3) are landed. The remaining work is **operator-action**: submit the four reference runs (LAMMPS-CPU + Atomiverse-CPU on `main`; LAMMPS-GPU + Atomiverse-CUDA on `gpugh`) on PDC Dardel and execute `tools/compare_msu.py` + `tools/plot_msu_compare.py` against the resulting CSVs — i.e., the Phase D (CPU arm) re-verify pass after the recent benchmark-fidelity fixes and the Phase I CUDA validation.

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
Phase 0/A/B/C ✅ (engine blockers, FF correctness, system construction)
                │
                ├─→ Phase G ✅ (functional gaps — CUDA components + fidelity defects;
                │              chain `msu-pe-cuda-components-01..06`)
                │
                ├─→ Phase H.1–H.3 ✅ (async correctness — chain `msu-pe-cuda-components-07`)
                │
                ├─→ Phase D-CPU arm ✅ tooling (msu-pe-phase-d-* + msu-pe-align-03;
                │                              operator-run gated)
                │
                ▼
   Phase I ✅ tooling (CUDA single-GPU 5-stage; chain `msu-pe-cuda-components-08`;
                       operator-run gated)
                │
                ▼
   Phase D-CUDA arm  →  2×2 + CPU↔CUDA cross-check (operator-run gated)
                │
                ▼
   Phase E ✅ E.8 per-stage timing; E.2/E.3/E.4 cluster-gated
                │
                ▼
   Phase F ✅ tooling (lmp_cuda, PDC gpugh sbatch, multi-engine validator)
```

All code-side work in the dependency graph is landed; what remains is the **operator-action** D-CPU re-verify pass (post benchmark-fidelity fixes) and the Phase I CUDA validation runs on PDC Dardel. The CPU arm of Phase D is delegated to the `msu-pe-phase-d-*` sub-chain and `msu-pe-align-03`; the CUDA arm is delegated to chain link `msu-pe-cuda-components-08`.

---

## Phase G — Functional-gap implementation ✅ landed

All G.1–G.9 sub-items landed via the `msu-pe-cuda-components-01..06` chain. Summary of what shipped:

- **G.1 `ZeroMomentum` (CPU + CUDA).** CPU `cpu::md::ZeroLinearMomentum` (with `using ZeroMomentum = ZeroLinearMomentum;` alias) in `include/atomiverse/cpu/md/init.h` — bound at `end_of_step @ +25`, multi-rank-correct via 4-double `MPI_Allreduce` of `{Σm, Σm·v_d}`. CUDA `cuda::md::ZeroMomentum` in `include/atomiverse/cuda/md/init.cuh`. Benchmark stages 2/3/4 wired (`benchmarks/md/polyethylene_tension.cpp`).
- **G.2 Per-axis `Pxx/Pyy/Pzz`.** `polyethylene_tension.cpp` thermo emits real diagonal-virial columns; `tools/run_msu_atomiverse.sh` no longer fakes (per the inline `Phase G.2 cleared` note in the script).
- **G.3 Energy decomposition.** thermo CSV emits `E_pair/E_bond/E_angle/E_dihedral`; `frame.scalars.energy_{pair,bond,angle,dihedral}` reduced across ranks in the same 18-double `MPI_Allreduce` as the rest of the scalars.
- **G.4 CUDA categorical API.** `cuda::md::MD` mirrors CPU `set_pair/set_bond/set_angle/set_dihedral/set_neighbor/set_integrator/set_thermostat/set_barostat/set_observer/set_record/set_thermo`.
- **G.5–G.7 CUDA integrators.** `cuda::md::NVELimit`, `cuda::md::FixLangevin`, `cuda::md::NoseHooverNPT` (iso + aniso y/z).
- **G.8 CUDA LJ special-bonds.** `PairLJCut::set_special_bonds` + kernel exclusion masking.
- **G.9 CUDA benchmark path.** `--backend {cpu|cuda}` in `polyethylene_tension.cpp`; CUDA arm uses G.4 categorical API.

Source of truth: chain spec files `msu-pe-cuda-components-{01..06}.md`; commits `02588b4`, `4fa13ec`, `85e6881`, `3d399e7`, `40c3fab`, `4c3e109`. Exit gate (a)–(d) of original Phase G satisfied.

---

## Phase H — Async correctness ✅ H.1–H.3 landed, H.4 9/10 ACs

H.1 (CPU `Thermo` → `AsyncComponentBase`), H.2 (CUDA driver hot-path `defer_stream` desync via `cuda::AsyncRecordBase` + zero-`cudaStreamSynchronize` hot-path guard), H.3 (CUDA `VirialStress` deferrable when observed-only) all landed via chain `msu-pe-cuda-components-07` (commit `3a4ad64`). H.4 exit-gate micro-benchmark passes the 25×-frequency wall-clock < 5% target; one acceptance criterion (ac-007 `compute-sanitizer memcheck`) deferred — the tool is absent from `~/atv-build-env`.

---

## Phase I — CUDA single-GPU 5-stage correctness validation

Goal: mirror the CPU Phase A–C structure for the CUDA backend on a single GPU — prove the CUDA force field is correct, the full 5-stage protocol completes, and CUDA agrees with both the Atomiverse CPU run and LAMMPS GPU mode. Blocked by G (components) and H (async). All cross-engine comparisons are statistical block-averages per `msu-pe-align-03`.

### I.1 CUDA force-field correctness

New `tests/cuda/md/test_pe_cuda_ff.cpp` (new): evaluate `E_pair/E_bond/E_angle/E_dihedral` for the C₂₀H₄₂ perturbed-geometry cross-check (the same geometry as `tests/cpu/md/test_pe_energy_xcheck.cpp`, with the closed-form perturbation in `SOURCES.md`) on CUDA. Because CUDA is single-precision the tolerance is statistical/loose, not the CPU 1e-6 eV bound: each energy term agrees with the CPU double-precision value to a relative tolerance consistent with float accumulation over 10 000 atoms (target ≤ 1e-3 relative). Confirms G.5–G.8 components produce physically correct forces, including special-bonds exclusions.

### I.2 CUDA full 5-stage protocol completes

CUDA path of `polyethylene_tension.cpp` runs the full 5-stage protocol on `msu_pe.data` on a single GH200 GPU to completion without NaN, writing the canonical CSV schema (thermo with per-axis P + energy decomposition, stress-strain). Operator action on the PDC `gpugh` partition. Density after stage 4 within 5% of the MSU reported amorphous-PE density at 100 K (~0.94 g/cm³).

### I.3 CUDA ↔ CPU cross-check

`tools/compare_msu.py` compares the Atomiverse-CUDA and Atomiverse-CPU runs on block-averaged ensemble means: ⟨T⟩, ⟨E_pot⟩ and the energy decomposition over stage 4, density, and the three stress-strain descriptors. Tolerance is the statistical bar widened for single vs double precision (target ~2–3% relative on means, ±5% on stress-strain descriptors). This isolates precision/RNG effects from genuine engine bugs before the LAMMPS comparison.

### I.4 CUDA ↔ LAMMPS-GPU statistical comparison

`tools/compare_msu.py` + `tools/plot_msu_compare.py --include-cuda` compare the Atomiverse-CUDA run against the `lmp_cuda` GPU-mode reference run on block-averaged ensemble means — the GPU-mode arm of the 2×2. Same statistical methodology (Flyvbjerg & Petersen 1989) and tolerances as the CPU↔LAMMPS arm.

**Exit gate for Phase I**: I.1 green; I.2 completes the 5-stage CUDA run with density within 5% of reference; I.3 CUDA↔CPU statistical comparison passes; I.4 CUDA↔LAMMPS-GPU statistical comparison passes.

---

## Phase D — Dual-backend LAMMPS comparison (2×2 + CPU↔CUDA cross-check)

Goal: demonstrate Atomiverse and LAMMPS, starting from the same `msu_pe.data` and running the same protocol, produce statistically matching observables on **both** backends. Phase D is now dual-backend and four-config.

### D-CPU arm — tooling ✅ landed, operator-action gated

The CPU arm code, scripts, and analysis tools are all landed (the original sub-chain specs `msu-pe-phase-d-01/03/04` + `msu-pe-align-03-statistical-thermo` shipped and were deleted on 2026-05-24 — git history retains them). Concrete artifacts:

- `tools/run_msu_lammps.sh --mode cpu` + `tools/run_msu_atomiverse.sh` (default CPU) produce reference CSVs.
- `tools/compare_msu.py` (was `validate_polyethylene.py` before renamings) — statistical block-average thermo parity per Flyvbjerg & Petersen 1989 block averaging.
- Equilibrium observables (ρ, ⟨E_pot⟩, RDF peak, ⟨R²⟩) and stress-strain descriptors (E, σ_y, hardening slope ±5%) are part of the validator's invariant set.

Remaining: the operator-action LAMMPS-CPU + Atomiverse-CPU runs on PDC Dardel `main` to re-verify after the recent benchmark-fidelity fixes (G.1–G.3, E.8 timing). The CUDA arm below stays consistent: same statistical methodology, same `compare_msu.py` surface, same exit-gate shape.

### D.5 The 2×2 comparison matrix

Four engine configs, three comparison tracks — all statistical block-average of ensemble means:

| Track | Engine A | Engine B | Methodology |
|---|---|---|---|
| CPU reference | Atomiverse-CPU | `lmp_cuda` CPU-mode | block-average (Flyvbjerg & Petersen 1989) |
| GPU reference | Atomiverse-CUDA | `lmp_cuda` GPU-mode | block-average — Phase I.4 |
| CPU↔CUDA cross-check | Atomiverse-CPU | Atomiverse-CUDA | block-average — Phase I.3 |

Compared observables in every track: energy (incl. the `E_pair/E_bond/E_angle/E_dihedral` decomposition from G.3), density ρ, equilibrium observables (RDF peak, ⟨R²⟩), and the stress-strain curve descriptors (E, σ_y, hardening slope). `tools/compare_msu.py` accepts the four engine CSVs and runs all three tracks; `tools/plot_msu_compare.py --include-cuda` overlays them.

**Exit gate for Phase D**: all three tracks pass per the block-average methodology; `tools/compare_msu.py` exits 0.

---

## Phase E — Performance benchmarks (CPU scaling unchanged + per-stage timing)

The CPU/MPI scaling work (E.1–E.7: single-rank breakdown, strong scaling, weak scaling, head-to-head vs LAMMPS-CPU, profile decomposition, hybrid MPI×OpenMP, multi-node) stays **CPU-only and unchanged** — see prior spec body. Multi-GPU / CUDA+MPI scaling is explicitly out of scope. One new sub-item is added:

### E.8 Per-stage wall-clock timing across four engine configs

✅ landed (commit `4c3e109`). Per-stage wall-clock timer in `benchmarks/md/polyethylene_tension.cpp` brackets each of the five stages and emits a per-stage timing CSV. Cross-config comparison (Atomiverse-CPU, Atomiverse-CUDA, `lmp_cuda` CPU-mode, `lmp_cuda` GPU-mode) is descriptive — not a hard pass/fail gate, except that Atomiverse-CUDA per-stage time must be finite and the run must complete.

**Exit gate for Phase E**: E.1–E.5 pass their CPU efficiency targets (unchanged, cluster-gated); E.8 emits a four-config per-stage timing table. E.6/E.7 remain stretch.

---

## Phase F — Reproducibility artifacts (dual-backend, PDC Dardel)

Everything an external reviewer needs to reproduce the dual-backend results on PDC Dardel.

### F.1 Input files under version control

Unchanged set under `benchmarks/md/reference/polyethylene/` (`force_field.yaml`, `msu_pe.data`, `lammps_msu_paper.in`, the reference CSVs, `SOURCES.md`, `README.md`), plus the GPU-mode and CUDA reference CSVs (engine-suffixed per the `SOURCES.md` Phase D.1 filename scheme, extended with a GPU/CUDA suffix).

### F.2 Scripts — `lmp_cuda` + PDC `gpugh`

- `tools/run_msu_lammps.sh` — extended to launch `lmp_cuda` in **both** CPU mode and GPU mode (not the old CPU/MPI-only path); PDC `sbatch` submission, `-A <project>` account flag, `main` partition for CPU mode and `gpugh` partition for GPU mode.
- `tools/run_msu_atomiverse.sh` — extended with a CUDA / GPU-mode launch path on the `gpugh` partition; stops faking per-axis P (G.2).
- `tools/run_msu_pdc.sh` (new) — PDC Dardel SLURM submission wrapper: `sbatch`/`srun`, `-A` account flag (concrete values: `naiss2025-1-25` for the `main` CPU partition, `naiss2025-1-25-gpu` for the `gpugh` GPU partition), partition selection (`main` vs `gpugh`), 24h wall-time cap awareness. The flag stays required with no default in the script CLI; the concrete account belongs in docs/reproduction commands, not hardcoded.
- `tools/compare_msu.py` + `tools/plot_msu_compare.py` — the multi-engine 2×2 validator (D.5 / E.8): three comparison tracks, energy-decomposition columns, per-stage timing, Flyvbjerg-Petersen block-average methodology.
- `tools/plot_scaling.py`, `tools/plot_profile.py` — unchanged.

### F.3 README

`benchmarks/md/reference/polyethylene/README.md` updated: the `lmp_cuda` build/usage (CPU mode + GPU mode), the PDC Dardel `gpugh` / `main` partition `sbatch` procedure with the `-A` account flag, the four-config 2×2 reproduction commands, and the dual-backend expected-results table. Reference: <https://support.pdc.kth.se/doc/run_jobs/job_scheduling/>.

**Exit gate for Phase F**: a fresh clone plus the documented operator commands reproduce all four engine runs on PDC Dardel and `tools/compare_msu.py` runs the three 2×2 tracks + the timing track green.

---

## Tasks

All implementation tasks below are ✅ landed; remaining work is operator-action SLURM runs on PDC Dardel.

- [x] Implement ZeroMomentum (CPU `cpu::md::ZeroLinearMomentum` + CUDA `cuda::md::ZeroMomentum`), CUDA NVELimit / FixLangevin / NoseHooverNPT, CUDA LJ special-bonds exclusions (Phase G.1/G.5/G.6/G.7/G.8)
- [x] Implement CUDA categorical setter API in `include/atomiverse/cuda/md/md.cuh` + `src/cuda/driver.cu` (Phase G.4)
- [x] Implement per-axis Pxx/Pyy/Pzz + E_pair/E_bond/E_angle/E_dihedral CSV columns and CUDA path in `benchmarks/md/polyethylene_tension.cpp` (Phase G.2/G.3/G.9)
- [x] Implement CPU Thermo on AsyncComponentBase + CUDA driver hot-path defer_stream desync + VirialStress deferrable (Phase H.1/H.2/H.3)
- [x] Async-correctness gate: 25×-frequency wall-clock < 5% (Phase H.4) — 9/10 ACs passed
- [x] CUDA force-field correctness against CPU double-precision values (Phase I.1)
- [x] Multi-engine comparison tooling — `tools/compare_msu.py` + `tools/plot_msu_compare.py --include-cuda` (Phase D.5/E.8)
- [x] `lmp_cuda` CPU/GPU-mode launch + PDC `gpugh` sbatch — `tools/run_msu_lammps.sh`, `tools/run_msu_atomiverse.sh`, `tools/run_msu_pdc.sh` (Phase F.2)
- [x] Per-stage wall-clock timer (Phase E.8)
- [ ] **Operator-action.** Submit Phase D-CPU re-verify reference runs (LAMMPS-CPU + Atomiverse-CPU on `main` partition) on PDC Dardel after the recent benchmark-fidelity fixes
- [ ] **Operator-action.** Submit Phase I CUDA validation runs (LAMMPS-GPU + Atomiverse-CUDA on `gpugh` partition) per `msu-pe-cuda-components-08`
- [ ] **Operator-action.** Run `tools/compare_msu.py` against all four engine CSVs and confirm the four invariants (T_stage4 within 15 K, ρ_stage4 within 0.10 g/cm³, E within 15%, σ_flow within 0.2 GPa)

## Testing strategy

- **Happy path** — `tests/cuda/md/test_cuda_md_components.cpp` and `tests/cpu/md/test_zero_momentum.cpp` exercise each new component on a small synthetic system: `ZeroMomentum` drives COM momentum to zero; `NVELimit` caps displacement at the configured bound; `FixLangevin` drives temperature to the setpoint with a fixed seed; `NoseHooverNPT` relaxes a box to target pressure in both iso and aniso y/z modes; CUDA `LJ` with special-bonds excludes 1-2/1-3/1-4 pairs. The CUDA benchmark path runs a short (≤1000-step) protocol without NaN.
- **Edge cases** — `NVELimit` on the overlapped `msu_pe.data` (initial E_pair ≈ 3.4×10⁸ kcal/mol) does not blow up; CUDA `LJ` with an empty special-bonds list behaves identically to the un-masked kernel; `tools/compare_msu.py` with a missing engine CSV exits non-zero with a clear error naming the absent track; the timing track on a CSV missing a stage column fails gracefully, not with a crash.
- **Domain validation** — CUDA force field: `test_pe_cuda_ff.cpp` confirms each of `E_pair/E_bond/E_angle/E_dihedral` agrees with the CPU double-precision reference to ≤1e-3 relative on the C₂₀H₄₂ perturbed geometry. Stochastic-protocol parity (Flyvbjerg & Petersen 1989 block averaging): the CUDA↔CPU and CUDA↔LAMMPS-GPU comparisons use block-averaged ensemble means (⟨T⟩, ⟨E_pot⟩, energy decomposition, density) over the equilibrated stages, never per-step; CUDA single precision widens the means tolerance to ~2–3% relative. Stress-strain descriptors (E, σ_y, hardening slope) agree to ±5%. Reference fidelity: an assertion confirms the benchmark applies a `fix momentum` equivalent in stages 2/3/4 and emits real (non-faked) `Pxx/Pyy/Pzz`.
- **Async correctness** — the H.4 micro-benchmark measures steps/s with recorders/observers off vs at high frequency; the gate requires high-frequency steps/s ≥ 95% of recorders-off steps/s on CUDA, proving the GPU never stalls on a recorder; the CPU counterpart confirms the same for the now-async `Thermo`.

## Out of scope

- Multi-GPU and CUDA+MPI scaling — explicitly out of scope per fixed scope decision 1; CUDA validation is single-GPU correctness only. The CPU/MPI Phase E (E.1–E.7 strong/weak scaling) stays CPU-only and unchanged. (Considered adding a CUDA scaling phase; rejected — the user fixed CUDA scope to single-GPU correctness.)
- Re-specifying the CPU Phase D sub-chain — the original `msu-pe-phase-d-01..04` and `msu-pe-align-03` specs shipped and were deleted on 2026-05-24 once their code-side ACs landed (git history retains them). The umbrella tracks the remaining operator-action work directly.
- Per-step or bit-level comparison against LAMMPS — physically impossible for this stochastic, single-precision-on-CUDA protocol; all cross-engine comparison is statistical block-average (Flyvbjerg & Petersen 1989).
- Matching MSU paper numbers directly — LAMMPS under identical inputs is the reproducible reference; MSU paper-number matching remains a downstream addendum, unchanged from the prior spec rationale.
- Changing the protocol physics — `lammps_msu_paper.in` is immutable ground truth; the benchmark is corrected to match it, never the reverse.

## Rationale

- **Why supersede in place instead of splitting?** The user explicitly asked to update this one umbrella file. The umbrella already delegates detailed implementable slicing to per-phase chains (`msu-pe-phase-d-*`); the new CUDA-component work likewise delegates to a `msu-pe-cuda-components-*` chain at impl time. Keeping one umbrella preserves the single roadmap view.
- **Why a dedicated functional-gap phase (G) before any CUDA validation?** The user said missing functionality must be implemented first. The CUDA backend literally cannot run the protocol — it has no `NVELimit`, `FixLangevin`, `NoseHooverNPT`, no special-bonds exclusions (which silently corrupt the polymer), and no categorical setter API. Phase G hard-blocks the CUDA arm; making it explicit prevents starting validation on a non-functional backend.
- **Why call out `fix momentum`, per-axis P, and energy decomposition as fidelity defects?** The benchmark does not yet strictly match `lammps_msu_paper.in`: it omits `fix momentum` (a physics deviation), fakes `Pxx/Pyy/Pzz` by copying isotropic P (wrong for the stress-strain curve), and omits the energy decomposition the user wants to compare. These are correctness defects, not enhancements.
- **Why statistical block-average for every CUDA comparison?** The MSU protocol is stochastic (seeded stage-1 Langevin) and CUDA additionally runs single precision. Two engines therefore never share a trajectory — they sample the same ensemble. The only physically meaningful cross-engine statement is agreement of block-averaged ensemble means (Flyvbjerg & Petersen 1989).
- **Why a dedicated async-correctness phase (H)?** The user's hard requirement is that the GPU never waits on a recorder. The CUDA `AsyncBridge` exists but is defeated by per-phase-slot `cudaStreamSynchronize(defer_stream)` in the driver hot path. An explicit phase with a measurable exit gate (step time independent of recorder frequency) turns "make it async" into a verifiable contract.
- **Why PDC Dardel as the cluster?** Fixed scope decision 4: the reference runs are ~24h-class operator actions; PDC Dardel's `main` (CPU) and `gpugh` (GH200 GPU) partitions map cleanly onto the 2×2 engine matrix, and the 24h job cap is why the runs are operator actions, not CI steps.

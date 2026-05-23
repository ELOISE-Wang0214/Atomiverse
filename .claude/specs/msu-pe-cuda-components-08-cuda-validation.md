---
title: MSU PE CUDA components — link 08 — CUDA vs LAMMPS-GPU statistical validation
status: code-complete
created: 2026-05-23
parent: msu-pe-validation.md (Phase I)
chain: msu-pe-cuda-components-* (link 8/8)
---

# MSU PE CUDA components — link 08 — CUDA vs LAMMPS-GPU statistical validation

## Summary

链节 01–07 落地后,把 CUDA 后端跑出来的 MSU PE 五阶段产物(thermo、stress-strain)与 LAMMPS `lmp_cuda` GPU 模式同协议产物做 **统计 block-average 等价性测试**:扩展 `canonicalize_msu_output.py` 接受 `engine=atomiverse` + `backend=cuda` 复合标签,扩展 `tools/compare_msu.py`(若不存在则新建)做跨后端 + 跨引擎比较,生成验证报告 `runs/msu_pe_cuda_validation_<ts>.md`,出图 `tools/plot_msu_compare.py` 加 CUDA-engine overlay。退出门:stage-4 ⟨T⟩、⟨ρ⟩ 跨引擎 block-average 在 3σ 内一致;stage-5 stress-strain 曲线在弹性段 (`ε < 0.02`) 模量 (`E = dσ/dε`) 跨引擎相对差 < 15%;塑性段 (`ε > 0.1`) 流动应力差 < 0.2 GPa。该链节是父规范 `msu-pe-validation` Phase I 的终极退出门 —— 通过即父规范从 `in-progress` 落到 `done`。

## Domain basis

统计等价性方法学 — block-average + Flyvbjerg-Petersen error analysis(*J. Chem. Phys.* 91 (1989) 461, doi:10.1063/1.457480):MD 序列存在自相关,直接算 σ/√N 会低估方差;先把序列切成长度 ≥ 5τ_auto 的 block,每个 block 内取平均得到一组近似独立样本,再算它们的样本标准差。两组样本均值在 ±3σ_combined 内一致即视为统计等价(σ_combined = √(σ_A²/N_A + σ_B²/N_B))。

MSU PE 协议的物理 invariants(链节 04 / 06 / 父规范 G.* 阶段已落):

- **Stage 2/3/4 ⟨T⟩** — NPT 控温收敛到设定值;stage 4 ⟨T⟩ = 100 K (设定),tolerance 来自 GH200 单精度 + 小盒子(8 链 × 30 CHₓ ≈ 240 原子)统计涨落,经验上 ±5 K 一倍 σ。
- **Stage 4 ⟨ρ⟩** — 非晶 PE 在 100 K, 0 atm 收敛到 ρ ≈ 0.92 ± 0.03 g/cm³ (Tschopp 2010,*Polymer*, Figure 4 报告 ρ = 0.91 g/cm³ at 1 atm, 300 K;100 K 略高;3σ 范围即 ±0.09 g/cm³)。
- **Stage 5 弹性模量 E** — 单轴 erate = 1×10⁻⁵ /fs,Dreiding-UA 在 100 K 实测 E ≈ 2.0 ± 0.4 GPa (Tschopp 2010 Figure 7,seed-dependent;链节早期 CPU 跑 12345 seed 给 2.18,7777 给 1.95)。CUDA 单精度 + 同 seed 经验上漂移 < 5%。
- **Stage 5 流动应力 σ_flow** — 大应变区 (ε > 0.1) σ_xx 平台值 ≈ 0.5 ± 0.1 GPa (Tschopp 2010 Figure 7)。

跨引擎门(基于 CPU vs LAMMPS-CPU 已知容差 + 单精度展宽):

- **⟨T⟩_stage4** : |⟨T⟩_cuda − ⟨T⟩_lammps_gpu| < 3·√(σ_T_cuda² + σ_T_lammps²) , 实际 ~ 15 K(σ_T ≈ 5 K)。
- **⟨ρ⟩_stage4** : |⟨ρ⟩_cuda − ⟨ρ⟩_lammps_gpu| < 0.10 g/cm³(3·σ_ρ)。
- **E_modulus_stage5** : |E_cuda − E_lammps_gpu| / E_lammps_gpu < 0.15 (15% 单精度 + seed 漂移合计)。
- **σ_flow_stage5** : |σ_flow_cuda − σ_flow_lammps_gpu| < 0.2 GPa。

单位:T [K],ρ [g/cm³],E [GPa],σ_xx [GPa]。

## Design

**总体策略:不重写已有 canonical 流水线,只扩展。** `canonicalize_msu_output.py` 已经处理 `engine ∈ {atomiverse, lammps}`,现增加 `--backend {cpu, cuda}` 维度;输出文件名带后缀(`thermo_atomiverse_cuda_1r.csv`)。`tools/compare_msu.py` 是新工具(目前 `plot_msu_compare.py` 是 plot-only,不出验证报告),做四路 (engine ∈ {atomiverse, lammps}) × (backend ∈ {cpu, cuda}) 比较;但本链节的退出门只关心 **`(atomiverse, cuda)` vs `(lammps, cuda)`**,其它三路(cpu/cpu, atomiverse/cpu vs lammps/cpu, atomiverse/cuda vs atomiverse/cpu)作为旁路在报告里列出。

**`canonicalize_msu_output.py` 扩展:**

- 新增 CLI flag `--backend {cpu, cuda, auto}` (default `auto`):`auto` 时从输入文件路径模式推断(`runs/atomiverse_cuda_*` → `cuda`,`runs/lammps_*_gpu_*` → `cuda`,其余 → `cpu`)。
- 写出 canonical CSV 时 header 加 metadata 行 `# engine=<E>,backend=<B>,seed=<S>` (注释起首,canonical reader 忽略)。
- 文件名:`thermo_<engine>_<backend>_<runs>r.csv` / `stress_strain_<engine>_<backend>_<runs>r.csv`。注意现有 `thermo_atomiverse_1r.csv` / `thermo_lammps_1r.csv` 无 backend 标签 → 视为 `cpu` 默认,保持向后兼容(`compare_msu.py` 读取时 backend 缺省映射为 cpu)。

**`tools/compare_msu.py` (new) — 验证报告生成器:**

```
usage: compare_msu.py <thermo_A> <thermo_B> <stress_A> <stress_B> \
                     [--label-A "atomiverse-cuda"] [--label-B "lammps-cuda"] \
                     --out-md <report.md> [--block-size 250] [--equilibration 25000]
```

- 读两套 (thermo, stress_strain) canonical CSV,按 stage 切片。
- 对每 stage 跑 Flyvbjerg-Petersen block-average,自动选 block size(从 256 步起倍增,直到 σ_block 平台);若用户传 `--block-size` 则用固定值。
- 计算四项 invariant:
  1. ⟨T⟩_stage4 ± σ
  2. ⟨ρ⟩_stage4 ± σ (ρ 从 box volume + 已知 atom mass total 反推:ρ = M_total / V_box)
  3. E_modulus_stage5(对 stage 5 前 5% 数据 — ε ∈ [0, 0.05] — 做最小二乘 σ_xx = E·ε,σ_E 用 jackknife)
  4. σ_flow_stage5(stage 5 后 50% 数据 — ε ∈ [0.5, ε_max] — 取均值 ± σ)
- 对每项算 `|A − B| − 3·σ_combined`;< 0 即 PASS,>= 0 即 FAIL。
- 写 markdown 报告:每项一行 (`| metric | A | B | Δ | 3σ | verdict |`),底部 PASS/FAIL overall;若 FAIL 给出最严重的违例。

**`tools/plot_msu_compare.py` 扩展 — 多 engine overlay:**

- 现版本支持 CPU AV vs LAMMPS;新增 `--include-cuda` flag 多画一组 (atomiverse-cuda, lammps-cuda) 曲线,不同颜色和线型,legend 区分 4 路。
- 输出 SVG 路径加后缀 `_with_cuda.svg`,保留旧无 cuda 版本不动。

**LAMMPS GPU 参考运行 — 跑通 lmp_cuda 在 GH200:**

- `tools/run_msu_lammps.sh` 已有 `--mode cpu` 默认;新增 `--mode gpu` 在脚本里 `mpirun lmp_cuda -pk gpu 1 -sf gpu -in lammps_msu_paper.in`(LAMMPS GPU package 的 PerKappaGPU 风格)。
- 单 GH200 跑 `--scale 1.0`(完整 171 800 步 stage 5)wall-clock 经验 ~ 40–60 min(链节早期 CPU lmp_cuda 单 rank 跑 41 min)。
- 同 seed (`904297`) 与 Atomiverse CUDA path 对齐,初速度种子一致以排除引擎间初态差。

**Atomiverse CUDA 参考运行:**

- 由链节 06 的 `tools/run_msu_atomiverse.sh --backend cuda` 落实;本链节只调用,不修改。
- 同 `--scale 1.0`,同 seed,产物 `runs/atomiverse_cuda_1r_<ts>/`。

**报告 + 提交流水线:**

1. 跑 LAMMPS GPU 参考运行 → `runs/lammps_cuda_1r_<ts>/`(操作员 action,SLURM)。
2. 跑 Atomiverse CUDA 运行 → `runs/atomiverse_cuda_1r_<ts>/`(操作员 action,SLURM)。
3. canonical:`canonicalize_msu_output.py runs/<ts> --engine atomiverse --backend cuda`;`canonicalize_msu_output.py runs/<ts> --engine lammps --backend cuda`。
4. compare:`compare_msu.py thermo_atomiverse_cuda_1r.csv thermo_lammps_cuda_1r.csv stress_strain_atomiverse_cuda_1r.csv stress_strain_lammps_cuda_1r.csv --out-md runs/msu_pe_cuda_validation_<ts>.md`。
5. plot:`plot_msu_compare.py --include-cuda --out msu_pe_validation_cuda_<ts>.svg`。
6. 报告写入 `benchmarks/md/reference/polyethylene/reports/`(版本控制),svg 同目录。

**显式不做的事:**

- 不动 cpu 路径下的 `compare_msu.py`(本链节首次创建,但只覆盖 cuda 子用例;cpu vs lammps-cpu 已在父规范 G 阶段早期验证完毕,作为旁路打印不当退出门)。
- 不实现 multi-seed ensemble(每 engine 跑一次,seed 对齐即可;父规范 Out of scope #2)。
- 不做温度扫描或应变率扫描(单点验证,单 seed,单 scale)。
- 不修复链节 01–07 暴露的任何回归 —— 那些是各自链节的职责;本链节只在前 7 个全部通过后才能跑。

## Files to create or modify

- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tools/canonicalize_msu_output.py` — 加 `--backend {cpu,cuda,auto}` 参数 + 文件名后缀逻辑 + metadata 注释行;auto-detect 从路径推断 backend。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tools/test_canonicalize_msu_output.py` — 加 backend auto-detect 测试 + 文件名 suffix 测试(2 个新测试,总数 21)。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tools/compare_msu.py` (new) — Flyvbjerg-Petersen block-average + 4-invariant verdict + markdown 报告。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tools/test_compare_msu.py` (new) — pytest:合成两份已知统计差异的数据集,断言报告 verdict 正确(PASS / FAIL 两个用例 + block-size auto 收敛测试)。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tools/plot_msu_compare.py` — 加 `--include-cuda` flag,多画一组 (atomiverse-cuda, lammps-cuda) 曲线。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/tools/run_msu_lammps.sh` — 加 `--mode gpu` 分支(`-pk gpu 1 -sf gpu`)。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/benchmarks/md/reference/polyethylene/reports/msu_pe_cuda_validation_<ts>.md` (new at run time) — 由 `compare_msu.py` 生成,不进 git index(但生成路径固定)。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/.claude/specs/INDEX.md` — 加 `msu-pe-cuda-components-08-cuda-validation`。
- `/cfs/klemming/home/j/jichenli/work/Atomiverse/.claude/specs/msu-pe-validation.md` — 父规范 Phase I checklist 收尾标记(此链节通过后,父规范从 `in-progress` 落到 `done` —— 但只在所有 ac 通过且操作员确认 LAMMPS GPU run 完成后落)。

## Tasks

- [x] Write failing test `tools/test_compare_msu.py` (new) covering: (a) synthetic PASS case (matched datasets), (b) synthetic FAIL case (datasets with 6σ shift in T_stage4), (c) auto block-size convergence (`σ_block` plateau detection)
- [x] Add 2 failing tests to `tools/test_canonicalize_msu_output.py`: `--backend cuda` auto-detect from `runs/atomiverse_cuda_*` path; output filename has `_cuda` suffix
- [x] Implement `--backend {cpu,cuda,auto}` in `tools/canonicalize_msu_output.py` + auto-detect + filename suffix + metadata comment
- [x] Implement `tools/compare_msu.py` with Flyvbjerg-Petersen block-average, 4-invariant verdict, markdown report
- [x] Implement `--include-cuda` overlay in `tools/plot_msu_compare.py`
- [x] Add `--mode gpu` branch to `tools/run_msu_lammps.sh` (`-pk gpu 1 -sf gpu`)
- [ ] Operator action: submit LAMMPS GPU reference run on Dardel gpugh / single GH200 with `--scale 1.0 --seed 904297` (ac-007 — gpu_required)
- [ ] Operator action: submit Atomiverse CUDA run with `tools/run_msu_atomiverse.sh 1 --backend cuda -- --mode full --scale 1.0 --seed 904297` (ac-007 — gpu_required)
- [ ] Canonicalize both run outputs into engine+backend-tagged CSVs (ac-008 — depends on ac-007)
- [ ] Run `tools/compare_msu.py` and capture markdown report (ac-009 — depends on ac-008)
- [ ] Run `tools/plot_msu_compare.py --include-cuda` and capture overlay SVG (ac-010 — depends on ac-008)
- [ ] Verify all 4 invariants PASS (T_stage4 within 15 K, ρ_stage4 within 0.10 g/cm³, E within 15%, σ_flow within 0.2 GPa) (ac-009)
- [x] Add `msu-pe-cuda-components-08-cuda-validation` entry to `.claude/specs/INDEX.md`
- [ ] On full PASS: flip parent `msu-pe-validation.md` Phase I status checkmark and surface to user for final approval to mark parent `done` (parent-gate, awaits ac-009 outcome)

## Testing strategy

- **Happy path — synthetic PASS:** `test_compare_msu.py::test_synthetic_pass` 构造两份 T_stage4 完全相同的合成 CSV,断言报告 verdict 全 PASS,markdown 输出有 4 行 `verdict: ✓`。
- **Happy path — synthetic FAIL:** `test_compare_msu.py::test_synthetic_fail` 构造 T_stage4 差 30 K(6σ)的两份 CSV,断言报告标 `T_stage4: ✗` 且 overall FAIL。
- **Edge — block size auto:** 给合成序列(已知 τ_auto = 200 步),断言 auto block-size 收敛到 ≥ 1000 步(5τ_auto)。
- **Edge — backend auto-detect:** `test_canonicalize_msu_output.py::test_backend_autodetect_cuda` 用 path `runs/atomiverse_cuda_1r_20260523/thermo.csv` 测 auto → cuda。
- **Domain validation:** *本链节自身是 domain validation* —— ac 由实际 LAMMPS GPU vs Atomiverse CUDA 运行的统计比较决定,无独立 domain unit test。
- **Integration:** 跑完两个完整运行后,`compare_msu.py` 报告必须 PASS;若 FAIL,debug 流程是 ① 检查 LAMMPS deck 与 Atomiverse benchmark 协议是否字对字一致(seed、dt、步数),② 拆分到单 stage 比较定位漂移源,③ 若是单精度问题,放宽对应 invariant 的 3σ 门(写入父规范 NOTES,不在本链节落地)。

## Out of scope

- Multi-seed ensemble averaging — single seed 904297 即可(父规范 fixed-scope #2)。
- 温度扫描 / 应变率扫描 — 单点验证 (100 K, erate=1e-5)。
- LAMMPS GPU vs LAMMPS CPU 内部一致性 — 假定 LAMMPS 自己已经通过 KOKKOS/GPU 测试套件;本链节只把它当黑盒参考。
- Atomiverse 多 GPU / MPI — 父规范 fixed-scope #1。
- 提升单精度到双精度 — 性能代价过高,3σ 门是为单精度设计的。若验证 FAIL 且确认是单精度,讨论双精度方案在独立链节(`atv-cuda-double-precision`)。
- 参考报告的 git commit — markdown 报告与 SVG 是 run-time artifact,放 `benchmarks/md/reference/polyethylene/reports/` 但是否进 git index 由用户决定(本链节不强制)。
- 性能比较 (Atomiverse CUDA vs LAMMPS GPU wall-clock) — 父规范 Phase F 的事;本链节只看正确性。

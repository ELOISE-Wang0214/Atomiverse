---
title: MSU-PE Phase D.1 — LAMMPS + Atomiverse-CPU reference runs
status: code-complete
created: 2026-05-22
parent: msu-pe-validation.md
---

# MSU-PE Phase D.1 — LAMMPS + Atomiverse-CPU reference runs

## 综述 (Summary)

本规范是 `msu-pe-phase-d-*` 链的第 1 环(共 4 环,父规范 `msu-pe-validation.md` 的 Phase D,**CPU 臂**)。它产出后续三个已落地子环——`msu-pe-align-03`(统计热力学)、`msu-pe-phase-d-03`(平衡态观测量)、`msu-pe-phase-d-04`(应力-应变)——共同消费的**参考运行产物**:LAMMPS-CPU 与 Atomiverse-CPU 两个引擎、各跑 1 与 8 个 MPI rank,产出 8 个规范模式参考 CSV(`thermo_<engine>_<N>r.csv` 与 `stress_strain_<engine>_<N>r.csv`,N ∈ {1,8}),外加阶段-4 热力学日志、最终 RDF 曲线与轨迹 XYZ。本规范不写任何生产 C++ 代码:它只在 `tools/` 与 `benchmarks/md/reference/polyethylene/` 内工作——审计 `lammps.in` 的物理保真度、把 `run_msu_lammps.sh` 切到 `lmp_cuda` CPU 模式、给两个运行脚本加 PDC Dardel `main` 分区的 `sbatch` 提交、并实现/验证“原始引擎输出 → 规范模式 CSV”的归一化(canonicalization)步骤。用户可见的成果:一套可审计、可复现、模式正确的参考产物,使下游三个子环从“被链 01 阻塞”转为“可机验”。

本规范诚实地区分两类交付:**(a) 现在即可实现并机验的代码**——deck 保真度审计、脚本切换与 `sbatch` 提交、归一化逻辑、以及一个短协议冒烟测试;**(b) 操作员运行时门**——四次 PDC Dardel 全协议参考运行(每次 ~24 小时级墙钟,是操作员动作而非 CI 步骤),以及落地 8 个规范 CSV + 阶段-4 热力学/RDF/轨迹产物。后者在代码完成后停泊(parked),与姊妹规范 `msu-pe-phase-d-04` 的 ac-002 同形。

## Chain dependencies

本规范是 `msu-pe-phase-d-*` 链的第 1 环(共 4 环)。链内编写/落地顺序严格为 **01 → 02 → 03 → 04**:

- **本环(01)无上游链依赖** —— `conflict_decision: independent`。它是链的起点,只读已存在的 `lammps.in` / `lammps_msu_paper.in` / `msu_pe.data` 与已存在的运行脚本。
- **下游消费者(不在此重新规范,仅需保持一致)**:
  - `msu-pe-align-03-statistical-thermo`(取代 D.2)读 `thermo_lammps_<N>r.csv` 与 `thermo_atomiverse_<N>r.csv`,N ∈ {1,8}。
  - `msu-pe-phase-d-03-equilibrium` 读阶段-4 热力学日志(volume + E_pot 每步)、最终 RDF 曲线、轨迹 XYZ,两引擎、两 rank 数。
  - `msu-pe-phase-d-04-stress-strain` 读 `stress_strain_lammps_<N>r.csv` 与 `stress_strain_atomiverse_<N>r.csv`;其 ac-002 是“载重门”,显式以本环的 `stress_strain_reference.csv` 落地为前提。
- **上游伞型依赖(显式,非本环职责)**:`msu-pe-validation.md` Phase **G.2**(per-axis `Pxx/Pyy/Pzz` 真实虚张应力路径)修复 `run_msu_atomiverse.sh` 当前用各向同性 P 伪造 `Pxx/Pyy/Pzz` 列的缺陷。**本环不修复该缺陷**。由此产生一个明确的产物分期决定(见 Design § 产物分期):本环现在产出 **LAMMPS** 的 4 个 CSV,Atomiverse 的 4 个 CSV 在伞型 Phase G.2 落地后由操作员补跑;若在 G.2 落地前运行 Atomiverse 参考,其 `Pxx/Pyy/Pzz` 列仍是各向同性占位值,下游消费 Atomiverse per-axis 应力的检查必须知悉这一点。

## Domain basis

参考运行的物理是不可变的上游 deck `lammps_msu_paper.in`(Tschopp 2010;M.A. Tschopp, M.F. Horstemeyer, M.K. Hossain et al., "Molecular dynamics simulations of deformation in amorphous polyethylene," *Polymer* 2010;CAVS ICME 来源 <https://www.cavs.msstate.edu/icme/code/lammps/atomistic_research/polymers/amorphous.php>)。体系:10 链 × 1000 个 united-atom CH₂ = 10000 原子,立方盒 80.0586 Å,每个 UA 珠质量 14.02 amu(`SOURCES.md`)。

力场——Dreiding-UA,逐字取自 deck:harmonic bond *k* = 350 kcal/mol/Å²、*r*₀ = 1.53 Å;harmonic angle *k* = 60 kcal/mol/rad²、θ₀ = 109.5°;multi/harmonic dihedral *A*₁..*A*₅ = 1.73, −4.49, 0.776, 6.99, 0.0 kcal/mol;LJ ε = 0.112 kcal/mol、σ = 4.01 Å、*r*_cut = 10.5 Å。单位:LAMMPS `units real` 即 kcal/mol;Atomiverse 用 eV(`1 kcal/mol = 1/23.060547830619026 eV`)。

五阶段协议(每个 fix、ensemble、timestep、步数都是逐字保真的基准真相):

| 阶段 | Ensemble / fix | T | P | 步数 | dt |
|---|---|---|---|---|---|
| 1 | `nve/limit 0.05` + `langevin 500 500 10.0 904297` | 500 K | — | 10000 | 1.0 fs |
| 2 | `npt ... iso` + `momentum 1 linear 1 1 1` | 500 K | 0 bar | 50000 | 0.5 fs |
| 3 | `npt temp 500 100 ... iso` + `momentum`(降温斜坡) | 500→100 K | 0 bar | 50000 | 0.5 fs |
| 4 | `npt ... iso` + `momentum` | 100 K | 0 bar | 50000 | 0.5 fs |
| 5 | `npt ... y ... z`(各向异性) + `deform 1 x erate 1e-5` | 100 K | 0 bar(y,z) | 171800 | 1.0 fs |

阶段 1 是载重的:`msu_pe.data` 是已构建但**未弛豫**的构型,存在非键重叠(初始 `E_pair ≈ 3.4×10⁸ kcal/mol`),位移受限的 `nve/limit` 积分器是强制的,普通 Velocity-Verlet 会让体系炸开。阶段 1 的 `fix langevin` 是**带种子的**(`904297`),因此是**随机**热浴;两引擎用不同 RNG 流(CUDA 还是单精度,但不在本 CPU 臂范围)永不共享轨迹,只采样同一系综——这是为何下游所有跨引擎对比都是统计分块平均(`msu-pe-align-03`,Flyvbjerg & Petersen 1989, *J. Chem. Phys.* 91, 461, doi:10.1063/1.457480)。本环只**产出**参考数据,不做对比,但归一化必须保持每列单位一致,使下游分块平均成立。

总步数 = 10000 + 50000 + 50000 + 50000 + 171800 = **331800 步**。`lammps.in` 头部注释写 ~330000,`SOURCES.md` § Phase D.1 写 ~360000——审计任务须以实际步数 331800 统一这两处描述(`SOURCES.md` 的 ~360000 是过估,须更正为与 deck 一致)。无论哪个数字,全协议都是 ~24 小时级墙钟,落在 PDC 24h 作业上限内,因此参考运行是操作员动作。

应力-应变:工程应变 ε = box_x/box_x0 − 1(无量纲);σ_xx 以 GPa 报告。LAMMPS deck 的 `variable sigxx equal "-pxx*1.01325e-4"` 把 atm 转 GPa 并取负号(拉伸为正)。规范 `stress_strain` 模式两列 `strain,sigma_xx`,`sigma_xx` 单位 GPa。

## Design

本规范的所有改动局限于 `tools/`(两个运行脚本 + 一个新建归一化脚本 + 其测试)与 `benchmarks/md/reference/polyethylene/`(deck 审计 + 产物落地)。**无生产 C++ 改动**。

### deck 保真度审计

`lammps.in` 在 `SOURCES.md` § Phase D.1 中声称其物理与 `lammps_msu_paper.in` 逐字相同,只在非物理输出钩子上不同。本环把这一声称从“散文断言”升级为**可机验的门**:一个脚本(`tools/audit_lammps_deck.sh`,新建)对两个 deck 的物理相关行做规范化 diff——剥离 `lammps.in` 独有的非物理钩子(`fix print` / `fix thermolog` / `fix ss5` / `fix rdfavg` / `dump traj` / `compute myrdf` / `variable p*` / `thermo_style` / `thermo` 频率行)、剥离 `lammps_msu_paper.in` 独有的非物理钩子(`log`、`write_restart`、`reset_timestep`、`variable simname/fname`、其 `def1`/`def2` 的 `fix print`、`variable p*`/`t*`、`thermo_style`/`thermo`),然后断言剩余的物理行——`units`/`boundary`/`atom_style`/`read_data`、`neighbor`/`neigh_modify`、所有 `bond_style`/`bond_coeff`/`angle_*`/`dihedral_*`/`pair_*`、每阶段的 `fix 1`/`fix 2`/`velocity`/`timestep`/`run`/`unfix`——逐行集合相等。脚本退出码 0 表示物理保真,非 0 表示发现偏差并打印具体行。该审计可在无集群环境下运行,是 `code` 类标准。

审计若发现真实偏差(例如 timestep 或 step count 不符),修复方向是改 `lammps.in`(它是可改的派生 deck),绝不改 `lammps_msu_paper.in`(不可变基准真相)。`reset_timestep` 的有无是非物理的:`lammps.in` 故意不发 `reset_timestep` 以让 step 列全局单调,这是设计意图,审计须把 `reset_timestep` 归入非物理钩子而非偏差。

### 运行脚本最终化

- `run_msu_lammps.sh` 当前以 `${LMP_BIN:-lmp}` 启动通用 `lmp`。改为默认 `lmp_cuda` 并以 **CPU 模式**运行(用户构建的单一 `lmp_cuda` 二进制兼跑 CPU 与 GPU 模式;CPU 臂用 CPU 模式——即不传 GPU package 参数,纯 MPI/OpenMP)。`LMP_BIN` 环境变量仍可覆盖。
- 两个脚本都加 PDC Dardel SLURM 提交路径。引入一个共享提交 wrapper `tools/run_msu_pdc.sh`(新建):它接收引擎(`lammps`|`atomiverse`)、rank 数、PDC 账户(`-A <project>`,必填,无默认——无 fallback;CPU 臂的具体值是 `naiss2025-1-25`,GPU 用 `naiss2025-1-25-gpu`,但脚本 CLI 不硬编码、由调用方传入)、分区(默认 `main`,CPU 臂只用 `main`)、墙钟上限(默认 24:00:00,PDC 上限),生成并 `sbatch` 一个 SLURM 作业脚本,作业体内以 `srun -n <RANKS>` 调用对应的 `run_msu_lammps.sh` / `run_msu_atomiverse.sh`。两个运行脚本保留其现有的本地 `mpirun`/`srun` 直跑路径(用于登录节点冒烟测试),`run_msu_pdc.sh` 是其上的 PDC 提交外壳。`run_msu_pdc.sh` 不读不写物理产物,只做提交编排。
- 注意:`run_msu_atomiverse.sh` 当前用各向同性 P 伪造 `Pxx/Pyy/Pzz`(`awk` 的 `$5","$5","$5`)。**本环不修该缺陷**——它是伞型 Phase G.2。本环只在脚本注释里把该 stand-in 标注为 G.2 待修,并在产物分期决定中记录其后果。

### 归一化(canonicalization)

“归一化”= 把原始引擎输出转换为下游消费的规范模式。规范模式:thermo 为 `step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz`;stress-strain 为 `strain,sigma_xx`。

- **LAMMPS 在线路径**:`lammps.in` 的 `fix thermolog` 与 `fix ss5` 已直接以规范模式流式输出 `thermo.lammps.csv` 与 `stress_strain.lammps.csv`;`run_msu_lammps.sh` 已把它们 `cp` 成 `thermo_lammps_<N>r.csv` / `stress_strain_lammps_<N>r.csv`。这条路径已正确,本环只需审计确认其 header 与列序与规范模式一致。
- **LAMMPS 离线路径(已存在的部分运行)**:`runs/lammps_n8_paper/` 已有一次完成的 8-rank `lammps_msu_paper.in` paper-deck 运行的原始阶段-5 输出——`PE_nc10_cl1000.def1.txt`(空格分隔,列 `strain Pxx_GPa Pyy_GPa Pzz_GPa lx ly lz temp`)与 `.def2.txt`(列 `strain epair ebond eangle edihed`)。本环新建 `tools/canonicalize_msu_output.py`,把 `def1.txt` 转换为规范 `stress_strain_lammps_8r.csv`(`strain` = def1 第 1 列,`sigma_xx` = def1 第 2 列即已转 GPa 并已取拉伸为正号的 Pxx)。这把一份既有的真实参考运行立即纳入规范模式,无需重跑该 8-rank LAMMPS 运行。
- **Atomiverse 路径**:benchmark 产出 `thermo.csv` 与 `stress_strain.csv`;`run_msu_atomiverse.sh` 已用内联 `awk` 重写成规范模式(thermo 的 `Pxx/Pyy/Pzz` 是各向同性占位,见上)。本环把这段内联 `awk` 重构进 `canonicalize_msu_output.py` 的一个 `atomiverse` 子命令,使两引擎共用一个、有测试覆盖的归一化实现,降低漂移风险;`run_msu_atomiverse.sh` 改为调用该脚本。
- `canonicalize_msu_output.py` 是纯 stdlib(`csv`/`argparse`)、无副作用、幂等的转换器,接口为 `canonicalize_msu_output.py <engine> <kind> <in> <out>`,`engine ∈ {lammps,atomiverse}`、`kind ∈ {thermo,stress-strain}`。它有 pytest 单元测试,用小型合成 fixture 验证列映射、header、空格-vs-逗号分隔的解析,以及对 `def1.txt`/`def2.txt` 实际格式的解析。

### 产物分期决定(显式)

因为 Atomiverse 的 per-axis P 在伞型 Phase G.2 落地前仍是伪造的,本环采取分期落地:

- **现在(本环代码落地即可)**:`stress_strain_lammps_8r.csv` 由 `runs/lammps_n8_paper/` 既有运行经 `canonicalize_msu_output.py` 产出并提交——这是 `msu-pe-phase-d-04` ac-002 载重门所需的 `stress_strain_reference.csv` 的 LAMMPS 侧首件实数据。
- **操作员运行后**:`thermo_lammps_<N>r.csv`、`stress_strain_lammps_<N>r.csv`(N ∈ {1,8})由四次 PDC 全协议运行中的 LAMMPS 两次产出并提交。
- **Phase G.2 落地后**:`thermo_atomiverse_<N>r.csv`、`stress_strain_atomiverse_<N>r.csv` 由 Atomiverse 两次运行产出并提交,带真实 per-axis P。若操作员在 G.2 前已跑 Atomiverse 参考,其产物可作占位提交,但须在 `SOURCES.md` 中标注 `Pxx/Pyy/Pzz` 为各向同性 stand-in。

### 文档一致性

`SOURCES.md` § Phase D.1 与 `README.md` 的描述须与本环最终化的脚本与文件名方案保持一致:`lmp_cuda` CPU 模式、PDC `sbatch` 程序、331800 步的实际步数(更正 `SOURCES.md` 的 ~360000 过估)、归一化脚本的存在。本环编辑这两个 Markdown 文件使其描述与代码一致。

## Files to create or modify

- `benchmarks/md/reference/polyethylene/lammps.in` — 仅在审计发现真实物理偏差时修改;否则不动(审计预期为绿)。
- `benchmarks/md/reference/polyethylene/SOURCES.md` — 更新 § Phase D.1:更正总步数为 331800、记录 `lmp_cuda` CPU 模式与 PDC `sbatch` 程序、记录归一化脚本与产物分期决定。
- `benchmarks/md/reference/polyethylene/README.md` — 更新 Phase D 节:`lmp_cuda` CPU 模式、PDC `main` 分区 `sbatch` 复现命令、归一化脚本。
- `tools/audit_lammps_deck.sh` (new) — `lammps.in` ≡ `lammps_msu_paper.in` 物理保真度审计门。
- `tools/canonicalize_msu_output.py` (new) — 原始引擎输出 → 规范模式 CSV 的转换器(LAMMPS `def1`/`def2`/`fix print` 与 Atomiverse `thermo.csv`/`stress_strain.csv`)。
- `tools/test_canonicalize_msu_output.py` (new) — `canonicalize_msu_output.py` 的 pytest 单元测试,含合成 fixture。
- `tools/run_msu_lammps.sh` — 默认 `lmp_cuda` CPU 模式;归一化改为调用 `canonicalize_msu_output.py`。
- `tools/run_msu_atomiverse.sh` — 归一化改为调用 `canonicalize_msu_output.py`(内联 `awk` 移除);保留 G.2 待修的伪造-P 注释。
- `tools/run_msu_pdc.sh` (new) — PDC Dardel SLURM 提交 wrapper(`sbatch`/`srun`、`-A <project>`、`main` 分区、24h 上限),包裹两个运行脚本。
- `benchmarks/md/reference/polyethylene/stress_strain_lammps_8r.csv` (new) — 由 `runs/lammps_n8_paper/PE_nc10_cl1000.def1.txt` 经归一化产出并提交(产物分期“现在”项)。

## Tasks

- [x] Write failing tests for canonicalize_msu_output column mapping and parsing (tools/test_canonicalize_msu_output.py)
- [x] Implement canonicalize_msu_output.py with lammps/atomiverse × thermo/stress-strain conversion (tools/canonicalize_msu_output.py)
- [x] Implement audit_lammps_deck.sh asserting lammps.in physics equals lammps_msu_paper.in (tools/audit_lammps_deck.sh)
- [x] Switch run_msu_lammps.sh to lmp_cuda CPU mode and route output through canonicalize_msu_output.py (tools/run_msu_lammps.sh)
- [x] Route run_msu_atomiverse.sh canonicalization through canonicalize_msu_output.py, keeping the G.2 faked-P note (tools/run_msu_atomiverse.sh)
- [x] Implement run_msu_pdc.sh PDC Dardel sbatch submission wrapper for both run scripts (tools/run_msu_pdc.sh)
- [x] Generate stress_strain_lammps_8r.csv from runs/lammps_n8_paper/PE_nc10_cl1000.def1.txt via canonicalize_msu_output.py
- [x] Update SOURCES.md and README.md for lmp_cuda CPU mode, PDC sbatch, 331800-step count, and the canonicalization script
- [x] Run full check + test suite

## Testing strategy

- **Happy path** — `canonicalize_msu_output.py lammps stress-strain` on a synthetic `def1.txt`-format fixture emits a two-column `strain,sigma_xx` CSV with the correct header and column mapping; `canonicalize_msu_output.py atomiverse thermo` on a synthetic benchmark `thermo.csv` emits the eight-column canonical schema. `audit_lammps_deck.sh` on the in-tree `lammps.in` / `lammps_msu_paper.in` exits 0. A shortened-protocol smoke run (`run_msu_lammps.sh` / `run_msu_atomiverse.sh` with `-- run` step counts reduced via the extra-args path on the login node) produces correctly-schema'd `thermo_*` and `stress_strain_*` CSVs.
- **Edge cases** — `canonicalize_msu_output.py` on an empty or header-only input emits a header-only output, not a crash; on a malformed row (wrong column count) it fails loudly naming the bad line, not silently. `audit_lammps_deck.sh` flags a deliberately mutated `lammps.in` (e.g. a changed `timestep`) as a non-zero physics deviation. `run_msu_pdc.sh` invoked without `-A <project>` exits non-zero with a clear required-flag error (no account fallback). `run_msu_lammps.sh` with `lmp_cuda` absent from PATH and no `LMP_BIN` override exits with the existing binary-not-found code.
- **Domain validation** — `canonicalize_msu_output.py` on the real `runs/lammps_n8_paper/PE_nc10_cl1000.def1.txt` produces a `stress_strain_lammps_8r.csv` whose first `strain` value is ≈ 0 (step-0 row) and whose `strain` column increases monotonically toward ε ≈ 1.0, with `sigma_xx` in GPa (tension-positive sign convention from the deck's `-pxx*1.01325e-4`). The committed `stress_strain_lammps_8r.csv` is consumable by `msu-pe-phase-d-04`'s stress-strain check. The deck-fidelity audit confirms every fix / ensemble / timestep / step count / FF coefficient of `lammps.in` matches `lammps_msu_paper.in` so the reference runs sample the immutable Tschopp 2010 protocol.
- **Operator runtime gates** (deferred, parked at code-complete) — the four full reference runs are ~24h-class PDC Dardel operator actions, not CI: `run_msu_pdc.sh lammps 1`, `run_msu_pdc.sh lammps 8`, `run_msu_pdc.sh atomiverse 1`, `run_msu_pdc.sh atomiverse 8` on the `main` partition with `-A <project>`. On completion the operator commits the eight canonical CSVs and the stage-4 thermo / RDF / trajectory artifacts for both engines at both rank counts.

## Out of scope

- 修复 `run_msu_atomiverse.sh` 伪造 per-axis `Pxx/Pyy/Pzz` 的缺陷——这是伞型 `msu-pe-validation.md` Phase **G.2**(真实虚张应力路径)。本环只标注该依赖并据此分期产物。(考虑过在本环顺手修;否决——它需要改生产 C++ benchmark,越出本环 `tools/`-only 范围,且伞型已专门立项。)
- `fix momentum` 等价物(Phase G.1)与能量分解 `E_pair/E_bond/E_angle/E_dihedral`(Phase G.3)——同属伞型 Phase G 的生产 C++ 修复,本环不碰。
- 任何对 LAMMPS 与 Atomiverse 引擎结果的**对比/校验**——那是下游 `msu-pe-align-03` / `msu-pe-phase-d-03` / `msu-pe-phase-d-04` 的职责。本环只产出参考数据,不下结论。(考虑过在本环加一个轻量一致性检查;否决——会与下游链环的对比逻辑重复并混淆“产数据”与“做对比”。)
- 任何对 `lammps_msu_paper.in` 或 `msu_pe.data` 的修改——它们是不可变基准真相;保真度只能单向 deck → `lammps.in`。
- CUDA / LAMMPS-GPU 参考运行——这是伞型 Phase **I**(I.2/I.4),在 `gpugh` 分区,显式不在本 CPU 臂范围。
- 生产 C++ 代码、CMake、CI 配置——本规范完全不触碰。

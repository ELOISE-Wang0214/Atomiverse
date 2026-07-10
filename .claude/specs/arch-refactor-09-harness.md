---
title: arch-refactor-09-harness — 文档收口：CLAUDE.md 架构章重写 + 术语漂移修复 + canon/notes 同步
status: approved
date: 2026-07-10
chain: arch-refactor (09/09)
depends_on: [arch-refactor-01-bindings, arch-refactor-02-halo, arch-refactor-03-neighbor, arch-refactor-04-engines, arch-refactor-05-context, arch-refactor-06-frame, arch-refactor-07-molrs, arch-refactor-08-python]
---

# arch-refactor-09-harness — 文档收口：CLAUDE.md 架构章重写 + 术语漂移修复 + canon/notes 同步

## Summary

本篇是 arch-refactor 链（9 篇）的链尾，全部为文档工作：在 01–08 的代码变更全部落地后，把仓库的"权威架构陈述"一次性对齐到最终态。核心交付是 CLAUDE.md `## Architecture` 章的全量重写（三支柱/三后端、Driver 五类 taxonomy 与 Gaussian job type 类比、方法接入模型"方法=inner Solver+Evaluator、物性=组件 Capability"、Component 核心基类集+opt-in 协议表+PhaseKind 绑定规则、Context service 分组表、Frame 多后端设计规范章、molrs 边界、校正后的 conventions），并补录长期缺失的现状：SamplingDriver+`cpu::mc::MC` 已实现却在 driver 表缺席、还被写成假设；同时修复 manual-discovery §9 记录的术语漂移残留（`component_role.h` "3 categories" 过期注释、`mainpage.dox` "three-category taxonomy"、"fourth category" 计数漂移等）、把本链新术语（PhaseKind、HaloParticipant、Neighbor 协议、service slice）补进 manual-system.md §2 术语 canon（Driver 作为 loop-owner canon 名已由用户裁决保留），最后按 notes→CLAUDE.md 知识流在 `.claude/notes/notes.md` 写一条链归档条目。完成后，任何 agent 或读者从 CLAUDE.md 读到的架构与 `dev` 分支代码事实逐条一致。

## Domain basis

不适用（`science.required: true` 但本篇不声明新物理）：Gaussian job type 类比（Evaluate=SP、Optimize=Opt、Response=Freq/振动、Dynamics=MD、Sampling=MC）是命名学定位，不引入方程或数值目标；所有科学内容（Metropolis 判据、FD Hessian 等）已在各自实现 spec 中给出依据，本篇只做文档陈述与代码事实的一致性校对。

## Design

不新增任何代码符号；重写后的 CLAUDE.md `## Architecture` 章结构（对应 scope a–h）：

1. **三支柱 + 三后端**（a）——Driver/Component/Frame 定义 + `cpu/`/`cuda/`/`torch/` 后端分层；目录表刷新至磁盘现状（补 `mc/`、`semiempirical/` 两个 domain 行，逐行核对 `include/atomiverse/<backend>/` 一级目录）。
2. **Driver 五类 taxonomy**（b）——按终止语义分类的 category base 表扩为五行，加 Gaussian job type 类比列；补录 `cpu::SamplingDriver`（`sampling::Phase`，phases `propose→energy→end_of_step`，inline Metropolis，`cpu/driver/sampling_driver.h:62`）与具体驱动 `cpu::mc::MC`（`cpu/mc/mc.h:52`）；删除 :239-242 把 Sampling 写成假设的句子，扩展指引改用尚未实现的类别（如 Path/NEB）举例；`## Implemented components → ### Drivers` 表（:452-461）加 Sampling 行，`### MD components` 旁补 MC 组件段（TranslateMove/Mover role），Key user-facing APIs 补一个 `cpu::mc::MC` snippet（与 `mc.h:41-51` @code 用法一致）。
3. **方法接入模型**（c）——用户钦定定位写透："方法（Gaussian SCF、PW、GPW、NDDO-AM1）= inner Solver（纯 kernel、非 Driver 子类）+ Evaluator（组件 wrapper）；物性（电荷/偶极/Mulliken/AM1-BCC 等）= 组件 Capability，永不出现在 Driver 层"。04 落地后 `GPWSolver` 为真实类，表述随之转实（当前 `cpu/gpw/` 尚无该类——重写时以 04 落地态为准）；`cpu/semiempirical/` NDDO/AM1 纳入方法行。
4. **Component 模型**（d）——01/02/03 后的 ComponentBase 核心基类最终集 + opt-in 协议表（`HaloParticipant`/`NeighborConsumer`/`NeighborProvider`/`PairEvaluator`/`AsyncComponentBase`）+ PhaseKind 绑定规则（跨 workflow 组件不声明 `default_bindings()`，强制 `.on(phase).add()`）；role 表补 `Mover`。
5. **Context 分层表**（e）——05 后的 core 字段 + `ctx.md`/`ctx.mc`/`ctx.opt`/`ctx.analysis` service slice，per-backend 差异。
6. **Frame 多后端设计规范章**（f）——per-backend 表示表（cpu per-atom LAMMPS 布局 / cuda 单卡 flat + 多卡 per-atom 双拓扑 / torch tensor；molrs flat Block 为边界 schema）、转换点地图（`materialize`/`write_back`/`to_cuda`/`to_cpu`/`torch::transfer` 职责）、write_back round-trip 契约表（06）、精度 tier 映射（`StorageReal`/`WorkReal`/`AccReal`、`TagInt` 等 ↔ cuda `real_t`/`work_t`/`prec_t`）、Frame purity 规则；**删除 :317-319 "identical field set" 不实表述**（CUDA 双拓扑现实）。
7. **molrs 边界**（g）——硬依赖、`molrs.hpp` 唯一 CXX include、engine-measures-molrs-analyzes、streaming accumulators；07 后 cuda RDF 归 molrs 的最终态。
8. **conventions**（h）——继承现文并校正（driver 模板句列出五个 category base；phase 命名补 Sampling 行 `propose, energy, end_of_step`）。

伴随修复（不属重写正文）：`facade/component_role.h:17-24` 注释扩为五类校验规则（Response：≥1 个 `Capability::Force` Evaluator、固定 6N 位移预算；Sampling：≥1 Mover + ≥1 `Capability::Energy`，与 `sampling_driver.h:32-33` 一致）；`cpu_driver.h:20`/`cuda_driver.cuh:20` "three category bases" 计数校正（CPU 五类、CUDA 按实际）；`mainpage.dox:96-98` "same three-category taxonomy" 改为按后端如实陈述；`ch14-vibrations.md:81,237`、`cpu-internals.md:27` "fourth category" 计数漂移更新；`ch10-classical-md.md:54` 表头 `Capabilities`→`Capability`；manual-discovery §9 其余条目（1/2/4/5/8，docs 现已无命中）以 grep 审计收口。`manual-system.md` §2 canon 表：loop-owner 行补 `SamplingDriver`，新增 PhaseKind/HaloParticipant/NeighborConsumer+NeighborProvider/Context service slice 行。`.claude/specs/manual-discovery.md` 是历史快照（"no rewrites"），**不修改**，所有 grep 门排除 `.claude/`。

## Files to create or modify

- CLAUDE.md
- include/atomiverse/facade/component_role.h
- include/atomiverse/cpu/driver/cpu_driver.h
- include/atomiverse/cuda/driver/cuda_driver.cuh
- docs/mainpage.dox
- docs/developer/architecture.md
- docs/developer/cpu-internals.md
- docs/tutorial/ch14-vibrations.md
- docs/tutorial/ch10-classical-md.md
- .claude/specs/manual-system.md
- .claude/notes/notes.md

## Tasks

- [ ] Rewrite CLAUDE.md `## Architecture` 章（scope a/c/d/e/f/g/h：三支柱+目录表刷新含 `mc/`+`semiempirical/`、方法接入模型 SCF/PW/GPW/NDDO-AM1=Solver+Evaluator 且物性=Capability、Component 核心集+五协议 opt-in 表+PhaseKind 强制 `.on()` 规则、Context `ctx.md/mc/opt/analysis` 分层表、Frame 多后端规范章并删除 "identical field set" 表述、molrs 边界、conventions 校正；顺手校正 frontmatter `notes_path` → `.claude/notes/notes.md`）(CLAUDE.md)
- [ ] Update CLAUDE.md Driver 现状补录（scope b：taxonomy 表扩五行加 Gaussian job type 类比列、删 :239-242 Sampling 假设句改用 Path/NEB 举例、`### Drivers` 表加 `cpu::mc::MC` 行、`### MD components` 旁补 MC 组件段、Key APIs 补 MC snippet）(CLAUDE.md)
- [ ] Fix 头文件 doc 注释计数漂移：`component_role.h:17-24` 校验列表扩为五类（补 ResponseDriver + SamplingDriver 规则，与 `sampling_driver.h:32-33` 一致）；`cpu_driver.h:20` 与 `cuda_driver.cuh:20` "three category bases" 改为如实计数 (include/atomiverse/facade/component_role.h, include/atomiverse/cpu/driver/cpu_driver.h, include/atomiverse/cuda/driver/cuda_driver.cuh)
- [ ] Fix Doxygen mainpage：`mainpage.dox:96-98` "same three-category taxonomy" 改为按后端如实陈述 CPU 五类/CUDA 覆盖，并在 workflows 节补 Sampling/MC 一句 (docs/mainpage.dox)
- [ ] Audit and fix manual-discovery §9 全部 8 条在 docs/ 的残留：`ch14-vibrations.md:81,237` 与 `cpu-internals.md:27` "fourth category" 计数更新、`ch10-classical-md.md:54` `Capabilities`→`Capability`、`architecture.md` 五类陈述核对；grep 确认条目 1/2/4/5/8（`cpu::observer::RDF`/`cpu::recorder::XYZ`/`cpu::Trajectory`/`scf::SCF` driver/`scf::Phase`）在 docs/+examples/+README.md 均无命中，有残留即修 (docs/developer/architecture.md, docs/developer/cpu-internals.md, docs/tutorial/ch14-vibrations.md, docs/tutorial/ch10-classical-md.md)
- [ ] Extend manual-system.md §2 terminology canon：loop-owner 行补 `SamplingDriver`；新增 `PhaseKind`、`HaloParticipant`、`NeighborConsumer`/`NeighborProvider`、Context "service slice"（`ctx.md/mc/opt/analysis`）四行 (.claude/specs/manual-system.md)
- [ ] Append notes 归档条目：`## [2026-07-10] arch-refactor 链收口` — 链摘要（01–08 一句一条）+ 关键决策（保留 Driver 命名、PhaseKind 判别、opt-in 协议表、Context service 分组、write_back 契约、cuda RDF 归 molrs）(.claude/notes/notes.md)
- [ ] Verify CLAUDE.md 全部代码 snippet 有编译实体对应（仓库无 snippet 编译门禁，按替代验证：MD↔`examples/07_md_lj.cpp`、geom-opt↔`examples/05_geom_opt.cpp`、PW↔`examples/08_pw.cpp`、新 MC snippet↔`tests/cpu/mc/` 等；以 `-DATV_BUILD_EXAMPLES=ON -DATV_BUILD_TESTS=ON` 构建确认对应文件编译通过且 API 调用形式一致）(CLAUDE.md)
- [ ] Run full check + test suite（`cmake --build build -j$(nproc)` + `ctest --test-dir build --output-on-failure`，确认头文件注释改动零行为影响）

## Testing strategy

- **Happy path（docs 一致性，全部 grep 可观察）**：taxonomy 表五行且含 `cpu::mc::MC`；目录表逐行对应 `include/atomiverse/<backend>/` 磁盘一级目录；`component_role.h` 注释五类校验；`grep -rn "three-category" docs/ include/ CLAUDE.md` 与 `grep -rni "fourth categor" docs/` 归零；`grep -c "identical field set" CLAUDE.md` = 0；§9 条目 grep 集（`cpu::observer::RDF|cpu::recorder::XYZ|cpu::Trajectory|scf::SCF|Capabilities`）在 docs/+examples/+README.md 归零。
- **Edge cases**：grep 门一律排除 `.claude/`（manual-discovery.md 为历史快照、`.claude/NOTES.md` 为 legacy 镜像，两者刻意保留旧措辞）；`Capabilities` 匹配限定为类型/表头用法，不误伤普通英文复数；`torch` 仍按 manual-system §2 canon 表述为 differentiable-potential library 而非第三 backend。
- **Runtime（唯一非 docs 判据）**：`component_role.h`/`cpu_driver.h`/`cuda_driver.cuh` 仅注释变更，全量构建 + ctest 必须与基线同绿；`-DATV_BUILD_EXAMPLES=ON` 下 snippet 对应的 examples/tests 目标编译通过。
- **Domain validation**：本篇无新物理；科学正确性由既有测试套件（MC NIST parity 等）回归覆盖即可。

## Out of scope

- manual-system.md §3 的 16 章教程本体撰写/改写（本篇只动其 §2 canon 表）。
- 任何代码行为变更——01–08 已全部落地；本篇对 `.h/.cuh` 的改动仅限 doc 注释。
- `.claude/specs/manual-discovery.md` 的修改（历史快照，"no rewrites"）与 `.claude/NOTES.md` legacy 镜像清理。
- Python 绑定文档与假可用面（08-python 已收口）；docs/api Doxygen 逐符号刷新（`docs` target 自动生成）。
- 新增 snippet 编译 CI 门禁（考虑过：为 CLAUDE.md 片段建 extract-and-compile 测试——收益/维护比不划算，以 examples/tests 对应验证替代；若未来漂移复发再立独立 spec）。
- benchmarks/README 与 benchmark 体系文档（归 manual-system Task #7）。

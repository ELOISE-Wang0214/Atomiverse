---
spec: arch-refactor-09-harness
created: 2026-07-10
criteria:
  - id: ac-001
    summary: "CLAUDE.md Driver taxonomy 表为五行且含 Gaussian job type 类比"
    type: docs
    evaluator_hint: "grep CLAUDE.md '## Architecture' driver taxonomy table"
    pass_when: "CLAUDE.md 的 Driver taxonomy 表恰含五个 category base 行（EvaluateDriver/DynamicsDriver/OptimizeDriver/ResponseDriver/SamplingDriver），Sampling 行同时含 sampling::Phase 与 cpu::mc::MC，且五个 Gaussian job type 类比对（Evaluate=SP、Optimize=Opt、Response=Freq、Dynamics=MD、Sampling=MC）全部出现在该节。"
    status: pending
  - id: ac-002
    summary: "Sampling 不再被表述为假设；Implemented Drivers 表含 MC 行"
    type: docs
    evaluator_hint: "grep -n 'Sampling for MC/HMC/REMD' CLAUDE.md"
    pass_when: "`grep -n 'e.g. Sampling for MC/HMC/REMD' CLAUDE.md` 无命中（扩展指引改用未实现类别如 Path/NEB 举例），且 `## Implemented components` 的 Drivers 表含一行同时出现 cpu::mc::MC 与 sampling::Phase。"
    status: pending
  - id: ac-003
    summary: "后端目录表与磁盘一级目录一致（含 mc/ 与 semiempirical/）"
    type: docs
    pass_when: "CLAUDE.md 目录分层表含 mc/ 与 semiempirical/ 两行，且表中每一行目录名都存在于 include/atomiverse/cpu/（或注明的对应后端）之下（ls include/atomiverse/cpu/ 逐行核对无缺失、无幽灵行）。"
    status: pending
  - id: ac-004
    summary: "方法接入模型成章：方法=Solver+Evaluator，物性=Capability"
    type: docs
    pass_when: "CLAUDE.md Architecture 章含方法接入小节，列出 SCF、PW、GPW、NDDO-AM1 四个方法行（各标注 inner Solver 类 + Evaluator 组件），并含明确陈述'物性（电荷/偶极等）是组件 Capability、永不出现在 Driver 层'；全文无任何暗示 GPWSolver 不存在/是占位的措辞。"
    status: pending
  - id: ac-005
    summary: "Component 模型含五协议 opt-in 表 + PhaseKind 绑定规则"
    type: docs
    pass_when: "CLAUDE.md 的 Component 小节含 opt-in 协议表，五个协议名 HaloParticipant、NeighborConsumer、NeighborProvider、PairEvaluator、AsyncComponentBase 全部出现，并含规则句：跨 workflow 组件不声明 default_bindings，必须经 .on(phase).add() 显式绑定；role 列表含 Mover。"
    status: pending
  - id: ac-006
    summary: "Context 分层表列出 core + 四个 service slice"
    type: docs
    pass_when: "CLAUDE.md 含 Context 分层表，ctx.md、ctx.mc、ctx.opt、ctx.analysis 四个 slice 名与 core 字段行全部出现，且注明 per-backend 差异。"
    status: pending
  - id: ac-007
    summary: "Frame 多后端规范章成立且 identical-field-set 表述删除"
    type: docs
    evaluator_hint: "grep -c 'identical field set' CLAUDE.md"
    pass_when: "`grep -c 'identical field set' CLAUDE.md` 为 0；Frame 章含四要素：per-backend 表示表（cpu per-atom LAMMPS / cuda 单卡 flat + 多卡双拓扑 / torch tensor / molrs flat Block 边界 schema）、转换点地图（materialize/write_back/to_cuda/to_cpu/torch::transfer）、write_back round-trip 契约表、精度映射（StorageReal/WorkReal/AccReal ↔ cuda real_t/work_t/prec_t）。"
    status: pending
  - id: ac-008
    summary: "component_role.h 注释列出全部五类 driver 的校验规则"
    type: docs
    pass_when: "include/atomiverse/facade/component_role.h 的文件头注释枚举五个 category base 的 build 期校验，其中含 SamplingDriver 行（≥1 Mover + ≥1 Capability::Energy）与 ResponseDriver 行（≥1 Capability::Force evaluator、固定位移预算）。"
    status: pending
  - id: ac-009
    summary: "类别计数漂移归零（three-category / fourth category）"
    type: docs
    evaluator_hint: "grep -rn 'three-category' docs/ include/; grep -rni 'fourth categor' docs/"
    pass_when: "`grep -rn 'three-category' docs/ include/ CLAUDE.md` 与 `grep -rn 'three category bases' include/` 均无命中，且 `grep -rni 'fourth categor' docs/` 无任何把 Response 表述为最新/最后一类的命中（.claude/ 目录除外）。"
    status: pending
  - id: ac-010
    summary: "manual-discovery §9 条目 1/2/4/5/8 残留 grep 归零"
    type: docs
    pass_when: "`grep -rn 'cpu::observer::RDF|cpu::recorder::XYZ|cpu::Trajectory' docs/ examples/ README.md` 无命中；docs/ 中无把 SCF 表述为 scf::SCF driver 或引用 scf::Phase 的命中；docs/ 中作为类型名/表头的 Capabilities（复数）无命中（含 docs/tutorial/ch10-classical-md.md:54 表头已改为 Capability）。"
    status: pending
  - id: ac-011
    summary: "manual-system.md §2 canon 补入本链新术语"
    type: docs
    pass_when: ".claude/specs/manual-system.md §2 表格：loop-owner 行的 category base 列表含 SamplingDriver，且表中新增行分别含 PhaseKind、HaloParticipant、NeighborConsumer/NeighborProvider、service slice（ctx.md/mc/opt/analysis）四个 canon 术语。"
    status: pending
  - id: ac-012
    summary: "notes.md 追加 arch-refactor 链归档条目"
    type: docs
    pass_when: ".claude/notes/notes.md 含一条 2026-07 日期、标题含 arch-refactor 的 ## 条目，正文同时提及：保留 Driver 命名、PhaseKind、opt-in 协议表、Context service 分组、write_back 契约五项决策。"
    status: pending
  - id: ac-013
    summary: "全量构建+测试绿，CLAUDE.md snippet 均有编译实体"
    type: runtime
    evaluator_hint: "cmake --build build && ctest --test-dir build; -DATV_BUILD_EXAMPLES=ON"
    pass_when: "按 CLAUDE.md Build & Test 配置（含 -DATV_BUILD_TESTS=ON -DATV_BUILD_EXAMPLES=ON）全量构建通过、`ctest --test-dir build --output-on-failure` 0 失败；且 CLAUDE.md Key-APIs 各 snippet（含新增 MC snippet）逐一对应一个本次构建中编译成功的 examples/ 或 tests/ 源文件（MC ↔ tests/cpu/mc/），API 调用形式与该文件一致。"
    status: pending
out_of_scope:
  - "manual-system §3 的 16 章教程本体"
  - "任何代码行为变更（.h/.cuh 仅 doc 注释）"
  - "manual-discovery.md 修改（历史快照）与 .claude/NOTES.md legacy 镜像清理"
  - "snippet 编译 CI 门禁（以 examples/tests 对应验证替代）"
---

# Acceptance — arch-refactor-09-harness

ac-001…ac-012 为 docs 型：每条都是"文档陈述 ↔ 代码事实"的单一 grep/对照可观察断言，
第三方无需重读 spec 即可判定。所有涉及旧措辞的 grep 门一律排除 `.claude/`
（manual-discovery.md 历史快照与 NOTES.md legacy 镜像刻意保留原文）。
ac-013 为本篇唯一 runtime 判据：仓库现无 CLAUDE.md snippet 编译门禁，按 spec
Testing strategy 的替代方案以 examples/tests 编译实体验证，同时确认三个头文件的
纯注释改动零行为影响。

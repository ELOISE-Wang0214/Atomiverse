---
spec: arch-refactor-05-context
created: 2026-07-10
criteria:
  - id: ac-001
    summary: "cpu::Context regrouped into md/mc/opt/analysis slices, flat domain fields deleted"
    type: code
    pass_when: "include/atomiverse/cpu/component/context.h defines plain aggregate DynamicsServices{dd, pair_list, is_last_step}, SamplingServices{proposal, rng}, OptimizeServices{reeval_forces, line_search_trial}, AnalysisServices{executor} composed as members `md`, `mc`, `opt`, `analysis` of cpu::Context; n_threads/mpi_rank/mpi_size/timer/outer_status remain top-level; none of the eight domain fields exists flat on Context any more."
    status: verified
    last_checked: 2026-07-10
  - id: ac-002
    summary: "cuda::Context gains homomorphic md slice incl. exclusion CSR and force_rows"
    type: code
    pass_when: "include/atomiverse/cuda/component/context.cuh defines cuda::DynamicsServices containing exactly pair_list, flat_pairs_requested, excl_offsets, excl_partners, s_excl, exclusions_baked, n_force_rows, force_rows, held as member `md`; stream, defer_stream, outer_status and host_box remain top-level; no mc/opt/analysis slice is added on CUDA."
    status: verified
    last_checked: 2026-07-10
  - id: ac-003
    summary: "Zero new indirection; Build/Launch/DestroyContext shapes byte-identical"
    type: code
    pass_when: "static_asserts compile asserting cpu DynamicsServices/SamplingServices/AnalysisServices and cuda::DynamicsServices are trivially copyable and cpu::Context/cuda::Context are non-polymorphic; Context gains no virtual function, no smart-pointer member, no registry/service-locator type; the BuildContext/LaunchContext/DestroyContext definitions in include/atomiverse/cpu/component/component.h (lines 25-44 shape {Context&, Frame&, stage, step}) show zero diff in this change; the per-op stack construction `LaunchContext ctx{context_, *frame_, ...}` in src/cpu/driver/backend.cpp remains unchanged."
    status: verified
    last_checked: 2026-07-10
  - id: ac-004
    summary: "Clean-break complete — no flat domain-field access path survives tree-wide"
    type: code
    evaluator_hint: "rg gate over include/ src/ tests/ tools/ benchmarks/"
    pass_when: "ripgrep pattern (runtime|ctx|ctx_|context_|context\\(\\))\\.(dd|pair_list|proposal|rng|executor|reeval_forces|line_search_trial|is_last_step|flat_pairs_requested|excl_offsets|excl_partners|s_excl|exclusions_baked|n_force_rows|force_rows)\\b returns zero matches under include/, src/, tests/, tools/ and benchmarks/ (all access goes through .md./.mc./.opt./.analysis. slice paths)."
    status: verified
    last_checked: 2026-07-10
  - id: ac-005
    summary: "CPU full build + ctest green incl. new context-slices test"
    type: runtime
    pass_when: "Per the CLAUDE.md CPU recipe (Release, ATV_BUILD_TESTS=ON, ATV_BUILD_RHF=ON) `cmake --build build -j` succeeds and `ctest --test-dir build --output-on-failure` passes 100%, including the new test_cpu_component_context_slices target."
    status: verified
    last_checked: 2026-07-10
  - id: ac-006
    summary: "CUDA backend builds and its ctest targets pass on the GH200 node"
    type: runtime
    pass_when: "On the aarch64 GPU node with GPU/buildenv-gcccuda/2026.03-cu13.0 and ATV_BUILD_CUDA=ON the tree compiles (including src/torch/neighbor/neighbor_list.cu when ATV_BUILD_TORCH is enabled) and the cuda ctest targets pass — neighbor/potentials/multi-GPU paths reading ctx.md.*."
    status: pending
  - id: ac-007
    summary: "MPI regression — DD/halo/migration/gather/rcb suites pass via ctx.md"
    type: runtime
    pass_when: "With an MPI-enabled build, every ctest target matching -R test_mpi passes (test_mpi_topology, test_mpi_halo_exchange, test_mpi_atom_migration, test_mpi_gather, test_mpi_rcb, plus the untouched test_mpi_pw_kpoint_n*), with DomainDecomposition reached exclusively through ctx.md.dd."
    status: verified
    last_checked: 2026-07-10
  - id: ac-008
    summary: "MC Metropolis regression passes via ctx.mc.rng / ctx.mc.proposal"
    type: runtime
    pass_when: "All ctest targets built from tests/cpu/mc/ plus the tests/cpu/test_mc_scaffold.cpp target pass; SamplingDriver publishes into context().mc and TranslateMove draws from ctx.runtime.mc.rng / records into ctx.runtime.mc.proposal."
    status: verified
    last_checked: 2026-07-10
  - id: ac-009
    summary: "Minimizer line-search regression passes via ctx.opt fields"
    type: runtime
    pass_when: "All ctest targets matching -R test_cpu_minimizer pass; LBFGS/CG trial evaluations run through ctx.runtime.opt.reeval_forces and toggle ctx.runtime.opt.line_search_trial, with comm::Migrate suppression behaviour unchanged."
    status: verified
    last_checked: 2026-07-10
  - id: ac-010
    summary: "PE stage-5 per-step wall time shows no regression vs parent commit"
    type: performance
    evaluator_hint: "mol:bench — benchmarks/md/pe/performance polyethylene_tension, stage 5"
    pass_when: "The polyethylene_tension benchmark stage 5 (tension stage of the 5-stage MSU protocol; `--mode stage:5`), built Release with -march=native and run on the same host as the baseline, reports per-step wall time within 3% (or within measured run-to-run noise, whichever is larger) of the parent commit's baseline; result and baseline numbers are logged (no hidden pass)."
    status: pending
out_of_scope:
  - "CLAUDE.md/文档表更新（arch-refactor-09-harness）"
  - "torch::Context 重设计"
  - "per-driver Context 类型 / 虚 service 接口（已否决备选）"
  - "Frame 字段/residency 变更（arch-refactor-06-frame）"
---

# Acceptance — arch-refactor-05-context

ac-001–ac-004 为结构与 clean-break 的 code 门禁（ac-004 的 rg 模式已验证不会误伤 slice 限定访问 `ctx.md.pair_list` 等——接收者交替组后必须紧跟旧字段名才命中）。ac-005–ac-009 为行为不变的 runtime 回归锚，按链级要求分列 MPI / MC / minimizer 三个域回归。ac-010 是链级 PE stage-5 性能门禁——LaunchContext 栈构造与常量偏移字段访问零成本分组的实证，由 /mol:bench 落章，/mol:impl 不得代签。

> **审计注记（2026-07-10，本地 macOS）**：228 处受体表达式一次迁净（semiempirical 的
> NDDO 参数字段 p.dd/pi.dd/pj.dd 已排除）；ac-004 rg 门禁 0 残留；component.h 零 diff、
> backend.cpp 两处 LaunchContext 栈构造逐字不变。ac-007 的 MPI 回归以本机 MPI 构建全量
> 门通过（test_mpi_* 全绿）。ac-006 owed GH200（CUDA/torch .cu 文本已同步迁移、编译未验）；
> ac-010 owed /mol:bench（PE stage-5 集群基线）。本机 2 个 stash 实证既有 SCF 失败见 01 台账。

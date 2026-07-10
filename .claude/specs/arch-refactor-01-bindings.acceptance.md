---
spec: arch-refactor-01-bindings
created: 2026-07-10
criteria:
  - id: ac-001
    summary: "core/binding.h defines PhaseKind + kind-first tagged DefaultBinding"
    type: code
    pass_when: "include/atomiverse/core/binding.h defines `enum class PhaseKind {MD, Minimizer, Evaluate, Vib, Sampling}`, `phase_kind_name(PhaseKind)` returning exactly those five names, and `DefaultBinding` as `{PhaseKind kind; int phase; int order}` with `kind` first and no default member initializer; the header includes nothing from facade/."
    status: verified
    last_checked: 2026-07-10
  - id: ac-002
    summary: "stage.h maps each PhaseEnum to its PhaseKind + make_binding helper"
    type: code
    pass_when: "include/atomiverse/facade/stage.h defines constexpr `kind_of(Phase)` in each of the md / minimizer / evaluate / vib / sampling namespaces returning the matching PhaseKind, plus `ATV::make_binding(PhaseEnum, int order = 0)` producing a DefaultBinding whose kind equals `kind_of(phase)`; its only new include is atomiverse/core/binding.h."
    status: verified
    last_checked: 2026-07-10
  - id: ac-003
    summary: "CPU cross-workflow default-bind add() throws at registration"
    type: runtime
    evaluator_hint: "ctest -R test_cpu_driver_binding_kind"
    pass_when: "test_cpu_driver_binding_kind passes: `cpu::md::MD{}.add(cpu::semiempirical::AM1Evaluator{})` throws std::runtime_error; `cpu::Evaluate{}.add(cpu::md::VelocityVerlet{0.01})` and `cpu::minimizer::Minimize{}.add(cpu::md::VelocityVerlet{0.01})` throw with messages containing both 'MD' and the driver kind name ('Evaluate' / 'Minimizer'); an Evaluate-kind probe component added to cpu::md::MD throws with a message containing both 'Evaluate' and 'MD'. Both lvalue and rvalue add() branches are exercised."
    status: verified
    last_checked: 2026-07-10
  - id: ac-004
    summary: "CUDA cross-workflow default-bind add() throws at registration"
    type: runtime
    evaluator_hint: "ctest -R test_cuda_driver_binding_kind (GH200 GPU node)"
    pass_when: "test_cuda_driver_binding_kind passes on the GPU node: `cuda::md::MD{}.add(cuda::scf::SCFEvaluator{2, 0.1, 1e-6, 100})` throws std::runtime_error, and `cuda::Evaluate{}.add(cuda::md::VelocityVerlet{...})` throws with a message containing both 'MD' and 'Evaluate'."
    status: pending
  - id: ac-005
    summary: "AM1Evaluator and cuda SCFEvaluator no longer declare default bindings"
    type: code
    pass_when: "Neither include/atomiverse/cpu/semiempirical/am1_evaluator.h nor include/atomiverse/cuda/scf/scf_evaluator.cuh contains a `default_bindings` override; their former bare `add()` call sites in tests/cpu/semiempirical/test_am1_scaffold.cpp, tests/cuda/scf/test_scf.cu and tests/cuda/scf/test_rhf.cu use `.on(evaluate::Phase::evaluate).add()` or `single_point`."
    status: verified
    last_checked: 2026-07-10
  - id: ac-006
    summary: "all remaining default_bindings bodies construct via make_binding"
    type: code
    pass_when: "`rg -U 'default_bindings[^}]*static_cast' include/ src/ tests/ tools/ benchmarks/` returns 0 matches (every non-empty default_bindings body uses make_binding), and the full tree compiles — the kind-first aggregate layout makes any un-migrated two-field initializer a compile error, so a successful build is the completeness gate."
    status: verified
    last_checked: 2026-07-10
  - id: ac-007
    summary: "5-driver plan-shape regression unchanged (order + auto-injection)"
    type: runtime
    evaluator_hint: "ctest -R 'test_cpu_driver_binding_kind|test_md_plan_equivalence|test_md_set_api|test_evaluate_driver'"
    pass_when: "Minimal Evaluate / MD / Minimize / Hessian / MC builds in test_cpu_driver_binding_kind assert, via Driver::plan(), the same phase order and auto-injected component sets as before the change; the pre-existing test_md_plan_equivalence, test_md_set_api, test_evaluate_driver, test_cpu_driver_role_validation, vib and mc suites all pass."
    status: verified
    last_checked: 2026-07-10
  - id: ac-008
    summary: "hot path untouched — plan.h and backend.cpp absent from the diff"
    type: code
    pass_when: "`git diff --name-only` for the landing commit contains neither include/atomiverse/facade/plan.h nor src/cpu/driver/backend.cpp; the only facade file with runtime changes is facade/driver.h (registration-time add() validation)."
    status: verified
    last_checked: 2026-07-10
  - id: ac-009
    summary: "full build + default ctest gate green after the one-commit migration"
    type: runtime
    evaluator_hint: "cmake --build build -j && ctest --test-dir build --output-on-failure"
    pass_when: "On the CPU login node the full build succeeds and the default ctest gate passes with zero failures (including test_step_timing and test_timing_report); the CUDA test subset passes on the GPU node. No [[deprecated]] token and no compatibility shim is introduced anywhere in the diff."
    status: pending
out_of_scope:
  - "ComponentBase slimming / HaloParticipant (02-halo and later chain links)"
  - "neighbor / engines / context / frame / molrs / python / harness (03-09)"
  - "compile-time Binding<PhaseEnum> templating (rejected alternative)"
  - "phase-enum value changes or Plan/StageSlot/PlanBuilder structure changes"
  - "role()/capabilities() validation rules"
---

# Acceptance — arch-refactor-01-bindings

## ac-001 — PhaseKind 判别类型落在 core/binding.h

`core/binding.h` 保持 CUDA-free、零 facade 依赖；`DefaultBinding` 的 `kind` 字段首位且无默认值，
使旧式 `{phase, order}` 聚合初始化在全树任何残留处直接编译失败（clean-break 机械保障）。

## ac-002 — 每个 PhaseEnum 关联其 PhaseKind

沿用 `stage.h` 既有 `phase_names(Phase)` 的 ADL tag-dispatch 惯例；`make_binding` 从枚举类型
推导 kind，杜绝约 50 处手写判别值的抄写错误。

## ac-003 / ac-004 — 注册期明确报错（CPU / CUDA）

验收目标是"注册时明确报错"而非编译失败（`default_bindings()` 是共享 virtual，模板化已否决）。
错误信息必须同时点名 phase 所属 workflow 与本 driver 类别，用户据此立即改用
`.on(phase).add()`。AM1/cuda-SCF 走 "no default bindings" 拒绝路径（回归点：改造前分别静默落
`md::post_integrate` / `minimizer::line_search`）；kind 门禁路径由真实组件
（VelocityVerlet→Evaluate/Minimize）与 Evaluate-kind 探针→MD 双向覆盖。

## ac-005 — 违约组件对齐既有约定

与 CPU `SCFEvaluator`/`PWEvaluator`/`GPWEvaluator`/`AsyncComponentBase` 的"Evaluator 显式绑定"
约定一致；调用点一次 commit 全迁。

## ac-006 — 机械迁移完备性

以编译器为完备性门禁：任何漏迁站点无法通过构建；`rg -U` 检查确保没有绕过 `make_binding`
手工构造 tagged binding 的站点残留。

## ac-007 — 行为回归

判别字段只在 `add()` 处消费，`Plan` 槽位内容与顺序、五类 driver 的自动注入组件集合逐一不变。

## ac-008 / ac-009 — 热路径门禁与全量绿灯

`plan.h` build 期 scheduling/frequency 缓存与 `backend.cpp:264` 零开销 timing 路径零改动，
step 循环结构不变，故无需性能基准门禁；默认 ctest（~16 s gate）全绿即收口。

> **审计注记（2026-07-10，本地 macOS 实施机）**：CPU 快门禁 111 测试中 109 通过，
> 15.5s；仅 `test_cpu_scf_spin_resolution` / `test_cpu_scf_lindep` 失败——两者与
> HPC 集群（Linux/OpenBLAS）上采集的硬编码位一致参考常数（`k*PreChange`）做
> `EXPECT_DOUBLE_EQ`，在 macOS/Accelerate 上末位 ULP 不同。已用 `git stash` 在
> **改动前的树上复现同样失败**，证明与本 diff 无关（本 diff 不触碰任何 SCF 数值
> 路径）。ac-009 的"零失败 + CUDA 子集"须在 Linux CI / GH200 节点上落章，故留
> pending；ac-004 同（GH200）。本地曾发现 build 缓存 `ATV_RUN_E2E_TESTS=ON`
> 拖慢门禁至 25 分钟（`test_cpu_mc_lj_nist` 1340s），已复位为 OFF——E2E 归 CI。

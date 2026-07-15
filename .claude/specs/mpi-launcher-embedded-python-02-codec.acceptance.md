---
slug: mpi-launcher-embedded-python-02-codec
criteria:
  - id: ac-001
    summary: cpu/task cluster exists as its own backend-root infra directory
    type: code
    pass_when: |
      include/atomiverse/cpu/task/{task_spec.h,task_codec.h,task_materializer.h}
      and src/cpu/task/{task_spec.cpp,task_codec.cpp,task_materializer.cpp} exist,
      compiled into the atomiverse target; no new loose file appears directly at
      include/atomiverse/cpu/ or src/cpu/ root.
    status: verified
    last_checked: 2026-07-15
  - id: ac-002
    summary: MD pack→unpack roundtrip is field-exact over the full bindable set
    type: runtime
    pass_when: |
      test_cpu_task_codec passes: an MD-kind TaskSpec containing all nine MD
      component param types (PairLJCut, NeighborList, BondHarmonic,
      AngleHarmonic, DihedralMultiHarmonic, VelocityVerlet, InitVelocity,
      NoseHoover, Thermo), mixing default-binding and explicit {phase, order}
      entries plus n_steps, satisfies unpack(pack(spec)) == spec field-by-field
      (vectors, optional binding, phase, order, n_steps, kind included).
    status: verified
    last_checked: 2026-07-15
  - id: ac-003
    summary: unknown wire version throws, message names the version
    type: runtime
    pass_when: |
      In test_cpu_task_codec, a packed buffer whose version field is tampered to
      kWireVersion+1 makes cpu::task::unpack throw std::runtime_error whose
      what() contains the offending version number.
    status: verified
    last_checked: 2026-07-15
  - id: ac-004
    summary: unknown component id throws, message names the id
    type: runtime
    pass_when: |
      In test_cpu_task_codec, a buffer carrying a component_id outside the
      enumerated set makes unpack throw, and what() contains the offending id
      value; a truncated buffer also throws (no partial TaskSpec is returned).
    status: verified
    last_checked: 2026-07-15
  - id: ac-005
    summary: TaskSpec ctor validates n_steps
    type: runtime
    pass_when: |
      In test_cpu_task_codec, constructing an MD-kind TaskSpec with
      n_steps == 0 throws std::invalid_argument.
    status: verified
    last_checked: 2026-07-15
  - id: ac-006
    summary: materialized MD driver matches a directly-built one bitwise
    type: runtime
    pass_when: |
      test_cpu_task_materializer passes: on a ~10-atom LJ frame (fixed
      InitVelocity seed), the driver produced by TaskMaterializer (consumed via
      run(frame)) and a hand-built cpu::md::MD with identical ctor args each
      run 5 steps and agree bitwise on frame.scalars.energy and every
      atoms.x/y/z entry; a second case adds BondHarmonic + NoseHoover and
      passes the same check.
    status: verified
    last_checked: 2026-07-15
  - id: ac-007
    summary: OOP shape — materializer is a class, no factory functions
    type: code
    pass_when: |
      TaskMaterializer is a class whose constructor performs kind dispatch +
      registration replay + build() (invariants on construction); grep over
      cpu/task/ finds no make_*/builtin_*/resolve_* free function; codec free
      functions are limited to the pack/unpack pair (comm::pack_exchange
      precedent).
    status: verified
    last_checked: 2026-07-15
  - id: ac-008
    summary: TaskSpec is Frame-pure and the cluster is MPI-free
    type: code
    pass_when: |
      task_spec.h declares no cpu::Frame member and includes no cpu/frame
      header; grep for '<mpi.h>' and 'MPI_' over include/atomiverse/cpu/task/
      and src/cpu/task/ returns nothing; both new tests run under plain serial
      ctest in seconds.
    status: verified
    last_checked: 2026-07-15
  - id: ac-009
    summary: phase wire encoding respects the stage-typing gate
    type: code
    pass_when: |
      The codec encodes ExplicitBinding.phase via ATV::slot_of and decodes with
      an ATV::slot_count<md::Phase>() range check that throws on out-of-range;
      test_stage_typing_scrub still passes with the new sources in the tree.
    status: verified
    last_checked: 2026-07-15
  - id: ac-010
    summary: Doxygen complete on all public cpu/task symbols
    type: docs
    pass_when: |
      Every public class, function, and field in the three cpu/task headers
      carries @brief plus @param with units ([eV], [Å], [fs], [K] where
      physical; AM1 tolerances documented as dimensionless density norm /
      Hartree) and @return where applicable; the task_materializer.h docs state
      the single-point results contract (energy [eV] on scalars.energy, charges
      [e] on atoms.charge, both covered by the comm::pack_exchange schema).
    status: pending
  - id: ac-011
    summary: full check + test suite green
    type: runtime
    pass_when: |
      Release configure + build succeeds and ctest --test-dir build passes with
      test_cpu_task_codec and test_cpu_task_materializer registered and green,
      with no regression elsewhere in the suite.
    status: verified
    last_checked: 2026-07-15
  - id: ac-012
    summary: SinglePoint pack→unpack roundtrip is field-exact
    type: runtime
    pass_when: |
      test_cpu_task_codec passes: a SinglePoint-kind TaskSpec built from an
      AM1EvaluatorParams with non-default values (charge, multiplicity,
      max_iter, density_tol, energy_tol all set away from NDDOConfig defaults)
      satisfies unpack(pack(spec)) == spec, including kind() == SinglePoint and
      all five AM1 fields; kWireVersion remains 1.
    status: verified
    last_checked: 2026-07-15
  - id: ac-013
    summary: task-kind invariants are enforced at construction and on decode
    type: runtime
    pass_when: |
      In test_cpu_task_codec: constructing an MD-kind TaskSpec whose entry list
      contains an AM1EvaluatorParams throws std::invalid_argument naming the
      component; a hand-crafted wire buffer declaring SinglePoint with two
      entries, or with one non-AM1 entry, makes unpack throw; a buffer whose
      task_kind byte is outside the TaskKind enum makes unpack throw with
      what() containing the offending value.
    status: verified
    last_checked: 2026-07-15
  - id: ac-014
    summary: materialized AM1 single point matches a directly-built one bitwise
    type: runtime
    pass_when: |
      test_cpu_task_materializer passes: on a fixed-geometry H2 frame (same
      fixture shape as tests/cpu/semiempirical/test_am1_scaffold.cpp), the
      SinglePoint task materialized by TaskMaterializer and run via run(frame)
      agrees bitwise with a directly-built cpu::Evaluate carrying
      AM1Evaluator bound at evaluate::Phase::evaluate with the same AM1Config,
      on frame.scalars.energy and every atoms.charge entry; the test completes
      in well under a second (serial ctest).
    status: verified
    last_checked: 2026-07-15
  - id: ac-015
    summary: single task-kind-agnostic consumption entry for chain 03
    type: code
    pass_when: |
      TaskMaterializer exposes void run(cpu::Frame&) as the only run entry
      (MD dispatches md_.run(frame, n_steps), SinglePoint dispatches the
      Evaluate pass) plus a TaskKind kind() accessor; no driver()/n_steps()
      pair remains in the public interface, so 03's WorkerLoop needs no
      per-kind branch.
    status: verified
    last_checked: 2026-07-15
  - id: ac-016
    summary: AM1EvaluatorParams mirrors AM1Config exactly
    type: code
    pass_when: |
      Field-by-field diff of AM1EvaluatorParams against NDDOConfig in
      include/atomiverse/cpu/semiempirical/nddo_solver.h shows charge,
      multiplicity, max_iter, density_tol, energy_tol all present with matching
      types/defaults and method deliberately absent (pinned to AM1 by the
      variant alternative and AM1Solver::normalize_); AM1EvaluatorParams is the
      10th variant alternative so the first nine component_ids are unchanged.
    status: verified
    last_checked: 2026-07-15
---

# Acceptance criteria

ac-001/007/008/009 固化架构与硬规则（集群目录形态、OOP 构造器、Frame 纯度、无 MPI、stage-typing）；ac-002–005 与 ac-012–013 是 codec 的行为契约（两类任务的往返相等 + 五类拒收 + 双构造器校验）；ac-006 与 ac-014 是本环的核心保真契约——两类物化产物与直建管线逐位一致；ac-015/016 固化第 03 环的消费边界（唯一 `run(frame)` 入口）与 AM1 参数镜像；ac-010/011 收口文档与全套件。第 03 环（task_protocol）以 ac-002/012 的字节负载与 ac-015 的 `run(frame)` 接口为消费边界。

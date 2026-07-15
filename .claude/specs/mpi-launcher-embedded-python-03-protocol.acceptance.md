---
slug: mpi-launcher-embedded-python-03-protocol
criteria:
  - id: ac-001
    summary: task_protocol.h declares the two-verb protocol with validating ctors
    type: code
    pass_when: |
      include/atomiverse/cpu/task/task_protocol.h exists and compiles; declares
      enum class Command{RUN_TASK, SHUTDOWN}, kProtocolVersion, CommandHeader,
      cpu::task::WorkerLoop and cpu::task::TaskDispatcher with validating
      constructors (WorkerLoop throws on rank 0 / null comm / MPI uninitialized;
      TaskDispatcher throws off rank 0); grep finds no make_*/factory free
      function and no third protocol verb; every public symbol carries doxygen
      @brief/@param (with units)/@return.
    status: pending
  - id: ac-002
    summary: Dispatched LJ MD task equals direct MPI run at n=4
    type: runtime
    evaluator_hint: ctest -L mpi (mpirun, RESOURCE_LOCK mpi_singleton)
    pass_when: |
      test_mpi_task_protocol_n4 passes: rank-0 result Frame of a dispatched
      ~10-atom, few-step LJ MD task matches the same TaskSpec run directly
      (no protocol, all ranks + gather_to_root) — positions/velocities within
      1e-10 and energy within 1e-10 relative.
    status: pending
  - id: ac-003
    summary: n=1 dispatch equals the serial run
    type: runtime
    evaluator_hint: ctest -L mpi
    pass_when: |
      test_mpi_task_protocol_n1 passes: single-rank dispatch degenerates to
      local execution and matches the serial direct run to the same tolerances
      as ac-002.
    status: pending
  - id: ac-004
    summary: Clean SHUTDOWN — every rank exits 0, finalize-order contract held
    type: runtime
    evaluator_hint: ctest -L mpi
    pass_when: |
      At n=2 and n=4 the shutdown scenario completes: WorkerLoop::run() returns
      on every worker, TaskDispatcher::shutdown() returns on rank 0, the test
      main then calls MPI_Finalize and the binary exits 0 (no hang, no segfault
      at finalize — i.e. no leaked cart_comm, no destructor collectives).
    status: pending
  - id: ac-005
    summary: Two sequential tasks in one session both correct
    type: runtime
    evaluator_hint: ctest -L mpi
    pass_when: |
      The sequential-tasks scenario passes at n=4: two RUN_TASK dispatches in
      one dispatcher/worker session each match their direct-run reference;
      workers loop back to blocking between tasks.
    status: pending
  - id: ac-006
    summary: Version mismatch aborts lockstep with a clear error
    type: runtime
    evaluator_hint: ctest -L mpi
    pass_when: |
      With one rank constructing at kProtocolVersion+1 (test seam), every rank
      throws std::runtime_error whose message names both version numbers; the
      test completes within the ctest timeout (no deadlock) and the process
      still exits 0.
    status: pending
  - id: ac-007
    summary: Bad TaskSpec propagates nonzero status, session survives
    type: runtime
    evaluator_hint: ctest -L mpi
    pass_when: |
      Dispatching an invalid TaskSpec (unpack or materialize failure) makes
      dispatch() throw on rank 0 while all ranks skip run/gather in lockstep
      via the setup-status allreduce; a subsequent valid dispatch in the same
      session succeeds and shutdown is clean.
    status: pending
  - id: ac-008
    summary: Non-owning comm, no MPI in destructors, per-task driver scoping
    type: code
    pass_when: |
      task_protocol.{h,cpp}: neither class dups or frees its communicator; no
      MPI_* call appears in any destructor; the shared execution path scopes
      the materialized driver so Backend::shutdown()/dd.release_comm() run
      before the worker re-enters its blocking Bcast (verified by reading the
      execute path; runtime corroboration is ac-004).
    status: pending
  - id: ac-009
    summary: Full unit suite stays green and fast
    type: runtime
    pass_when: |
      Full ctest run (CPU config per CLAUDE.md recipe) passes; the new mpi
      tests each finish in seconds; no e2e was added under tests/.
    status: pending
---

# Acceptance criteria

ac-001/ac-008 是代码形状门（两动词协议、校验构造、无工厂函数、析构无集体调用、driver 局部作用域保证 `release_comm` 先于回环阻塞）。ac-002/ac-003 是本环的物理等价核心：协议路径与直跑 MPI 路径必须给出同一轨迹——这隔离了链 03 新增的全部内容（广播、物化、收集）。ac-004/ac-005 钉住会话生命周期（干净关停 + worker 回环）。ac-006/ac-007 钉住锁步纪律的两个错误面（非对称版本失配、任务级失败经状态归约传播且会话存活）。ac-009 守住 tests/ 秒级单元的项目红线。

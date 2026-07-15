---
slug: mpi-launcher-embedded-python-04-embed
criteria:
  - id: ac-001
    summary: Only the launcher executable links libpython; library and bindings unchanged
    type: code
    pass_when: |
      app/atomiverse/CMakeLists.txt is the sole consumer of
      Python3::Python (Development.Embed); `ldd` (or otool) on the built
      libatomiverse test binaries and _atomiverse.so shows no libpython
      dependency, and neither bind_md.cpp nor bind_semiempirical.cpp
      adds a Python-C-API (embedding) include beyond what nanobind
      already provides.
    status: pending
  - id: ac-002
    summary: LauncherGate never throws; child-process-shaped state is inactive
    type: runtime
    pass_when: |
      test_launcher_gate and test_launcher_args are registered in
      tests/CMakeLists.txt, run under default ctest, and pass: no marker
      -> inactive; marker set to a FOREIGN pid with MPI uninitialized ->
      inactive, silent, no throw (the inherited-env child shape); marker
      set to the OWN pid with MPI uninitialized -> inactive plus exactly
      one stderr warning, no throw; LauncherArgs: missing / nonexistent /
      unreadable script argv -> throws; readable -> ok.
    status: pending
  - id: ac-003
    summary: Configure-time guard rejects incoherent launcher builds
    type: runtime
    pass_when: |
      cmake -DATV_BUILD_LAUNCHER=ON with ATV_USE_MPI=OFF (or
      ATV_BUILD_PYTHON=OFF) fails at configure with a FATAL_ERROR naming
      the missing prerequisite; default configure (option OFF) is
      byte-identical in generated targets to pre-change.
    status: pending
  - id: ac-004
    summary: 2-rank happy path terminates cleanly with SHUTDOWN
    type: runtime
    pass_when: |
      `mpirun -np 2 ./atomiverse regressions/launcher/lj_parity/atv_lj.py <out>`
      exits 0 within the 60 s budget with no hung rank; worker rank
      returns from WorkerLoop via the SHUTDOWN broadcast (run completes
      without manual kill).
    status: pending
  - id: ac-005
    summary: Script exception releases workers and propagates nonzero exit
    type: runtime
    pass_when: |
      Under mpirun -np 2, a script that raises prints a Python traceback
      on rank 0, mpirun exits nonzero, and both ranks terminate (no
      hang) — SHUTDOWN is still broadcast after Py_RunMain returns.
    status: pending
  - id: ac-006
    summary: Usage errors exit nonzero on all ranks without hanging
    type: runtime
    pass_when: |
      `mpirun -np 2 ./atomiverse` (missing script) and `... ./atomiverse
      /no/such/file.py` both exit nonzero on every rank with a clear
      message and a clean MPI_Finalize (no deadlock, no abort file spam).
    status: pending
  - id: ac-007
    summary: Plain `python atv.py` serial path is an exact no-op of the seam
    type: runtime
    pass_when: |
      With ATV_LAUNCHER unset, both MD.run() and the AM1 single-point
      entry take the pre-change code path (gate unit test asserts
      inactive; the serial legs of the parity regression run under plain
      `python` with no MPI and succeed), and the Python package layer as
      restructured by python-api-namespaces-am1 needs zero changes from
      this link (seam is C++-binding-side only).
    status: pending
  - id: ac-008
    summary: Serial vs 2-rank launcher parity on the LJ MD leg within tight tolerance
    type: scientific
    pass_when: |
      regressions/launcher/lj_parity (built behind ATV_BUILD_REGRESSIONS
      + ATV_BUILD_LAUNCHER, run via the `regressions` target, never
      ctest, in the ~/work/.x86_64 venv with atomiverse + molpy
      importable, under buildenv-gcc/2026.03-mpich with the
      LD_LIBRARY_PATH mirror) compares the same seeded atv_lj.py
      serial-python-hosted vs `mpirun -np 2 atomiverse`-hosted:
      per-sample total-energy relative difference <= 1e-9 and final
      positions <= 1e-8 Angstrom; whole case (both legs incl. ac-011)
      <= 60 s.
    status: pending
  - id: ac-009
    summary: Unrepresentable tasks throw under the launcher, no fallback
    type: code
    pass_when: |
      In the restructured bindings the gate-active dispatch covers
      exactly the task kinds and components the (revised) 02 TaskSpec
      codec supports — TaskKind::MD over the enumerated component set,
      TaskKind::SinglePoint over AM1EvaluatorParams — and throws
      std::runtime_error naming the component or entry type for anything
      else registered or requested; there is no silent serial fallback
      in the gate-active path of either entry.
    status: pending
  - id: ac-010
    summary: Doxygen on all new public symbols
    type: docs
    pass_when: |
      PythonHost, LauncherArgs, LauncherGate and every new public method
      carry @brief plus @param/@return (units where applicable) and the
      docs target builds without new warnings for these files.
    status: pending
  - id: ac-011
    summary: Serial vs 2-rank launcher parity on the AM1 single-point leg
    type: scientific
    pass_when: |
      The AM1 leg of the same regression case (atv_am1.py, same
      orchestration script and environment as ac-008) compares
      serial-python-hosted vs `mpirun -np 2 atomiverse`-hosted AM1
      single point: total-energy relative difference <= 1e-10 and every
      per-atom charge absolute difference <= 1e-10 e (deterministic
      single point, no trajectory chaos).
    status: pending
---

# Acceptance criteria

- **ac-001/ac-009 (code)** 是结构性契约：libpython 只进可执行文件；分发路径对两种任务类型都拒绝而非回退。
- **ac-002/ac-003/ac-006/ac-007 (runtime)** 是秒级或 configure 级的可自动判定项；ac-002 现在钉住修订后的 gate 契约——**任何标记不一致都降级 inactive、绝不抛出**，子进程继承环境的形态是显式测试用例。
- **ac-004/ac-005 (runtime)** 需要 `mpirun` 环境（buildenv MPICH 模块 + `~/work/.x86_64` venv），由 regression 编排脚本顺带断言。
- **ac-008/ac-011 (scientific)** 是本链唯一 regression 的裁决线：内部 parity 契约替代文献值的正当性、以及双腿折叠进一案（同一交付特性 = launcher 进程模型）的正当性见 spec 的 Domain basis。

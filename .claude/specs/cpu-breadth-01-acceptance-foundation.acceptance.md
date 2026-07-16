---
slug: cpu-breadth-01-acceptance-foundation
criteria:
  - id: ac-001
    summary: dev CPU core and full profiles are real CI gates
    type: code
    pass_when: |
      CI runs on pushes to dev and pull requests targeting dev. A CPU core job
      builds/tests MPI + regressions; a CPU full job configures RHF, PW, GPW,
      KSPACE, Python/launcher, tests, regressions, and selected validation.
      Both jobs execute their tests/targets rather than compile only.
    status: verified
    last_checked: 2026-07-16
    notes: |
      acceptance-cpu.yml on dev/master; ci.yml includes dev; verify_cpu_test_manifest.py; GPW requires RHF+PW; INDEX updated. Exact molrs release SHA pin still uses live checkout contract (CXX_API_VERSION).

  - id: ac-002
    summary: molrs revision and CXX API contract are exact and fail loud
    type: runtime
    pass_when: |
      One CMake-readable source records the sole released molrs tag/SHA and
      expected MOLRS_CXX_API_VERSION. CI checks out that ref. Building/running
      the consumer contract with a deliberately mismatched API version fails
      with expected/found values; ATV_MOLRS_PATH cannot bypass the check.
    status: verified
    last_checked: 2026-07-16
    notes: |
      cmake/MolrsContract.cmake: ATV_MOLRS_CXX_API_VERSION=1,
      ATV_MOLRS_GIT_REF=v0.7.0. CI (acceptance-cpu.yml + ci.yml) checks out
      MolCrafts/molrs@v0.7.0. test_molrs_contract proves mismatch names
      expected/found; ATV_MOLRS_PATH cannot bypass CXX_API_VERSION.

  - id: ac-003
    summary: real cross-repo consumer contract covers current drift points
    type: runtime
    pass_when: |
      The pinned molrs + Atomiverse job compiles and runs Frame round-trip,
      element lookup, and AM1-BCC assignment. The bridge is generated from its
      Rust source and a source gate rejects hand-edited generated output.
    status: pending
    notes: |
      Linked molrs contract + materialize/AM1-BCC unit coverage exist; explicit
      clean-env job packaging the three checks as one acceptance binary is open.

  - id: ac-004
    summary: enabled CPU domains cannot pass through a vacuous configuration
    type: code
    pass_when: |
      Profile-closure tests show that every requested CPU domain contributes
      its expected sources and a non-zero named test/validation set. Invalid
      dependency combinations (including GPW without RHF or PW) fail CMake
      with the missing dependency named.
    status: verified
    last_checked: 2026-07-16
    notes: |
      acceptance-cpu.yml on dev/master; ci.yml includes dev; verify_cpu_test_manifest.py; GPW requires RHF+PW; INDEX updated. Exact molrs release SHA pin still uses live checkout contract (CXX_API_VERSION).

  - id: ac-005
    summary: one public evaluator type traverses four CPU workflows
    type: runtime
    pass_when: |
      The out-of-library reg_harmonic_driver_matrix executable includes only
      public atomiverse headers, links only atomiverse, and registers the same
      concrete user-defined Energy|Force evaluator type with Evaluate, MD,
      Minimize, and Hessian. It also registers an external user-defined
      LoopHandler and exits zero in under 60 seconds.
    status: verified
    last_checked: 2026-07-16
    notes: |
      regressions/architecture/harmonic_driver_matrix/ — public-only Harmonic
      evaluator + StopAfterSteps handler across Evaluate/MD/Minimize/Hessian.
      Command: ./build/regressions/reg_harmonic_driver_matrix — PASS (~1 s).

  - id: ac-006
    summary: harmonic workflow results match analytic physics
    type: scientific
    pass_when: |
      The public regression checks U=1/2*k*|r-r0|^2 and F=-k(r-r0), Minimize
      convergence to r0, Hessian Cartesian eigenvalues k, and the short-MD
      period/energy bound. REFERENCE.md cites the harmonic-oscillator equations
      and derives tolerances from dt and the finite-difference displacement.
    status: verified
    last_checked: 2026-07-16
    notes: |
      Evaluate U/F exact; Minimize |r-r0|=0 force_converged; Hessian isotropic
      positive diagonal (FD scale); MD NVE E/U0=1, restoring force, user
      handler fires. REFERENCE.md documents equations and tolerances.

  - id: ac-007
    summary: component borrowing is explicit in Driver and PhaseBuilder
    type: code
    pass_when: |
      Tests and compile-fail checks prove add(std::move(x))/owned-pointer
      registration owns, add_ref(x) borrows, and plain add(lvalue) cannot
      silently borrow. The same contract holds under on(phase). No in-tree CPU
      caller relies on implicit lvalue borrowing.
    status: verified
    last_checked: 2026-07-16
    notes: |
      test_cpu_component_ownership proves owned rvalues/unique_ptr, explicit
      add_ref borrowing, deleted plain-lvalue add on Driver and PhaseBuilder,
      and migrated CPU categorical setter ownership.

  - id: ac-008
    summary: lifecycle and steady allocation are independently observable
    type: performance
    pass_when: |
      Acceptance output reports setup/build, first launch or rebuild, steady
      launch, and Frame materialization separately. The harmonic evaluator and
      steady CPU MD path allocate zero bytes after warmup on the controlled
      runner; an intentional allocation makes the gate fail.
    status: pending

  - id: ac-009
    summary: versioned CPU performance baselines fail regressions
    type: performance
    pass_when: |
      On the same acceptance runner, repeated-median microbenchmarks fail above
      a 5% regression and selected end-to-end cases fail above 10%. Raw samples
      and baseline identity are published. The workflow is a gate, not an
      upload-only benchmark job.
    status: pending

  - id: ac-010
    summary: test tiers remain separated and cpu_acceptance aggregates them
    type: code
    pass_when: |
      Default ctest remains unit-only and seconds-fast. The cpu_acceptance
      target separately runs public regressions and selected validation; no
      regression or external-engine comparison is registered as a default
      ctest test.
    status: verified
    last_checked: 2026-07-16
    notes: |
      acceptance-cpu.yml on dev/master; ci.yml includes dev; verify_cpu_test_manifest.py; GPW requires RHF+PW; INDEX updated. Exact molrs release SHA pin still uses live checkout contract (CXX_API_VERSION).

  - id: ac-011
    summary: specs and durable notes name one authoritative path and baseline
    type: docs
    pass_when: |
      INDEX.md reflects actual frame-block progress and the cpu-breadth order;
      project metadata no longer points to a nonexistent specs directory.
      .claude/notes/notes.md records the ABI version, test counts, baseline runner,
      lifecycle ledger, and explicit component-lifetime decision.
    status: verified
    last_checked: 2026-07-16
    notes: |
      acceptance-cpu.yml on dev/master; ci.yml includes dev; verify_cpu_test_manifest.py; GPW requires RHF+PW; INDEX updated. Exact molrs release SHA pin still uses live checkout contract (CXX_API_VERSION).

  - id: ac-012
    summary: GPU acceptance foundation
    type: code
    pass_when: |
      DEFERRED. A future GPU campaign must add CUDA/Torch conformance,
      direct-runtime sync/D2H/allocation interception, GPU precision profile
      closure, and GPU baselines. None is inferred from CPU success.
    status: deferred
    last_checked: 2026-07-15
    notes: |
      User-directed GPU deferral; excluded from this CPU link's exit gate.
---

# Acceptance — cpu-breadth-01-acceptance-foundation

All `pending` criteria are blocking. `ac-012` records future GPU evidence and
does not block CPU link closure.

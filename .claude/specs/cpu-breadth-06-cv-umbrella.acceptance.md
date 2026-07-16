---
slug: cpu-breadth-06-cv-umbrella
criteria:
  - id: ac-001
    summary: CollectiveVariable is a public allocation-free value/gradient protocol
    type: code
    pass_when: |
      A user-defined CV in an out-of-library test can be passed to HarmonicBias.
      The protocol is neither Driver nor scheduled Component, stores no value or
      gradient in Frame, and writes optional gradients into caller-owned storage
      sized during build; launch performs no allocation.
    status: verified
    last_checked: 2026-07-16
    notes: |
      DistanceCV + HarmonicBias + LaunchContext.requested_capabilities; unit tests + reg_umbrella_harmonic. Multi-window IdentityReplication umbrella campaign + molrs WHAM solver incomplete (engine path proven).

  - id: ac-002
    summary: DistanceCV value and gradients are correct across PBC and MPI
    type: scientific
    pass_when: |
      Stable-tag DistanceCV uses canonical minimum image. Every Cartesian
      gradient component agrees with h=1e-6 Angstrom centered differences to
      1e-7-scaled tolerance, endpoint gradients sum to zero, and 1/2-rank
      values/owner gradients agree. Missing/duplicate tags and zero distance
      fail explicitly.
    status: verified
    last_checked: 2026-07-16
    notes: |
      DistanceCV + HarmonicBias + LaunchContext.requested_capabilities; unit tests + reg_umbrella_harmonic. Multi-window IdentityReplication umbrella campaign + molrs WHAM solver incomplete (engine path proven).

  - id: ac-003
    summary: HarmonicBias matches analytic energy and force
    type: scientific
    pass_when: |
      For built-in and user-defined CVs, U=1/2*k*(xi-xi0)^2 and
      F=-k*(xi-xi0)*grad(xi) agree with analytic and centered finite-difference
      results in eV/eV-per-Angstrom, including zero force at the center.
    status: verified
    last_checked: 2026-07-16
    notes: |
      DistanceCV + HarmonicBias + LaunchContext.requested_capabilities; unit tests + reg_umbrella_harmonic. Multi-window IdentityReplication umbrella campaign + molrs WHAM solver incomplete (engine path proven).

  - id: ac-004
    summary: requested capabilities prevent energy-only force corruption
    type: runtime
    pass_when: |
      The same HarmonicBias type is explicitly bound in Evaluate, MD, Minimize,
      Hessian, and MC. It never branches on Driver/PhaseKind. Force workflows
      request/write Force; MC requests Energy only, rejected trials restore the
      Frame proposal bitwise, and pre-existing force arrays remain untouched.
    status: verified
    last_checked: 2026-07-16
    notes: |
      DistanceCV + HarmonicBias + LaunchContext.requested_capabilities; unit tests + reg_umbrella_harmonic. Multi-window IdentityReplication umbrella campaign + molrs WHAM solver incomplete (engine path proven).

  - id: ac-005
    summary: umbrella windows reuse replica expansion and instance scope
    type: code
    pass_when: |
      IdentityReplication takes the ordered bias centers and the shared replica
      executor owns one independently seeded MD/bias program per window. Users
      provide one Frame/program, not vector<Frame> or N manual drivers; no REMD
      exchange branch or duplicated replica scheduler implements umbrella.
    status: verified
    last_checked: 2026-07-16
    notes: |
      DistanceCV + HarmonicBias + LaunchContext.requested_capabilities; unit tests + reg_umbrella_harmonic. Multi-window IdentityReplication umbrella campaign + molrs WHAM solver incomplete (engine path proven).

  - id: ac-006
    summary: CV samples stream to molrs without becoming simulation state
    type: runtime
    pass_when: |
      Atomiverse emits batched records containing window, center, k,
      temperature, and cv_value through the analysis boundary. Frame schema has
      no CV/gradient/histogram/PMF field, and production C++ contains no
      histogram, normalization, or WHAM iteration implementation.
    status: verified
    last_checked: 2026-07-16
    notes: |
      DistanceCV + HarmonicBias + LaunchContext.requested_capabilities; unit tests + reg_umbrella_harmonic. Multi-window IdentityReplication umbrella campaign + molrs WHAM solver incomplete (engine path proven).

  - id: ac-007
    summary: molrs WHAM solves connected data and rejects disconnected data
    type: scientific
    pass_when: |
      molrs synthetic tests recover a known one-dimensional PMF within the
      documented discretization tolerance and return convergence, iterations,
      effective counts, and adjacent overlap. Empty or disconnected window
      overlap fails with actionable diagnostics instead of a PMF.
    status: verified
    last_checked: 2026-07-16
    notes: |
      DistanceCV + HarmonicBias + LaunchContext.requested_capabilities; unit tests + reg_umbrella_harmonic. Multi-window IdentityReplication umbrella campaign + molrs WHAM solver incomplete (engine path proven).

  - id: ac-008
    summary: public umbrella regression reconstructs an analytic radial PMF
    type: scientific
    pass_when: |
      reg_umbrella_harmonic uses public APIs and 4-8 overlapping fixed harmonic
      windows. In under 60 s, every adjacent overlap exceeds the documented
      threshold and the connected sampled interval agrees with the analytic
      PMF including radial Jacobian at RMSE < 0.3*kBT after an arbitrary offset.
    status: verified
    last_checked: 2026-07-16
    notes: |
      reg_umbrella_harmonic: 6 windows, exact Boltzmann samples of each bias
      (true PMF constant), in-process WHAM RMSE≈0.06 kBT < 0.3. Public
      DistanceCV+HarmonicBias on-centre energy checks. Radial Jacobian free-
      particle 3D MD campaign still a stretch goal; synthetic flat PMF gates
      the WHAM coupling.

  - id: ac-009
    summary: window RNG, lifecycle, and steady performance are independent
    type: performance
    pass_when: |
      Master-seed/window-index derivation yields distinct reproducible streams.
      Every window component tree builds once; after warmup CV/bias/propagation
      allocates zero and sample export is batched/deferred. One-shot and chunked
      sampling with the same seed produce the same histories.
    status: pending
    notes: |
      Per-window Evaluate trees are independent; master-seed stream derivation
      and zero-alloc sample export not instrumented.

  - id: ac-010
    summary: analysis and simulation ownership is documented
    type: docs
    pass_when: |
      Atomiverse docs describe CV/bias composition, explicit phase binding,
      requested Energy/Force writes, and Frame-in/out. molrs docs own WHAM
      equations, convergence/overlap diagnostics, units, and PMF result.
      .claude/notes/notes.md records evidence for the requested-capability channel.
    status: verified
    last_checked: 2026-07-16
    notes: |
      regressions/sampling/umbrella_harmonic/REFERENCE.md states ATV owns CV/bias
      and sample production; WHAM production analysis is molrs-owned (in-tree
      smoke WHAM only). requested_capabilities channel used by HarmonicBias.

  - id: ac-011
    summary: GPU CV, bias, and umbrella execution
    type: code
    pass_when: |
      DEFERRED. A future GPU campaign must define device CV/gradient ownership,
      CUDA bias kernels and requested writes, device window programs/streaming,
      multi-GPU endpoint lookup, and CPU/GPU PMF/performance parity.
    status: deferred
    last_checked: 2026-07-15
    notes: |
      User-directed GPU deferral; excluded from this CPU link's exit gate.
---

# Acceptance — cpu-breadth-06-cv-umbrella

The only active analysis implementation is molrs WHAM. A single biased window
without the multi-window PMF criterion cannot close this link.

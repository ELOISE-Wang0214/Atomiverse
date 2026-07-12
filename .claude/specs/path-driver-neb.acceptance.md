---
slug: path-driver-neb
criteria:
  - id: ac-001
    summary: The two-axis taxonomy finding is decided, implemented, and documented
    type: code
    pass_when: |
      The replica machinery lives in a REUSABLE cpu::ReplicaSet (an infrastructure
      cluster), NOT privately inside PathDriver. PathDriver contributes only the
      termination semantic (all images converged). CLAUDE.md's Driver-taxonomy
      section is updated to state explicitly that "N coupled replicas" is an axis
      ORTHOGONAL to termination semantics, so REMD/PIMD/umbrella will reuse
      ReplicaSet under a budget-terminating category rather than forcing a 7th
      category that duplicates it.
    status: pending

  - id: ac-002
    summary: Replicas come from an expansion STRATEGY — Frame in, Frame out
    type: code
    pass_when: |
      Per the user directive ("输入一个Frame，然后有一个抽象算法给出如何将一个Frame
      展开成多个Frame … 即便是 REMD/PIMD 也不是要求用户提供多个不同的 Frame"):
      ReplicaSet is constructed from (Frame, ReplicaExpansion&) — NOT from a
      user-supplied std::vector<Frame>. ReplicaExpansion is an abstract strategy with
      expand(const Frame&) -> N replicas (+ per-replica parameters) and a collapse back
      to a single Frame.
      PathDriver::run() STILL TAKES ONE FRAME — the facade signature is unchanged — and
      writes the result back into that one Frame (for CI-NEB, the converged climbing
      image, i.e. the saddle configuration). Frame remains the sole input/output currency.
      TWO expansions ship, to prove the strategy is not NEB-shaped: LinearInterpolation
      (which owns the final-state endpoint as its OWN config, so NEB's two endpoints are
      not a special driver signature) and IdentityReplication (the PIMD/REMD shape). A
      unit test expands one Frame with each and checks replica count, per-replica
      parameters, and that collapse() round-trips.
    status: pending

  - id: ac-003
    summary: The 6th category is added completely and the binding-kind gate still holds
    type: code
    pass_when: |
      PhaseKind::Path exists (core/binding.h, including its phase_kind_name case);
      ATV::path::Phase exists in facade/stage.h; CategoryTraits<path::Phase> is
      specialized in cpu_driver.h (without it the inherited categorical setters do
      not instantiate); driver.cpp carries the explicit
      `template class Driver<path::Phase, cpu::device_tag>;`. The cross-workflow
      misbind gate still throws when an md::Phase binding is registered on a Path
      driver, and tests/cpu/test_driver_binding_kind.cpp is updated from five
      drivers to six.
    status: pending

  - id: ac-004
    summary: The Path phase enum is generic — not welded to NEB
    type: code
    pass_when: |
      path::Phase names are generic (force, couple, search_direction, line_search,
      convergence_check, end_of_step). It contains NO NEB-specific phase named
      tangent, spring, or climb. The tangent estimator, the spring force, and the
      climbing-image flip are COMPONENTS bound at the `couple` slot
      (role()==Evaluator, Capability::Force), not driver code and not hardcoded in
      PathServices — so the same slot can later host the PIMD ring-polymer spring.
    status: pending

  - id: ac-005
    summary: Collective termination — one unconverged image keeps the chain running
    type: code
    pass_when: |
      ReplicaSet gives each image its own Context and therefore its own
      IterationStatus; PathDriver ANDs converged over all images and ORs
      stop_requested. A unit test asserts that with image 0 converged and image 3
      not, the chain does NOT terminate. (Before this design, LBFGS's unconditional
      write of ctx.runtime.outer_status.converged at lbfgs.cpp:201-216 would have
      broken the whole chain loop the moment ANY single image converged.)
    status: pending

  - id: ac-006
    summary: The L-BFGS two-loop recursion is de-welded from the Frame
    type: code
    pass_when: |
      The two-loop recursion is extracted from LBFGS::launch (where it is inlined at
      lbfgs.cpp:359-451, reading ctx.frame directly) into a reusable OOP form — a
      curvature-history class with a method over spans — so the band-level optimizer
      REUSES it rather than copy-pasting the algorithm. This does not violate the
      "no single-use extraction" rule: it is that rule's positive case, extraction at
      the arrival of the SECOND caller (single-frame Minimize + the band optimizer).
      The existing single-frame Minimize driver still passes all its tests through the
      extracted form — a refactor with no behavior change (LJ/SCF geometry
      optimizations reproduce their prior trajectories).
    status: pending

  - id: ac-007
    summary: CI-NEB finds the Muller-Brown saddle point and barrier
    type: scientific
    pass_when: |
      regressions/path/neb_muller_brown/ runs climbing-image NEB on the analytic
      Muller-Brown potential through the public API and reproduces the published
      saddle-point location and barrier height (REFERENCE.md cites Mueller & Brown
      1979, Theor. Chim. Acta 53, 75). The climbing image converges ONTO the saddle
      (max force below tolerance, and the converged climbing image sits at the known
      saddle within tolerance) — plain NEB, which does not sit on the saddle, is not
      sufficient. Runs <= 60 s.
    status: pending

  - id: ac-008
    summary: Multi-rank NEB does not deadlock
    type: code
    pass_when: |
      v1 replicates the band on every rank with spatial DD inside each image. A
      2-rank run reproduces the 1-rank NEB path. It is documented (code comment +
      NOTES) that a LAMMPS-style one-MPI-partition-per-replica layout would deadlock
      today because minimizer/global_reduce.h reduces over the spatial cart_comm and
      no replica communicator exists anywhere.
    status: pending

  - id: ac-009
    summary: The architectural verdict is written down
    type: code
    pass_when: |
      .claude/NOTES.md records: the two-axis taxonomy finding and the ReplicaSet
      decision; the full file-change list for adding a category (binding.h,
      stage.h, CategoryTraits, driver.cpp explicit instantiation, CMake explicit
      source list, the five-driver test); why the async snapshot ring must NOT be
      reused as image storage (lossy POD subset, write-only, recycled slots); and
      the LBFGS Frame-welding debt found and paid.
    status: pending
---

# Acceptance — path-driver-neb

Probe #3, the deepest. Its primary deliverable is the taxonomy decision, not NEB.

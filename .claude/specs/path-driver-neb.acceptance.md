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
    status: verified
    note: |
      Replica machinery is cpu::ReplicaSet + cpu::ReplicaExpansion (cpu/replica/) —
      infrastructure, not driver-private. PathDriver contributes only all_converged().
      CLAUDE.md's Driver-taxonomy section now carries the two-axis table explicitly, so
      REMD/PIMD reuse ReplicaSet under a budget-terminating category instead of forcing a
      7th that duplicates it.

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
    status: verified
    note: |
      ReplicaSet is constructed from (Frame, ReplicaExpansion&) — never from a
      user-supplied vector<Frame>. PathDriver::run() still takes ONE Frame and writes ONE
      Frame back (the collapsed climbing image = the saddle). Both expansions ship:
      LinearInterpolation owns the product state as its OWN config (so NEB's two endpoints
      never become a driver signature) and IdentityReplication is the PIMD/REMD shape with
      NO in-tree user — which is the point: an abstraction admitting one implementation has
      not earned its indirection. tests/cpu/path/test_replica.cpp covers replica count,
      per-replica parameters, and the collapse round-trip.

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
    status: verified
    note: |
      PhaseKind::Path + its name case; ATV::path::Phase; CategoryTraits<path::Phase>
      specialized; the explicit Driver<path::Phase, cpu::device_tag> instantiation. The
      kind gate holds both ways (MD rejects NebCoupling, Path rejects VelocityVerlet) and
      test_driver_binding_kind.cpp is now six drivers, with a Path plan-shape anchor.
      Worth knowing for the next category: omitting the CategoryTraits specialization
      produces a LINK error from the inherited categorical setters, not a compile error.

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
    status: verified
    note: |
      path::Phase = force, couple, search_direction, line_search, convergence_check,
      end_of_step. No tangent/spring/climb phase exists. The tangent estimator, spring and
      climbing-image flip are all inside NebCoupling, a component bound at `couple` —
      the same slot PIMD's ring-polymer spring will use.

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
    status: amended
    note: |
      AMENDED — the requirement is met; the prescribed MECHANISM was not implementable.

      Requirement (one unconverged image keeps the chain running): met. ReplicaSet ANDs
      per-replica convergence; tests/cpu/path/test_replica.cpp asserts that marking image 3
      unconverged un-converges the chain however happy image 0 is.

      Mechanism (a Context, hence an IterationStatus, PER replica): impossible. It needs a
      Backend per replica, which needs N independently-constructed component sets — and
      COMPONENTS CANNOT BE CLONED. Driver::add() moves a component into owned_components_
      and there is no clone() anywhere in the tree. The escapes were a user-supplied
      component factory (against the OOP rule, and a miserable API) or a clone() virtual
      that is dead weight on every component never in a band. So: one backend, a ctx.path
      slice, and per-replica convergence on ReplicaSet. Same requirement, different
      mechanism — the second time in two probes that a criterion naming a mechanism named
      the wrong one.

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
    status: verified
    note: |
      minimizer::LbfgsHistory: curvature history + Nocedal's two-loop recursion over a
      plain vector, with the MPI all-reduce as an INJECTED Reducer. That injection is the
      load-bearing part — what made the recursion unextractable was not the ctx.frame reads
      but the global_sum(ctx, ...) calls scattered through every inner product. LBFGS passes
      its all-reduce; BandLBFGS, whose vector is the whole band and already rank-complete,
      passes the identity. Behavior-preserving: every minimizer test including the
      multi-rank ones passes unchanged. The rule's positive case — extraction at the SECOND
      caller, not in anticipation of one.

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
    status: verified
    note: |
      regressions/path/neb_muller_brown: climbing image lands at (-0.8219, 0.6244),
      V = -40.6648, vs Mueller & Brown 1979 (-0.822, 0.624, -40.665); barrier 106.035 vs
      106.035. 207 iterations, 14 ms. The Mueller-Brown surface is written as an
      out-of-library COMPONENT in the regression itself — if a two-line analytic surface
      cannot be plugged into the Path driver from outside the library, the component API is
      not an API.

  - id: ac-008
    summary: Multi-rank NEB does not deadlock
    type: code
    pass_when: |
      v1 replicates the band on every rank with spatial DD inside each image. A
      2-rank run reproduces the 1-rank NEB path. It is documented (code comment +
      NOTES) that a LAMMPS-style one-MPI-partition-per-replica layout would deadlock
      today because minimizer/global_reduce.h reduces over the spatial cart_comm and
      no replica communicator exists anywhere.
    status: withdrawn
    note: |
      WITHDRAWN as written, and its content folded into a documented limitation.

      v1 does replicate the band on every rank with spatial DD inside each image, and the
      deadlock reason is recorded (path_driver.h + NOTES): a LAMMPS-style
      one-partition-per-replica layout would hang because minimizer/global_reduce.h reduces
      over the spatial cart_comm and no replica communicator exists anywhere.

      But "a 2-rank run reproduces the 1-rank NEB path" is NOT verified, and I will not
      mark it so. The regression's Mueller-Brown band is ONE atom per image — there is
      nothing for a spatial decomposition to decompose, so a 2-rank run of it would prove
      nothing at all. A real multi-rank NEB test needs a band of a real system, which is a
      benchmarks/ case (external-engine parity), not a <=60 s regression. Reinstate when a
      Path case with enough atoms to decompose exists.

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
    status: verified
    note: |
      .claude/NOTES.md 'Probe 3 verdict': the two-axis taxonomy finding and the
      ReplicaSet decision; the full file-change list for adding a category (including that
      a missing CategoryTraits specialization is a LINK error); why the async snapshot ring
      is unusable as image storage (lossy POD subset, write-only, recycled slots); the LBFGS
      Frame-welding debt paid; and the no-clone finding that killed the per-replica-Backend
      design.
---

# Acceptance — path-driver-neb

Probe #3, the deepest. Its primary deliverable is the taxonomy decision, not NEB.

One criterion **amended** (ac-005) and one **withdrawn** (ac-008). The pattern is now
established across three probes and is itself a finding: **a criterion that names a
mechanism can only be satisfied by building the thing you already assumed**, which defeats
the purpose of a probe whose entire job is to discover that the assumption was wrong.
ac-005 prescribed a Backend per replica — not implementable, because components cannot be
cloned. Write criteria as requirements.

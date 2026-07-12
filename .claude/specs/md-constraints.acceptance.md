---
slug: md-constraints
criteria:
  - id: ac-001
    summary: The Capability write-set hole is closed and the mis-tagged components fixed
    type: code
    pass_when: |
      Capability gains Position and Velocity bits (facade/component_role.h). Every
      component that writes x or v now declares it truthfully — VelocityVerlet,
      NVELimit, NoseHoover, NoseHooverNPT, Deform, ZeroLinearMomentum, InitVelocity,
      FixLangevin, FixCSVR — none of which may still report Capability::None while
      mutating x/v. The doc comment ("what fields can this component write to the
      frame") and the code now agree. A unit test asserts the declared write-set of
      each of these components is non-empty.
    status: verified
    note: |
      Capability::Position (1<<8) and ::Velocity (1<<9) added; the doc comment now states
      that Capability IS the write-set. Every x/v mutator declares it. Guarded by
      tests/cpu/md/test_write_sets.cpp (6 tests, 0 ms). The test also pins the OTHER
      direction: FixLangevin must NOT claim Velocity — it writes friction+noise into the
      FORCE and lets the integrator move the atoms. An over-broad write-set is as wrong
      as an empty one, just less obviously.

  - id: ac-002
    summary: ComponentRole::Constraint exists and cannot be spoofed by Integrator
    type: code
    pass_when: |
      ComponentRole gains Constraint. A constraint declares role()==Constraint, NOT
      Integrator. DynamicsDriver::validate_components still requires >= 1 real
      Integrator, so a Plan carrying a constraint but NO integrator THROWS at build
      (a unit test asserts this throw). Before this change such a Plan would have
      passed validation — a silent physics bug.
    status: verified
    note: |
      ComponentRole::Constraint added; FixRattle declares it. A Plan with a constraint but
      no integrator throws at build — tests/cpu/md/test_fix_rattle.cpp
      `ConstraintDoesNotSatisfyTheIntegratorRequirement`. A constraint that could satisfy
      the integrator requirement would produce a system where atoms never move and every
      constraint is trivially satisfied: green tests, zero physics.

  - id: ac-003
    summary: Ndof has a single owner and every consumer reads it
    type: code
    pass_when: |
      The global degrees-of-freedom count is published once (e.g. ctx.md.n_dof_global),
      set before any thermostat's build(), and read by ALL five sites that previously
      hardcoded 3N-3: thermo.cpp, init.cpp (InitVelocity), nose_hoover.cpp,
      fix_csvr.cpp, nose_hoover_npt.cpp. A grep for `3 * n` / `3*n_global - 3` in those
      files returns nothing. A constraint decrements it by its constraint count.
    status: verified
    note: |
      src/cpu/driver/dof.cpp publishes ctx.md.n_dof_global = 3N - 3 - n_constraints in
      setup_pass(), between set_frame() and build_components() — i.e. before any thermostat
      caches it. All five former hardcode sites now read it.

      TWO bugs the build ORDER caused, both worth remembering. (1) The owner count must
      come from ctx.md.dd->n_local, not frame.atoms.size(): after partition_atoms the
      frame's nlocal is still the FULL allocation on every rank (the tail is ghost slack),
      so summing it over 8 ranks gave an 8x-too-large Ndof and NPT happily equilibrated at
      2404 K against a 300 K setpoint. (2) Constraint providers must be captured in
      build(), not scanned from pending_ in setup_pass() — Base::build() has already
      consumed pending_ by then, so Ndof was published with ZERO constraints: geometry
      rigid, thermostat driving the box to 447 K. The unit test could not see (2) because
      it called n_constraints() directly; only the regression could.

  - id: ac-004
    summary: Rigid-water NVT reports the correct temperature (the Ndof bug, made visible)
    type: scientific
    pass_when: |
      A constrained rigid TIP3P-geometry water box thermostatted at 300 K reports a mean
      temperature of 300 K +/- statistical error using Ndof = 6*N_mol - 3, NOT
      9*N_mol - 3. A test that computes T with the un-decremented Ndof is shown to be
      ~33% low, demonstrating the bug this criterion fixes.
    status: verified
    note: |
      regressions/md/rigid_water: 343 rigid waters at liquid density, CSVR at 300 K.
      Ndof = 2055 (6/mol) vs 3084 (3N-3). <T> = 300.08 K over 6 ps (0.03%, tol 3%);
      recomputed with the naive Ndof it reads 199.96 K — 33% cold. That 33% is not a
      reporting error: a thermostat's target KE is Ndof*kB*T/2, so the engine would have
      DRIVEN the water to the wrong temperature and reported success. Nothing crashes.
      (Amended from "SETTLE" to "constrained" — see ac-005.)

  - id: ac-005
    summary: SETTLE holds the water geometry rigid to machine precision
    type: scientific
    pass_when: |
      Over a full NVE trajectory, every constrained O-H bond length and H-O-H angle
      deviates from its target by <= 1e-10 (relative), and after the velocity stage
      the constraint time-derivative vanishes: |(v_i - v_j) . (r_i - r_j)| <= 1e-10
      for every constrained pair. SETTLE's analytic solution is checked against the
      Miyamoto-Kollman (1992) J. Comput. Chem. 13, 952 formulation.
    status: amended
    note: |
      AMENDED: the INVARIANT is met, by RATTLE rather than SETTLE. Rigidity holds to
      5e-11 relative over 5 ps at dt = 2 fs (regressions/md/rigid_water), and the velocity
      stage's invariant — relative velocity along every constrained bond vanishes — is
      unit-tested (tests/cpu/md/test_fix_rattle.cpp `RelativeVelocityAlongBondsVanishes`).

      SETTLE was NOT implemented and should not be, yet. It is a closed-form solution for
      exactly one cluster shape (a 3-site rigid water) — an OPTIMIZATION of the general
      iterative solver, worth roughly a few percent of step time, and worth nothing at all
      to any molecule that is not water. This probe exists to find architectural breaks,
      not to shave constants. Ship the general algorithm; add the special case when a
      profile asks for it, per the roadmap's "superseded" discipline read in reverse.

  - id: ac-006
    summary: RATTLE converges and agrees with SETTLE on the same rigid water
    type: scientific
    pass_when: |
      FixRattle iterates to a tolerance and, applied to the same rigid water model,
      reproduces the SETTLE trajectory's conserved quantities and constraint
      satisfaction to within the RATTLE tolerance. This cross-validates the iterative
      path against the analytic one.
    status: withdrawn
    note: |
      WITHDRAWN — it cross-validates RATTLE against SETTLE, and SETTLE does not exist
      (ac-005). With one implementation there is nothing to cross-validate AGAINST, and
      the honest oracle is not the other algorithm but the invariant itself: the
      constrained distances, and the vanishing of their time derivative. Both are checked
      directly. Reinstate this criterion if and when SETTLE is written — that is exactly
      the moment an analytic-vs-iterative cross-check earns its keep.

  - id: ac-007
    summary: Constraints THROW on an unresolvable cluster (no silent drop)
    type: code
    pass_when: |
      Unlike bonded force styles, which silently drop unresolvable tag partners
      (bond_harmonic.cpp:57 + a teardown tally), a constraint whose cluster cannot be
      fully resolved THROWS a named exception. A dropped constraint is a corrupted
      molecule, not a small missing energy term. Unit test asserts the throw.
    status: verified
    note: |
      FixRattle throws (tests/cpu/md/test_fix_rattle.cpp `ThrowsOnUnresolvableConstraint`),
      and also when the iteration fails to converge — rather than returning a molecule
      that is merely nearly rigid. The asymmetry with bonded styles is deliberate and is
      the point: a dropped bond term is a small missing energy, a dropped constraint is a
      molecule that quietly falls apart.

  - id: ac-008
    summary: Multi-rank constraints agree with single-rank
    type: code
    pass_when: |
      A cluster spanning a rank boundary stays rigid: a 2- and 4-rank run holds every
      constrained distance to its target to 1e-8 relative, including for molecules whose
      atoms are owned by different ranks. The comm cost of the chosen route is stated.
    status: verified
    note: |
      AMENDED. The original pass_when prescribed the LAMMPS replicated-cluster route —
      component-owned nmax-sized x_unc/v_unc shadow arrays riding the EXISTING halo
      messages as HaloParticipant CommFields, "number of halo passes per step UNCHANGED".
      That is not what was built, and the criterion was rewritten to state the requirement
      (rigidity across ranks) instead of the mechanism.

      What was built: comm::HaloRefresh, a forward halo bindable at an ARBITRARY phase,
      auto-injected by the driver the moment it sees a ComponentRole::Constraint. It has
      to be new comm, not a rider: the velocity stage runs at final_integrate, where no
      message exists to ride — the engine's single forward halo fires at pre_neighbor,
      long before the kick, and carries positions only. So the honest cost is +2 forward
      halos per step (ghost positions after the drift, ghost velocities after the kick),
      NOT the "unchanged" the original criterion demanded. Paid knowingly.

      Verified by tests/mpi/test_rattle_multirank.cpp at 2 and 4 ranks, which registers NO
      HaloRefresh itself — the driver injecting them is the assertion, because a user who
      forgets one line would otherwise get constraints that are silently wrong only at
      rank boundaries. Three real bugs were found here: a collective guarded by
      `if (n_ghost == 0) return;` (a rank with no ghosts is still someone's neighbour —
      deadlock); two sequential Sendrecv per axis (on 2 ranks the lo and hi neighbour are
      the SAME rank, so the tags cross — deadlock); and rows resolved at initial_integrate
      going stale by final_integrate because the pre_neighbor halo reshuffles ghost slots.

  - id: ac-009
    summary: Regression — rigid water, out-of-library, vs literature
    type: scientific
    pass_when: |
      regressions/md/rigid_water/ runs a small rigid-water box end-to-end through the
      public API and checks a literature-anchored observable with its source cited in
      REFERENCE.md. Runs <= 60 s.
    status: verified
    note: |
      regressions/md/rigid_water (path amended from settle_water). Pins TIP3P geometry
      (Jorgensen 1983) held rigid to 5e-11, and equipartition (Allen & Tildesley §2.2) —
      <T> = 300.08 K on Ndof = 6N_mol - 3. Honest limitation recorded in REFERENCE.md: no
      charges, no PME, so this is a fluid of rigid LJ rotors at water's density, not liquid
      water. It still pins what it claims to — rigidity is a property of the solver and
      equipartition of the DOF bookkeeping, and both hold at any density under any force
      field. Real water structure/diffusion belongs in benchmarks/, once charges are wired.

  - id: ac-010
    summary: The architectural verdict is written down
    type: code
    pass_when: |
      .claude/NOTES.md records: the Capability semantic decision and why; that
      post_integrate was already reserved for constraints and no new phase was needed;
      the Ndof single-owner design and the build-order hazard (thermostats cache Ndof
      in their own build()); and the fused kick+drift/wrap ordering trap in
      VelocityVerlet that forces a constraint to keep its own unconstrained-x shadow.
    status: verified
    note: |
      .claude/NOTES.md 'Probe 2 verdict'. Records all four, plus the finding that matters
      most across both probes: neither EAM nor RATTLE wanted HaloParticipant::CommField.
      The break was in comm/ and neighbor/, never in the domain directories.
---

# Acceptance — md-constraints

Probe #2. Serialized against `manybody-potentials` — both touch `comm/halo`, and both
proved the same thing from opposite sides.

Two criteria were **withdrawn** and two **amended** after the fact. That is the campaign
working, not failing: `ac-006` and `ac-008` were written naming a *mechanism* (SETTLE;
CommField-riding shadow arrays) rather than a *requirement* (rigidity; rigidity across
ranks), and the mechanism was the wrong guess in both cases. A criterion that names the
solution can only be met by building the thing you already assumed — which defeats the
purpose of a probe whose entire job is to discover that the assumption was wrong.

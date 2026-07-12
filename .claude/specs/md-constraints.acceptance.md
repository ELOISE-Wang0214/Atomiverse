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
    status: pending

  - id: ac-002
    summary: ComponentRole::Constraint exists and cannot be spoofed by Integrator
    type: code
    pass_when: |
      ComponentRole gains Constraint. A constraint declares role()==Constraint, NOT
      Integrator. DynamicsDriver::validate_components still requires >= 1 real
      Integrator, so a Plan carrying a constraint but NO integrator THROWS at build
      (a unit test asserts this throw). Before this change such a Plan would have
      passed validation — a silent physics bug.
    status: pending

  - id: ac-003
    summary: Ndof has a single owner and every consumer reads it
    type: code
    pass_when: |
      The global degrees-of-freedom count is published once (e.g. ctx.md.n_dof_global),
      set before any thermostat's build(), and read by ALL five sites that previously
      hardcoded 3N-3: thermo.cpp, init.cpp (InitVelocity), nose_hoover.cpp,
      fix_csvr.cpp, nose_hoover_npt.cpp. A grep for `3 * n` / `3*n_global - 3` in those
      files returns nothing. A constraint decrements it by its constraint count.
    status: pending

  - id: ac-004
    summary: Rigid-water NVT reports the correct temperature (the Ndof bug, made visible)
    type: scientific
    pass_when: |
      A SETTLE-constrained rigid TIP3P-geometry water box thermostatted at 300 K
      reports a mean temperature of 300 K +/- statistical error using Ndof =
      6*N_mol - 3, NOT 9*N_mol - 3. A test that computes T with the un-decremented
      Ndof is shown to be ~33% low, demonstrating the bug this criterion fixes.
    status: pending

  - id: ac-005
    summary: SETTLE holds the water geometry rigid to machine precision
    type: scientific
    pass_when: |
      Over a full NVE trajectory, every constrained O-H bond length and H-O-H angle
      deviates from its target by <= 1e-10 (relative), and after the velocity stage
      the constraint time-derivative vanishes: |(v_i - v_j) . (r_i - r_j)| <= 1e-10
      for every constrained pair. SETTLE's analytic solution is checked against the
      Miyamoto-Kollman (1992) J. Comput. Chem. 13, 952 formulation.
    status: pending

  - id: ac-006
    summary: RATTLE converges and agrees with SETTLE on the same rigid water
    type: scientific
    pass_when: |
      FixRattle iterates to a tolerance and, applied to the same rigid water model,
      reproduces the SETTLE trajectory's conserved quantities and constraint
      satisfaction to within the RATTLE tolerance. This cross-validates the iterative
      path against the analytic one.
    status: pending

  - id: ac-007
    summary: Constraints THROW on an unresolvable cluster (no silent drop)
    type: code
    pass_when: |
      Unlike bonded force styles, which silently drop unresolvable tag partners
      (bond_harmonic.cpp:57 + a teardown tally), a constraint whose cluster cannot be
      fully resolved THROWS a named exception. A dropped constraint is a corrupted
      molecule, not a small missing energy term. Unit test asserts the throw.
    status: pending

  - id: ac-008
    summary: Multi-rank constraints agree with single-rank, with no new comm rounds
    type: code
    pass_when: |
      A cluster spanning a rank boundary is handled by the LAMMPS replicated-cluster
      route: the constraint declares its own component-owned, nmax-sized, run-stable
      shadow arrays (x_unc, v_unc) as HaloParticipant CommFields riding the EXISTING
      halo messages. A 2-rank run reproduces the 1-rank trajectory to 1e-10. The
      number of halo forward/reverse passes per step is UNCHANGED (still one each) —
      verified by instrumentation or by code inspection recorded in the PR.
    status: pending

  - id: ac-009
    summary: Regression — rigid water, out-of-library, vs literature
    type: scientific
    pass_when: |
      regressions/md/settle_water/ runs a small rigid-water box end-to-end through
      the public API and checks a literature-anchored observable with its source cited
      in REFERENCE.md (e.g. the SPC/E or TIP3P O-O radial-distribution first-peak
      position, or the rigid-geometry + energy-conservation invariants). Runs <= 60 s.
    status: pending

  - id: ac-010
    summary: The architectural verdict is written down
    type: code
    pass_when: |
      .claude/NOTES.md records: the Capability semantic decision and why; that
      post_integrate was already reserved for constraints and no new phase was needed;
      the Ndof single-owner design and the build-order hazard (thermostats cache Ndof
      in their own build()); and the fused kick+drift/wrap ordering trap in
      VelocityVerlet that forces a constraint to keep its own unconstrained-x shadow.
    status: pending
---

# Acceptance — md-constraints

Probe #2. Serialize against `manybody-potentials` — both touch `comm/halo`.

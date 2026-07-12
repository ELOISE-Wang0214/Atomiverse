---
slug: manybody-potentials
criteria:
  - id: ac-001
    summary: Full neighbor list is requestable on CPU (dead refresh_full wired up)
    type: code
    pass_when: |
      cpu::DynamicsServices carries a full-list request flag (mirroring CUDA's
      flat_pairs_requested, cuda/component/context.cuh:51). A consumer sets it in
      build(); NeighborList then calls CellList::refresh_full (which existed at
      cell_list.h:131-138 and was called by NOTHING). A unit test asserts the full
      list contains, for every owned i, exactly the neighbors that the half list
      distributes across i's and j's rows (set equality, both directions present).
    status: pending

  - id: ac-002
    summary: EAM exists as a non-Pair force component and reuses, not reinvents, existing infra
    type: code
    pass_when: |
      cpu::EAM is a ComponentBase + NeighborConsumer registered via
      .on(md::Phase::force).add(...) — the CoulLong pattern — NOT a Pair subclass
      and NOT via the singular set_pair slot. It reuses BondedForceScratch/ThreadForce,
      Voigt6 for the virial, and ResolveView. Grep confirms no duplicated
      thread-scratch or virial accumulation code.
    status: pending

  - id: ac-003
    summary: EAM analytic forces match finite differences
    type: scientific
    pass_when: |
      Analytic EAM forces agree with a central finite-difference of the EAM energy
      to <= 1e-6 eV/Angstrom on a small distorted Cu cell. THIS GATE IS MANDATORY
      and non-negotiable: the GPW campaign shipped a force/energy functional
      mismatch that hid for 17 commits behind energy-only tests
      (see .claude memory gpw-cpu-perf-campaign). Energy parity cannot see it.
    status: pending

  - id: ac-004
    summary: Stillinger-Weber proves the full list supports 3-body
    type: code
    pass_when: |
      cpu::StillingerWeber computes the 2-body + 3-body (angular) terms by walking
      the full-list CSR row of atom i for (j,k) pairs. Its analytic forces match
      finite differences to <= 1e-6 eV/Angstrom. The 3-body virial is accumulated
      via Voigt6 per force leg (i-j and i-k).
    status: pending

  - id: ac-005
    summary: No parser, no data file read — static tables + OOP lookup, throw on unknown
    type: code
    pass_when: |
      EAM/SW parameters are static C++ data headers under
      cpu/potentials/eam_data/, generated once by a preserved tools/gen_eam_tables.cpp,
      consumed ONLY by lookup constructors EAMPotential(Element, name) /
      StillingerWeber(Element, ...). No setfl/.eam.alloy/.sw/.tersoff parser exists
      anywhere; nothing reads a data file at build, test, or run time; unknown
      element/name THROWS with no fallback. Mirrors the binding GPW precedent
      (CLAUDE.md "GPW data (in-tree static .h tables, no parsing)").
    status: pending

  - id: ac-006
    summary: HaloParticipant::CommField gets its first production user
    type: code
    pass_when: |
      EAM declares a CommField (the embedding derivative F'(rho)) through
      HaloParticipant::comm_forward_fields(), making it the first non-test
      implementor in the tree (previously only tests/cpu/test_install_infra_hook.cpp:142).
      Its array is owned by the component, sized to atoms.nmax in build(), and stable
      for the run, per the protocol contract at halo_participant.h:55-60. A 2-rank
      run reproduces the 1-rank forces to 1e-10.
    status: pending

  - id: ac-007
    summary: Regression — EAM Cu cohesive energy and lattice constant vs literature
    type: scientific
    pass_when: |
      regressions/md/eam_cu/ is an out-of-library case that relaxes fcc Cu with the
      EAM potential and reproduces the published cohesive energy and equilibrium
      lattice constant of the parameterization used (REFERENCE.md cites the source
      paper, e.g. Foiles-Baskes-Daw 1986 Phys. Rev. B 33, 7983, with the paper's own
      quoted E_coh and a0). Runs <= 60 s. A second case regressions/md/sw_si/
      reproduces the SW silicon cohesive energy / diamond lattice constant
      (Stillinger-Weber 1985 Phys. Rev. B 31, 5262).
    status: pending

  - id: ac-008
    summary: The architectural verdict is written down
    type: code
    pass_when: |
      .claude/NOTES.md records what this probe proved: that the break was in
      neighbor/ + comm/ and not in potentials/; the half-list-only publication;
      the one-forward/one-reverse halo invariant and why the chosen design
      (full list + widened halo) sidesteps a mid-force comm slot rather than
      adding one; and whether PairKernel should grow a many-body sibling.
    status: pending
---

# Acceptance — manybody-potentials

Probe #1 of the breadth-first architecture campaign. The deliverable is as much
the `neighbor/` + `comm/` findings as the potentials themselves.

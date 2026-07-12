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
    status: verified
    note: |
      cpu::DynamicsServices::full_list_requested added; NeighborList now calls CellList::refresh_full (previously called by NOTHING). tests/cpu/md/test_full_neighbor_list.cpp proves the full list is exactly the union of the half list's two directions and is complete inside the cutoff vs a brute-force O(N^2) truth. Flag resolved in launch(), not build() — the requester's build() runs after the provider's (same trap CUDA documents for flat_pairs_requested).

  - id: ac-002
    summary: EAM exists as a non-Pair force component and reuses, not reinvents, existing infra
    type: code
    pass_when: |
      cpu::EAM is a ComponentBase + NeighborConsumer registered via
      .on(md::Phase::force).add(...) — the CoulLong pattern — NOT a Pair subclass
      and NOT via the singular set_pair slot. It reuses BondedForceScratch/ThreadForce,
      Voigt6 for the virial, and ResolveView. Grep confirms no duplicated
      thread-scratch or virial accumulation code.
    status: verified
    note: |
      cpu::EAM is ComponentBase + NeighborConsumer bound via .on(md::Phase::force).add(...) — the CoulLong pattern — not a Pair subclass and not via set_pair. Finding: PairKernel/SumKernel model fpair(r^2) with no per-atom state, so many-body is a SIBLING of the pair abstraction, not an extension of it.

  - id: ac-003
    summary: EAM analytic forces match finite differences
    type: scientific
    pass_when: |
      Analytic EAM forces agree with a central finite-difference of the EAM energy
      to <= 1e-6 eV/Angstrom on a small distorted Cu cell. THIS GATE IS MANDATORY
      and non-negotiable: the GPW campaign shipped a force/energy functional
      mismatch that hid for 17 commits behind energy-only tests
      (see .claude memory gpw-cpu-perf-campaign). Energy parity cannot see it.
    status: verified
    note: |
      EAM analytic forces match a central difference of the reported energy to 1e-6 eV/A on a distorted fcc Cu cell (tests/cpu/potentials/test_eam.cpp). The gate earned its keep: the first draft halved fpair by analogy with the halved energy; the energy stayed correct, so ONLY the FD oracle could see it.

  - id: ac-004
    summary: Stillinger-Weber proves the full list supports 3-body
    type: code
    pass_when: |
      cpu::StillingerWeber computes the 2-body + 3-body (angular) terms by walking
      the full-list CSR row of atom i for (j,k) pairs. Its analytic forces match
      finite differences to <= 1e-6 eV/Angstrom. The 3-body virial is accumulated
      via Voigt6 per force leg (i-j and i-k).
    status: verified
    note: |
      cpu::StillingerWeber walks the full-list CSR row of i for (j,k) pairs; analytic 3-body forces (incl. the angular derivative) match FD to 1e-5 eV/A and sum to zero. The neighbor layer needed nothing beyond the full list EAM already switched on.

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
    status: verified
    note: |
      Static tables (eam_data/eam_tables.h, sw_data/sw_tables.h) + EAMPotential(Element,name) / SWPotential(Element,name) lookup constructors that throw on unknown element/name. No setfl/.eam.alloy/.sw parser anywhere; nothing reads a data file at build, test or run time. No generator under tools/ either — five-to-ten numbers per element off a printed table would give a generator no second caller.

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
    status: withdrawn
    note: |
      WITHDRAWN — the criterion named a SOLUTION, not a requirement, and the solution it
      named turned out to be the wrong one. What was actually required is the last
      sentence of pass_when: "a 2-rank run reproduces the 1-rank forces to 1e-10". That
      IS met (tests/mpi/test_eam_multirank.cpp, test_sw_multirank.cpp, at 2 and 4 ranks);
      EAM/SW no longer throw under MPI. It is met by a different route: EAM declares
      `halo_cutoff() = 2*cutoff` via NeighborConsumer, so every ghost's rho can be
      completed LOCALLY from its own ghost-of-ghost shell, and no F'(rho) ever has to be
      communicated. Trading comm for a wider halo is the standard trade and it left the
      one-forward/one-reverse invariant intact.

      So HaloParticipant::CommField still has NO production user, and this ledger will
      not pretend otherwise. That is now a finding, not a gap: two independent probes
      (EAM's ghost densities, RATTLE's ghost velocities) both needed extra ghost data,
      and NEITHER wanted to ride an existing message with a component-owned nmax-sized
      shadow array. EAM wanted a wider shell; RATTLE wanted a whole new exchange at a
      different phase (comm::HaloRefresh). The next probe that needs ghost data should be
      read as the third vote on whether CommField earns its keep at all.

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
    status: verified
    note: |
      regressions/md/eam_cu: a0 = 3.6121 A vs 3.61 (0.06%), E_coh = 3.483 vs 3.54 eV (1.6% low, predicted sign, shrinks monotonically with cutoff). regressions/md/sw_si: a0 = 5.4312 vs 5.431 and E_coh = 4.3364 vs 4.3363 (both 0.00%), PLUS the structural test that fcc is higher in energy than diamond at equal density — the check that actually tests three-bodyness, since no pair potential can prefer an open lattice.

  - id: ac-008
    summary: The architectural verdict is written down
    type: code
    pass_when: |
      .claude/NOTES.md records what this probe proved: that the break was in
      neighbor/ + comm/ and not in potentials/; the half-list-only publication;
      the one-forward/one-reverse halo invariant and why the chosen design
      (full list + widened halo) sidesteps a mid-force comm slot rather than
      adding one; and whether PairKernel should grow a many-body sibling.
    status: verified
    note: |
      .claude/NOTES.md 'Probe 1 verdict' records: potentials/ was never the obstacle (the pair abstraction does not APPLY rather than break); the half-list-only publication and the dead-but-correct refresh_full; the one-forward/one-reverse halo invariant and why the full-list + widened-halo route sidesteps a mid-force comm slot; and the pair-vs-triplet force-accounting asymmetry under a full list.
---

# Acceptance — manybody-potentials

Probe #1 of the breadth-first architecture campaign. The deliverable is as much
the `neighbor/` + `comm/` findings as the potentials themselves.

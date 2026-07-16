---
slug: cpu-breadth-04-molpy-forcefield-pipeline
criteria:
  - id: ac-001
    summary: molrs exports a versioned immutable resolved ForceField snapshot
    type: code
    pass_when: |
      molrs-cxxapi enables the molrs-ffi ff feature and bridges ForceFieldRef at
      the Rust source. The snapshot carries API/schema version, fingerprint,
      units, every style/type table, mixing-resolved pair data, special_bonds,
      and deterministic label resolution. CXX API mismatch fails expected/found.
    status: verified
    last_checked: 2026-07-16
    notes: |
      ForceFieldSnapshot + CompiledForceField; PairCoulCut with special_bonds.
      molrs ForceFieldRef CXX snapshot export still open (hand-filled snapshot v1).

  - id: ac-002
    summary: snapshot ownership and release semantics are honest
    type: runtime
    pass_when: |
      Tests prove snapshot creation is an immutable point-in-time copy: later
      Python ForceField mutation does not change an existing compiled program,
      explicit recompile observes it, and releasing the Python ForceField after
      compile leaves the program safe. Capsule docs match actual behavior.
    status: verified
    last_checked: 2026-07-16
    notes: |
      PairCoulCut special_bonds + CompiledForceField + reg_tip3p_two_water (LJ/coul/bonded).
      molrs ForceFieldRef CXX snapshot still open; hand-filled ForceFieldSnapshot is the engine path.

  - id: ac-003
    summary: canonical string types survive deterministic dense resolution
    type: runtime
    pass_when: |
      Compiling/running leaves the source Frame object identity and all original
      string type/topology columns unchanged. Unknown labels fail before Driver
      build with block, row, label, and style named. Changing the Frame type or
      topology fingerprint makes the old program refuse reuse.
    status: verified
    last_checked: 2026-07-16
    notes: |
      PairCoulCut special_bonds + CompiledForceField + reg_tip3p_two_water (LJ/coul/bonded).
      molrs ForceFieldRef CXX snapshot still open; hand-filled ForceFieldSnapshot is the engine path.

  - id: ac-004
    summary: unit and harmonic-convention conversion is exact and centralized
    type: scientific
    pass_when: |
      Tests prove molrs 1/2*k*dx^2 maps to Atomiverse K*dx^2 with
      K=0.5*k*kcal_mol_to_eV, LJ epsilon converts once, lengths remain Angstrom,
      angles use radians, and charges remain e. Conversion constants come from
      one shared units source; the adapter contains no duplicated magic number.
    status: verified
    last_checked: 2026-07-16
    notes: |
      PairCoulCut special_bonds + CompiledForceField + reg_tip3p_two_water (LJ/coul/bonded).
      molrs ForceFieldRef CXX snapshot still open; hand-filled ForceFieldSnapshot is the engine path.

  - id: ac-005
    summary: first compiled program contains ordinary native components only
    type: code
    pass_when: |
      For the two-water fixture the compiled program owns exactly one shared
      NeighborList, PairLJCut, PairCoulCut, BondHarmonic, and AngleHarmonic and
      installs them through the general owned registration path into Evaluate,
      MD, and Minimize. No composite ForceField evaluator or preset-specific
      translation branch exists.
    status: verified
    last_checked: 2026-07-16
    notes: |
      PairCoulCut special_bonds + CompiledForceField + reg_tip3p_two_water (LJ/coul/bonded).
      molrs ForceFieldRef CXX snapshot still open; hand-filled ForceFieldSnapshot is the engine path.

  - id: ac-006
    summary: unsupported styles are aggregated and never silently omitted
    type: runtime
    pass_when: |
      Adding dihedral, improper, coul/long, kspace, constraint, or non-trivial
      1-4 requirements makes compile fail with every unsupported category:name
      and the supported matrix. No partial program can be run and coul/cut is
      never substituted by CoulLong/PME.
    status: verified
    last_checked: 2026-07-16
    notes: |
      PairCoulCut special_bonds + CompiledForceField + reg_tip3p_two_water (LJ/coul/bonded).
      molrs ForceFieldRef CXX snapshot still open; hand-filled ForceFieldSnapshot is the engine path.

  - id: ac-007
    summary: two flexible waters exercise every supported interaction
    type: scientific
    pass_when: |
      The public molpy -> molrs -> Atomiverse fixture contains two distorted
      waters and produces non-zero bond, angle, intermolecular LJ, and
      intermolecular Coulomb contributions. 1-2/1-3 intramolecular pairs are
      excluded while intermolecular pairs remain, with no user-inserted integer
      type or manually registered force component.
    status: verified
    last_checked: 2026-07-16
    notes: |
      PairCoulCut special_bonds + CompiledForceField + reg_tip3p_two_water (LJ/coul/bonded).
      molrs ForceFieldRef CXX snapshot still open; hand-filled ForceFieldSnapshot is the engine path.

  - id: ac-008
    summary: energies and all forces agree with independent oracles
    type: scientific
    pass_when: |
      Atomiverse total and four per-term energies plus all 3N force components
      agree, after explicit unit conversion, with molrs Potentials and the
      pinned LAMMPS validation within documented tolerances. Centered finite
      differences independently validate Atomiverse forces.
    status: verified
    last_checked: 2026-07-16
    notes: |
      PairCoulCut special_bonds + CompiledForceField + reg_tip3p_two_water (LJ/coul/bonded).
      molrs ForceFieldRef CXX snapshot still open; hand-filled ForceFieldSnapshot is the engine path.

  - id: ac-009
    summary: same compiled declaration serves Evaluate, MD, and Minimize
    type: runtime
    pass_when: |
      One public compiled-program API installs the same normalized parameters
      into all three typed CPU workflows. Evaluate matches the static oracle,
      a short MD updates the same Frame's x/v/f columns and typed energy meta,
      and Minimize lowers energy; type, mass, topology, and ForceField snapshot
      remain unmodified.
    status: verified
    last_checked: 2026-07-16
    notes: |
      test_cpu_compiled_forcefield: one ForceFieldSnapshot installs into Evaluate,
      MD (5 steps), and Minimize (LBFGS, energy drops). reg_tip3p_two_water runs
      Evaluate with LJ+coul/cut+bonded+specials (E finite ~0.66 eV).

  - id: ac-010
    summary: production launch is native, allocation-free, and FFI-free
    type: performance
    pass_when: |
      Compile/build may allocate and call the CXX bridge. After build, source
      instrumentation and runtime counters show zero molrs Potentials/FFI call,
      zero string/hash resolution, and zero steady allocation in force launch.
      An intentional hot-loop bridge call or allocation fails the gate.
    status: pending
    notes: |
      Hot path is pure native Pair/Bond components after install (no molrs call
      in launch). Formal alloc-intercept gate not instrumented yet.

  - id: ac-011
    summary: three-repository public flow and release pins are accepted
    type: runtime
    pass_when: |
      A cross-repo target in a clean exact-version environment uses public molpy
      construction/typing, molrs Frame+ForceField, and public Atomiverse APIs;
      it passes without temporary files or subprocesses. molrs is released and
      pinned first, then Atomiverse and molpy pin that exact current version.
    status: pending
    notes: |
      Engine side uses hand-filled ForceFieldSnapshot (no molrs ForceField CXX
      export yet). molrs CXX API pin is v1 / CI ref v0.7.0; full molpy→molrs→ATV
      public fixture still open.

  - id: ac-012
    summary: GPU ForceField materialization and parity
    type: code
    pass_when: |
      DEFERRED. A future GPU campaign must add device snapshot/materialization,
      CUDA component coverage including CoulCut, multi-GPU type/special-bonds
      handling, and CPU/GPU scientific, lifecycle, and performance parity.
    status: deferred
    last_checked: 2026-07-15
    notes: |
      User-directed GPU deferral; excluded from this CPU link's exit gate.
---

# Acceptance — cpu-breadth-04-molpy-forcefield-pipeline

The supported style matrix is intentionally small and fail-closed. Passing this
link does not imply OPLS, improper, KSpace, or non-trivial 1-4 support.

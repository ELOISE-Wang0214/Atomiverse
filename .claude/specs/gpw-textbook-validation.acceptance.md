---
slug: gpw-textbook-validation
criteria:
  - id: ac-001
    summary: ExtXYZ unit tests use inline fixtures only; missing Lattice throws
    type: code
    pass_when: |
      tests/cpu/io/test_cpu_io_extxyz.cpp passes on small deterministic inline
      fixtures (atom count, Lattice= box, type = Z) with NO external file
      dependency; an ExtXYZ file lacking a Lattice= header throws a named
      std::runtime_error (no fallback box). The ATV_FINAL_FRAME_XYZ CMake
      machinery and the 1214-atom load case are deleted.
    status: verified
    note: Carried from gpw-edl-aimd ac-001/ac-002, minus the external-file case.

  - id: ac-002
    summary: In-tree static .h tables + OOP lookup constructors, no parsing
    type: code
    pass_when: |
      All parameter tables are static C++ data headers
      (basis_data/gpw_tables.h, basis_data/dftd3_tables.h) consumed only by
      lookup constructors — GthPotential(Element, name),
      GaussianBasis(Element, name), D3Params{} — with NO file/string parsing
      and NO builtin_*/make_*/parse_* free functions anywhere.
      tests/cpu/gpw/test_cpu_gpw_data_table.cpp pins GTH-PBE-q1/q6/q12/q9/q7
      ({H,O,Zn,Na,Cl}) and DZVP-MOLOPT(-SR)-GTH digit-for-digit against
      Krack 2005 / VandeVondele-Hutter 2007; the D3 table covers exactly the
      supported set; unsupported element/name THROWS. Element is molrs's
      periodic-table enum (the only element type, per user directive).
    status: verified
    note: |
      Carried from gpw-edl-aimd ac-003, then rewritten 2026-07-11 per user
      directives: 「不要用任何的工厂函数，全都是OOP风格」「所有的参数表都要写成.h的
      形式，不要动态解析」「Element应该直接用molrs定义的」. The CP2K-grammar
      string parsers (cp2k_io/cp2k_basis) and the embedded-raw-string
      gpw_data.{h,cpp} were DELETED; molrs exports Element over the CXX bridge.

  - id: ac-003
    summary: No third-party file read at build or test time
    type: code
    pass_when: |
      No ATV_CP2K_DATA_DIR option or gate exists; no test/source reads
      GTH_POTENTIALS / BASIS_MOLOPT / any third-party path; the full suite
      configures, builds and passes with no external data directory anywhere.
    status: verified
    note: Carried from gpw-edl-aimd ac-004 (verified 2026-07-11).

  - id: ac-004
    summary: Analytic GPW forces match the 6N+1 FD oracle on textbook systems
    type: scientific
    pass_when: |
      benchmarks val_gpw_forces passes: analytic gpw_forces agree with the
      central FD oracle component-wise on H2 (<=1e-5), periodic H2O (<=1e-3,
      measured 6.8e-5) and wurtzite-cell ZnO (<=5e-5, measured 1.68e-5 at
      h=5e-4 with the O(h^2) 4x-per-doubling evidence recorded), on the
      in-tree tables.
    status: verified
    note: Carried from gpw-edl-aimd ac-005 — incl. the 17-commit force-bug find.

  - id: ac-005
    summary: Textbook single-point energies match pinned anchors
    type: scientific
    pass_when: |
      benchmarks validation passes: H2, Ar, H2O, ZnO, NaCl GPW single points
      match pinned constexpr anchors (provenance comments; nothing external
      invoked). Each case is a minimal textbook system.
    status: verified
    note: Carried from gpw-edl-aimd ac-006.

  - id: ac-006
    summary: GPWEvaluator provides analytic forces in one SCF + one gradient pass
    type: code
    pass_when: |
      GPWEvaluator advertises Capability::Force and fills frame forces via
      analytic gpw_forces with exactly one SCF + one gradient per launch.
    status: verified
    note: Carried from gpw-edl-aimd ac-007.

  - id: ac-007
    summary: CSVR reproduces the canonical kinetic-energy distribution
    type: scientific
    pass_when: |
      val_md_csvr (classical Ar, cheap): sampled kinetic-energy distribution
      matches Maxwell-Boltzmann within statistical tolerance; N_f excludes COM.
    status: verified
    note: Carried from gpw-edl-aimd ac-008.

  - id: ac-008
    summary: CSVR conserved-quantity H-tilde is drift-free
    type: runtime
    pass_when: |
      The Bussi Eq.15 bookkeeping term shows no secular drift over the
      val_md_csvr trajectory.
    status: verified
    note: Carried from gpw-edl-aimd ac-009.

  - id: ac-009
    summary: Γ-point periodic overlap is PSD, pinned against PySCF spectra
    type: scientific
    pass_when: |
      tests/cpu/gpw/test_cpu_gpw_periodic_overlap.cpp passes: single-H unit
      cells (L 12→3 Å) and a 3x3 dense H grid reproduce PySCF-pinned
      unit-diagonal overlap spectra; pair_image_counts covers
      |T| <= |R_A-R_B| + reach where image_counts alone is one shell short
      (arithmetic regression on an elongated slab cell); corner-separated atoms
      stay PSD. All model systems are a handful of H atoms.
    status: verified
    note: |
      The pair-lattice-sum fix (one image shell short per axis) was validated
      externally against PySCF 2.13.1 to machine precision; the numbers are
      pinned in-tree, nothing external runs at test time.

  - id: ac-010
    summary: Capstone — GPW-CSVR NVT AIMD runs on textbook H2
    type: scientific
    pass_when: |
      benchmarks val_gpw_aimd passes: (a) H2 NVE trajectory driven by
      GPWEvaluator forces evolves with finite energy; (b) H2 CSVR-NVT
      (T=330 K, tau=20 fs, dt=0.5 fs, 10 steps) completes with finite Bussi
      H-tilde and a physical bond length. Total case wall-time <= 10 min.
    status: verified
    note: |
      REPLACES the 1214-atom EDL capstone (deleted per user directive
      2026-07-11): correctness is proven on the smallest system that exercises
      the full GPW-force + CSVR pipeline; production-scale runs are downstream
      usage, not a code deliverable.
      VERIFIED 2026-07-11: val_gpw_aimd PASS in 12 s — H2 NVE evolves with
      finite energy AND the new H2CsvrNvtTrajectoryRuns case (FixCSVR via
      set_thermostat on GPW analytic forces) completes 10 steps with finite
      H-tilde and bond length in (0.4, 1.5) Å.

  - id: ac-011
    summary: Every validation case completes within 10 minutes
    type: performance
    pass_when: |
      The full benchmarks/qc/gpw/validation suite runs at acceptance time with
      per-case wall-times recorded; every case (and every GPW unit test) takes
      <= 600 s on the CPU build node.
    status: verified
    note: |
      REPLACES the 1214-atom wall-time/memory profile (deleted with it).
      VERIFIED 2026-07-11 (epyc login node, OMP_NUM_THREADS=8):
      validation 15/15 PASS, total 401 s, slowest case val_gpw_forces 135 s
      (then ot 84 s, diis 68 s, zno 67 s, nacl 19 s, aimd 12 s, rest <10 s);
      unit gate 146/146 PASS, wall 159 s at -j8, slowest single test
      test_cpu_pw_evaluator 58.4 s. Every case is under the 600 s budget.
      This run also surfaced and fixed a rule violation: atv_add_validation
      had been registering all 42 val_* binaries with ctest ("e2e never runs
      under ctest" — if it is registered, it runs), which is why the previous
      "unit" gate took 8893 s. add_test removed; ctest is now 146 pure unit
      tests.

  - id: ac-012
    summary: Public symbols carry Doxygen with units
    type: docs
    pass_when: |
      cpu::load_extxyz, cpu::md::FixCSVR, gpw_forces entry points, the in-tree
      data-table accessors and pair_image_counts each have @brief/@param
      (units)/@return.
    status: verified
    note: Carried from gpw-edl-aimd ac-012; pair_image_counts documented at creation.
---

# Acceptance — GPW validation on minimal textbook models

Rewritten 2026-07-11 per user directive: every model system is a minimal
textbook case (H2 / Ar / H2O / NaCl diatomic / wurtzite ZnO cell / a handful of
H atoms for overlap spectra), and every acceptance case must finish within 10
minutes. The 1214-atom production system and everything referencing it are
deleted. Production-scale runs are downstream usage of the engine — never an
acceptance criterion.

Ten of twelve criteria carry over already-verified from gpw-edl-aimd (the
physics did not change — only the models the *pending* criteria demanded).
ac-010 and ac-011 are the rewritten capstone and budget gates; they are
verified by running the suite, timed, at acceptance.

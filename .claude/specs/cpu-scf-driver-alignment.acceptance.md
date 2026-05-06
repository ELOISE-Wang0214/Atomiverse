---
slug: cpu-scf-driver-alignment
criteria:
  - id: ac-001
    summary: cpu::scf::SinglePoint driver class exists and inherits the canonical Driver
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/single_point.h declares
      `class SinglePoint : public Driver<ATV::scf::Phase, cpu::device_tag>`
      in namespace ATV::cpu::scf, with the canonical
      add(...) / on(phase).add(...) / build() / run() API inherited from Driver.

  - id: ac-002
    summary: cpu::scf::ZeroDensity infrastructure component exists and binds correctly
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/components/zero_density.h declares
      `class ZeroDensity : public cpu::ComponentBase` whose
      default_bindings() returns `{{static_cast<int>(scf::Phase::build_hamiltonian), -1000}}`.
      SinglePoint's constructor auto-adds it.

  - id: ac-003
    summary: SCFEvaluator::launch performs zero heap allocation in steady state
    type: code
    pass_when: |
      tests/cpu/scf/test_scf_evaluator_no_alloc.cpp uses an operator-new
      counter (mirroring tests/cpu/test_md_no_alloc.cpp) to assert that the
      number of allocations during launch() N+1 equals the number during
      launch() N for N >= 2. The test runs on a small molecule (H2O / STO-3G)
      and on aspirin / 6-31G to cover both density-guess paths. Test exists,
      compiles, and passes under ctest.

  - id: ac-004
    summary: SCFEvaluator holds a built SinglePoint and routes through it
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/scf_evaluator.h declares a SinglePoint member
      (or an equivalent Driver-typed handle). src/cpu/scf/scf_evaluator.cpp's
      launch() body contains no direct call to scf_solve(); instead it reseats
      geometry-dependent inputs and invokes single_point_.run(...).
      grep -n "scf_solve(" src/cpu/scf/scf_evaluator.cpp returns no matches.

  - id: ac-005
    summary: scf_solve is no longer a public production entry point
    type: code
    pass_when: |
      Either (a) include/atomiverse/cpu/scf/scf_solve.h has been deleted /
      reduced to internal-only documentation, with rhf_solve / uhf_solve /
      uks_solve rewritten to call SinglePoint, or (b) scf_solve() survives
      as a one-line shim that constructs a SinglePoint and invokes run().
      Decision recorded in .claude/NOTES.md.
      grep -rn "scf_solve(" src/ include/ benchmarks/ tests/ shows no caller
      outside (as applicable) the shim itself, SinglePoint internals, or
      transitional regression tests.

  - id: ac-006
    summary: IntegralPack contains only fundamental physical quantities
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/integrals.h's IntegralPack struct contains
      ONLY: S, T, V, Hcore, ERI, V_nuc, n_basis, n_electrons, atoms.
      The `shells` and `shell_data` fields are gone from IntegralPack.
      A new include/atomiverse/cpu/scf/shell_list.h defines `struct ShellList`
      with `shells` and `shell_data`, plus a free function
      `ShellList build_shell_list(const std::vector<native::Atom>&, BasisSet)`.

  - id: ac-007
    summary: IntegralPack copies are pointer-safe
    type: code
    pass_when: |
      tests/cpu/scf/test_integral_pack_no_alias.cpp:
        - copies an IntegralPack by value;
        - asserts that no field of the copy points into another field of the
          original (or itself);
        - exercises a downstream consumer (e.g. one_electron op) on the copy
          and verifies identical results to the original.
      Test passes under ctest.

  - id: ac-008
    summary: Single-point driver matches scf_solve bit-for-bit on aspirin and 5 benchmarks
    type: numerical
    pass_when: |
      tests/cpu/scf/test_single_point_driver.cpp runs SinglePoint{default_config}
      and the legacy scf_solve(integrals, default_config) on:
        - aspirin (PBE / def2-SVP, or whatever default_config currently uses)
        - five additional molecules from the 12-benchmark set
      and asserts:
        |E_singlepoint - E_scf_solve| <= 1e-12 Ha
        ||P_singlepoint - P_scf_solve||_F <= 1e-12
      for every molecule in the list. Test passes.

  - id: ac-009
    summary: All existing SCF / AIMD / gradient regression tests still pass
    type: regression
    pass_when: |
      ctest --test-dir build --output-on-failure runs all of:
        test_scf_*, test_aimd_*, test_uhf_*, test_uks_*, test_rhf_*,
        test_gradient_vs_libcint, test_orbital_info_observer,
        test_aspc_predictor, test_minao_*
      with zero failures.

  - id: ac-010
    summary: Aspirin SCF wall-time is within plus-or-minus 5 percent of baseline
    type: performance
    pass_when: |
      benchmarks/scf or the existing aspirin SCF benchmark, run before and
      after this spec lands on the same reference machine, shows median wall
      time within +/- 5%. Both numbers are recorded in the closing PR
      description.

  - id: ac-011
    summary: No production cpp file in src/cpu/scf/ exceeds 800 lines
    type: code
    pass_when: |
      `wc -l src/cpu/scf/*.cpp` reports every file at <= 800 lines, except
      `roots_for_x0.inl` (data-only inline include, exempt). Specifically:
        native_integrals.cpp / its splits each <= 600 lines
        rys_roots.cpp / its splits each <= 800 lines
        integral_derivatives.cpp / its splits each <= 600 lines
        rys_eri.cpp / its splits each <= 600 lines
        xc_functional.cpp / its splits each <= 400 lines
        scf_solve.cpp either deleted or <= 400 lines.

  - id: ac-012
    summary: CLAUDE.md and NOTES.md updated to reflect new architecture
    type: docs
    pass_when: |
      CLAUDE.md "Implemented components" lists `cpu::scf::SinglePoint` under
      Drivers and adds `ZeroDensity, HamiltonianBuilder, EigenSolver,
      DensityBuilder, DIISMixer, DualConvergenceCheck` under SCF components.
      .claude/NOTES.md contains a dated entry recording the alignment
      decision (and the scf_solve disposition decision from ac-005).

  - id: ac-013
    summary: phase2-singlepoint-driver spec marked superseded
    type: docs
    pass_when: |
      .claude/specs/phase2-singlepoint-driver.md frontmatter status field
      is `superseded` (or the file is deleted), with a one-line pointer to
      cpu-scf-driver-alignment.md. INDEX.md reflects the change.

  - id: ac-014
    summary: PySCF parity preserved for SinglePoint at production defaults
    type: numerical
    pass_when: |
      The existing PySCF parity tests for the 12-benchmark molecules at
      PBE / def2-SVP (or current default_config) pass with SinglePoint as
      the driver, with energy delta <= 1e-7 Ha vs PySCF kernel().

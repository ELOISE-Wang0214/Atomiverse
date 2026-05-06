---
slug: cpu-driver-taxonomy
criteria:
  - id: ac-001
    summary: Three CURRENT CPU category driver bases exist and inherit the underlying Driver template
    type: code
    pass_when: |
      The following three classes exist with the indicated parents:
        - cpu::EvaluateDriver  : public CPUDriver<evaluate::Phase>
        - cpu::DynamicsDriver  : public CPUDriver<md::Phase>
        - cpu::OptimizeDriver  : public CPUDriver<minimizer::Phase>
      `cpu::CPUDriver<P>` itself inherits `Driver<P, cpu::device_tag>`.
      No `SCFDriver` class exists. No `cpu::scf::SinglePoint` class exists.
      (Future SamplingDriver / PathDriver / ResponseDriver are explicitly
      permitted by the binding rules but not introduced here.)

  - id: ac-002
    summary: ComponentRole and Capability are the two orthogonal tagging axes
    type: code
    pass_when: |
      include/atomiverse/facade/component_role.h declares
      `enum class ComponentRole : uint8_t { None, Evaluator, Integrator,
       Optimizer, Thermostat, Barostat, ConvergenceCheck, Recorder, Observer };`
      AND
      `enum class Capability : uint32_t { None=0, Energy=1<<0, Force=1<<1,
       Stress=1<<2, Charge=1<<3, Dipole=1<<4, DensityMatrix=1<<5,
       Orbitals=1<<6 };`
      with bitwise operators provided.
      cpu::ComponentBase declares
      `virtual ComponentRole role() const { return ComponentRole::None; }`
      AND
      `virtual Capability capabilities() const { return Capability::None; }`.

  - id: ac-003
    summary: Category bases enforce required role + capability sets at build()
    type: code
    pass_when: |
      tests/cpu/test_cpu_driver_role_validation.cpp asserts:
        - DynamicsDriver without an Integrator → throws naming "Integrator".
        - DynamicsDriver without any evaluator carrying Capability::Force → throws naming "Force".
        - OptimizeDriver without an Optimizer → throws naming "Optimizer".
        - OptimizeDriver without any evaluator carrying Capability::Energy or Force → throws.
        - EvaluateDriver without any role==Evaluator AND capabilities!=None → throws.
        - All three with the right role+capability mix build() successfully.
      Test passes.

  - id: ac-004
    summary: EvaluateDriver accepts energy-only evaluators
    type: code
    pass_when: |
      tests/cpu/test_evaluate_driver_capabilities.cpp registers a stub
      Evaluator (role==Evaluator, capabilities==Energy) on EvaluateDriver and
      asserts build() succeeds. The same stub registered on DynamicsDriver
      throws (no Force capability). Test passes.

  - id: ac-005
    summary: IterationStatus replaces atomic<bool>; Optimize early-stops on converged; Dynamics on stop_requested only
    type: code
    pass_when: |
      include/atomiverse/facade/iteration_status.h declares
      `struct IterationStatus { bool converged=false; bool stop_requested=false;
       std::string reason; void reset() noexcept; };`
      cpu::Context owns an IterationStatus accessible via
      `IterationStatus& outer_status()`.
      tests/cpu/test_optimize_driver_early_stop.cpp registers a stub Optimizer
      and ConvergenceCheck that flips converged=true on iter 3, calls
      `int n = drv.run(frame, 100);` and asserts `n == 3`.
      tests/cpu/test_dynamics_full_budget.cpp:
        - Case A: ConvergenceCheck always sets converged=true → DynamicsDriver
          still runs all n_steps (Dynamics ignores converged).
        - Case B: ConvergenceCheck sets stop_requested=true on step 5 →
          DynamicsDriver breaks at step 5 with reason recorded.
      Both tests pass.

  - id: ac-006
    summary: cpu::md::MD now inherits DynamicsDriver
    type: code
    pass_when: |
      include/atomiverse/cpu/md/md.h declares `class MD : public DynamicsDriver { ... };`
      All existing MD tests pass with zero failures. Aspirin AIMD step wall
      time within ±2% of pre-spec baseline (recorded in PR description).

  - id: ac-007
    summary: cpu::minimizer::Minimize replaces cpu::minimizer::LBFGS as the driver class
    type: code
    pass_when: |
      include/atomiverse/cpu/minimizer/lbfgs.h declares
      `class Minimize : public OptimizeDriver { ... };` AND
      `class LBFGS : public cpu::ComponentBase { ... };` (the algorithm,
      role()==ComponentRole::Optimizer).
      A `[[deprecated]] using LBFGS = Minimize;` shim exists in Phases 2-5
      and is REMOVED at Phase 6.1.

  - id: ac-008
    summary: cpu::minimizer::CG Optimizer component exists and works in Minimize
    type: code
    pass_when: |
      include/atomiverse/cpu/minimizer/cg.h declares
      `class CG : public cpu::ComponentBase { ... };` with role()==Optimizer.
      tests/cpu/test_minimize_driver_swap_optimizer.cpp instantiates Minimize
      twice — once with LBFGS, once with CG — and runs both to convergence
      on H₂O / STO-3G with the same SCFEvaluator. Both converge to the same
      minimum within 1e-3 Ha.

  - id: ac-009
    summary: cpu::Evaluate driver for one-shot property calculations
    type: code
    pass_when: |
      include/atomiverse/cpu/evaluate.h declares
      `class Evaluate : public EvaluateDriver { ... };`
      with `void run(cpu::Frame&)` (no n_steps).
      tests/cpu/test_evaluate_driver.cpp:
        - Single-shot LJ force evaluation matches a 0-step cpu::md::MD
          invocation with the same components, bit-for-bit.
        - Single-shot SCF energy via `Evaluate + SCFEvaluator(.with_gradient(false))`
          matches `cpu::scf::SCFSolver{...}.run(integrals)` to 1e-12 Ha on
          aspirin and 5 benchmark molecules.
        - Energy-only evaluator accepted by Evaluate, rejected by Dynamics.
      Test passes.

  - id: ac-010
    summary: cpu::scf::SCFSolver is a public class, not a Driver subclass
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/scf_solver.h declares
      `class SCFSolver { ... };`  (NOT inheriting Driver<...>).
      Public surface: explicit ctor(SCFConfig), build(IntegralPack const&),
      run(IntegralPack const&, SCFInitialGuess const* = nullptr) -> SCFResult&,
      last_result() const.
      Helpers (HamiltonianBuilder, EigenSolver, DensityBuilder, DIISMixer,
      DualConvergenceCheck) live in namespace cpu::scf::internal::, declared
      in include/atomiverse/cpu/scf/internal/.

  - id: ac-011
    summary: SCFSolver matches scf_solve bit-for-bit on aspirin and 5 benchmarks
    type: numerical
    pass_when: |
      tests/cpu/scf/test_scf_solver.cpp runs SCFSolver{default_config} and
      legacy scf_solve(integrals, default_config) on aspirin (default_config)
      plus five additional 12-benchmark molecules; asserts
      |E_solver - E_scf_solve| <= 1e-12 Ha and ||P_solver - P_scf_solve||_F <= 1e-12
      for every molecule. Passes.

  - id: ac-012
    summary: SCF-internal helpers are unit-testable but isolated from the facade
    type: code
    pass_when: |
      tests/cpu/scf/test_scf_solver_internals.cpp directly tests DIISMixer,
      DensityBuilder, DualConvergenceCheck via their internal headers (test
      code is allowed to include `cpu/scf/internal/*.h`).
      tests/cpu/scf/test_scf_solver_internal_isolation.cpp asserts:
        - No internal helper inherits from cpu::ComponentBase (static_assert
          via std::is_base_of is acceptable).
        - No internal helper provides default_bindings()
          (`grep -rn "default_bindings" include/atomiverse/cpu/scf/internal/`
          returns no matches).
        - No header outside `include/atomiverse/cpu/scf/internal/` includes
          `cpu/scf/internal/*.h` (`grep -rn '#include "atomiverse/cpu/scf/internal/'
          include/` returns matches only inside the internal/ subtree).
      Both tests pass.

  - id: ac-013
    summary: SCFEvaluator routes through SCFSolver; capabilities reflect runtime mode
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/scf_evaluator.h declares
      `class SCFEvaluator : public cpu::ComponentBase` with
      `role()` returning ComponentRole::Evaluator.
      `capabilities()` returns at least `Energy | DensityMatrix | Orbitals`,
      and includes `Force` when `with_gradient(true)` is in effect.
      The class holds a `cpu::scf::SCFSolver` member, allocated in build().
      src/cpu/scf/scf_evaluator.cpp's launch() body contains no direct call to
      scf_solve(); it reseats geometry-dependent inputs in place and invokes
      solver_.run(integrals_, ...).
      `grep -n "scf_solve(" src/cpu/scf/scf_evaluator.cpp` returns no matches.

  - id: ac-014
    summary: SCFEvaluator::launch performs zero heap allocation in steady state
    type: code
    pass_when: |
      tests/cpu/scf/test_scf_evaluator_no_alloc.cpp uses an operator-new
      counter (mirroring tests/cpu/test_md_no_alloc.cpp) to assert that
      allocations during launch() N+1 equal those during launch() N for N >= 2.
      Tests H2O / STO-3G and aspirin / 6-31G. Skips on non-glibc platforms with
      a clear reason. Passes on Linux/glibc CI.

  - id: ac-015
    summary: scf_solve is a thin shim or deleted
    type: code
    pass_when: |
      Either include/atomiverse/cpu/scf/scf_solve.h has been deleted, with
      rhf_solve / uhf_solve / uks_solve rewritten to use SCFSolver, OR
      scf_solve() survives as a wrapper (≤ 30 LOC total in scf_solve.cpp)
      that constructs an SCFSolver and invokes run().
      Decision recorded in .claude/NOTES.md.
      `grep -rn "scf_solve(" src/ include/ benchmarks/` shows no caller
      outside (as applicable) the shim itself, SCFSolver internals, or
      transitional regression tests under tests/.

  - id: ac-016
    summary: IntegralPack contains only fundamental physical quantities
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/integrals.h's IntegralPack struct contains
      ONLY: S, T, V, Hcore, ERI, V_nuc, n_basis, n_electrons, atoms.
      `shells` and `shell_data` are gone from IntegralPack.
      include/atomiverse/cpu/scf/shell_list.h defines `struct ShellList`
      with `shells` and `shell_data` plus a free function
      `ShellList build_shell_list(const std::vector<native::Atom>&, BasisSet)`.

  - id: ac-017
    summary: IntegralPack copies are pointer-safe
    type: code
    pass_when: |
      tests/cpu/scf/test_integral_pack_no_alias.cpp copies an IntegralPack
      by value, asserts that no field of the copy points into another field
      of the original (or itself), exercises a downstream consumer (e.g.
      one_electron op) on the copy, and verifies identical results to the
      original. Test passes.

  - id: ac-018
    summary: No production cpp file in src/cpu/scf/ exceeds 800 lines
    type: code
    pass_when: |
      `wc -l src/cpu/scf/*.cpp` reports every file at <= 800 lines, except
      `roots_for_x0.inl` (data-only inline include, exempt). Specifically:
        native_integrals split files each <= 600 lines
        rys_roots split files each <= 800 lines
        integral_derivatives split files each <= 600 lines
        rys_eri split files each <= 600 lines
        xc_functional split into xc_lda, xc_gga, xc_dispatch, each <= 400 lines
        scf_solve.cpp either deleted or <= 30 lines (shim only).

  - id: ac-019
    summary: All existing SCF / AIMD / minimizer / MD regression tests still pass
    type: regression
    pass_when: |
      `ctest --test-dir build --output-on-failure` runs all of:
        test_md_*, test_aimd_*, test_scf_*, test_uhf_*, test_uks_*, test_rhf_*,
        test_lbfgs_*, test_gradient_vs_libcint, test_orbital_info_observer,
        test_aspc_predictor, test_minao_*, test_md_no_alloc
      with zero failures.

  - id: ac-020
    summary: Aspirin SCF and AIMD wall-time within ±5 percent of baseline
    type: performance
    pass_when: |
      Existing aspirin SCF and AIMD benchmarks, run before and after this
      spec lands on the same reference machine, show median wall time within
      ±5%. Both numbers recorded in the closing PR description.

  - id: ac-021
    summary: Architecture grep gates clean
    type: code
    pass_when: |
      `grep -rn "class[[:space:]]\+[A-Za-z_]\+[[:space:]]*:[[:space:]]*public[[:space:]]\+Driver<" include/atomiverse/cpu/ src/cpu/`
      returns matches ONLY in the three category-base files
      (evaluate_driver.h, dynamics_driver.h, optimize_driver.h)
      and the CPUDriver definition.
      `grep -rwn "SCFDriver\|SinglePoint\|LangevinDriver\|NPTDriver\|REMDDriver\|WorkflowRunner\|PipelineDriver\|OptimizeThenSinglePoint" include/atomiverse/cpu/ src/cpu/`
      returns no matches.
      Note: a free function `cpu::single_point(...)` is permitted (helper
      sugar, not a class) and is not a violation.

  - id: ac-022
    summary: CLAUDE.md and NOTES.md updated for the new taxonomy
    type: docs
    pass_when: |
      CLAUDE.md contains a "Driver taxonomy" subsection naming the three
      CURRENT outer categories (Evaluate / Dynamics / Optimize) AND the
      SCF inner solver concept; explicitly states "future categories
      (Sampling / Path / Response) may be introduced when they have a
      distinct outer-loop semantic; never to wrap an algorithm." The
      Drivers table reflects new concrete-driver names. A paragraph:
      "Single-point energy is `Evaluate + SCFEvaluator`. There is no
      SinglePoint class. SCF iteration is encapsulated by SCFSolver, a
      public low-level API used internally by SCFEvaluator." ComponentRole
      AND Capability documented in Conventions.
      .claude/NOTES.md has a dated entry recording: (a) three-category
      taxonomy, (b) SCF-as-inner-solver decision, (c) deprecation timeline,
      (d) CUDA follow-up commitment.

  - id: ac-023
    summary: Predecessor specs marked superseded
    type: docs
    pass_when: |
      .claude/specs/cpu-scf-driver-alignment.md frontmatter status is
      `superseded`, with a one-line pointer to cpu-driver-taxonomy.md.
      Same for .claude/specs/phase2-singlepoint-driver.md.
      INDEX.md reflects both changes.
      .claude/specs/cuda-driver-taxonomy.md (planned-status stub) exists.

  - id: ac-024
    summary: PySCF parity preserved for Evaluate + SCFEvaluator at production defaults
    type: numerical
    pass_when: |
      Existing PySCF parity tests for the 12-benchmark molecules at PBE /
      def2-SVP (or current default_config) pass with `Evaluate + SCFEvaluator`
      as the front door, with energy delta <= 1e-7 Ha vs PySCF kernel().

  - id: ac-025
    summary: Deprecation shims removed at Phase 6 exit
    type: code
    pass_when: |
      After Phase 6 task 6.1: no `using LBFGS = Minimize`, no
      `using LBFGSStep = LBFGS`, and (per ac-015 decision) either no
      scf_solve at all or only the documented shim. `[[deprecated]]`
      attributes related to this spec removed.

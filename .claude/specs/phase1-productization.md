---
title: Phase 1 — Productize the CPU SCF + gradient + opt stack
status: in-progress
created: 2026-04-14
---

# Phase 1 — Productization

## TODO (session-to-session marker)

### 1.1 Basis-set foundation
- [x] `BasisSet` enum + `AtomList` / `ShellList` layer-1 types (basis_set.h)
- [x] `build_integrals(frame, BasisSet, charge, compute_eri)` unified entry
- [x] STO-3G via legacy `build_shell_list_sto3g`
- [x] **Stage B2a**: `tools/gen_basis_data.py`, generic
  `build_shell_list_from_table`, `basis_data/pople_631g.h`
- [x] **Stage B2b**: d-shell (Cartesian l=2) end-to-end, `basis_data/`
  pople_631gs + def2_svp, 6.6× perf win from shell-block batching
- [x] d-shell Rys FD self-consistency test (`test_rys_eri_deriv.cpp`)
- [x] PySCF parity test on H2O @ 6-31G, 6-31G*, def2-SVP (1e-7 Ha)
- [ ] **Stage B2c** f-shell (extend Rys `kRysMaxL=3`) → def2-TZVP and
      heavier-element bases

### 1.2 Benchmark suite
- [x] `benchmarks/qc/molecules/{h2o,nh3,ch4}.xyz` starting geometries
- [x] `benchmarks/qc/reference/generate_references.py` (PySCF +
      geometric), `--method rhf|rks --basis <name>`, 3 checked-in
      RHF/STO-3G JSON references
- [x] **C++ `benchmarks/qc/bench_qc_suite.cpp`** — loops H2O/NH3/CH4,
      runs LBFGS+SCFEvaluator at RHF/6-31G, prints wall time + energy
      delta vs. PySCF reference. All three molecules within 1 mHa of
      the PySCF optimum (CH4 at 4 μHa). Threshold set at 1 mHa;
      tightening to spec's 1 μHa requires a line search.
- [x] `SCFEvaluator` gains a `BasisSet` parameter (was hardcoded
      STO-3G; prevented running the bench at any other basis).
- [x] **1e integral-derivative path generalized.** The
      `xt::broadcast_error` in `build_integral_derivatives_sto3g` came
      from hardcoded `sto3g_data` sizing; replaced with
      `build_integral_derivatives(const IntegralPack&)` that reuses
      the SCF's resolved shell list. STO-3G legacy wrapper preserved
      as a thin shim. Test: `test_cpu_scf_def2svp_gradient.cpp`
      (translational invariance holds on H2O @ def2-SVP RHF).
- [x] **3-atom SCF gradient: exchange symmetrization bug** —
      canonical-enumeration exchange factor used `P(i,k)·P(j,l)` alone,
      but the 8-fold ERI symmetry requires the symmetrized product
      `½·[P(i,k)·P(j,l) + P(j,k)·P(i,l)]` (the `i↔j` permutation of a
      canonical `(ij|kl)` keeps the integral value but remaps the
      density product). Fix applied to BOTH libcint and native 2e
      branches. H2O STO-3G RHF gradient now matches PySCF to ~5 mHa/Bohr
      (was off by factor-3 + sign). H2O 6-31G RHF gradient matches FD.
      Tests: `test_def2svp_gradient.cpp::H2O_STO3G_analytic_vs_PySCF_reference`
      and `::H2O_631G_analytic_vs_FD_control`.
- [x] **d-shell gradient path** — two root causes fixed in sequence:
      (a) `cart_exponents` in `integral_derivatives.cpp` silently
      returned (0,0,0) for every d-shell component (only s/p were
      handled), so dS/dT/dV for d-shells were all computed at the
      wrong angular momentum; replaced with `fill_cart_comps` decode;
      (b) `native::md_e_coeffs` / `nuclear_prim` / `eri_prim` sized
      their Hermite/Boys tables for ai,bj ≤ 2 only, overflowing in
      the derivative path that raises ai by +1 (d → f-like temp).
      Buffers now sized from `kMaxL` (cart.h) to support the +1
      raise. H2O / 6-31G* and H2O / def2-SVP analytic gradients now
      match FD to 1e-2 eV/Å tolerance. Tests:
      `test_def2svp_gradient.cpp::H2O_631Gs_dS_dT_dV_vs_FD`,
      `::H2O_631Gs_analytic_vs_FD_minimal_d`,
      `::H2O_def2svp_analytic_vs_FD`.
- [x] **FD function takes BasisSet parameter** — previously hardcoded
      to STO-3G, which silently gave wrong FD when an analytic gradient
      was tested at any other basis (exactly how the pre-existing
      3-atom bug stayed hidden).
- [x] **basis_eval d-shell (l=2) extension** — generalized value,
      gradient, Hessian, and fused variants in `basis_eval.h` via the
      polynomial prefactor `dx^nx dy^ny dz^nz` decoded from `sh.cart`
      through `fill_cart_comps`. `xc_grid_engine.h` now delegates to
      `eval_basis_val_grad` so the XC Phi/dPhi tables inherit l>=2
      automatically. Cube writers drop their `l>1` guards.
      Tests: `test_cpu_scf_basis_eval_dshell` — FD check of grad+hess
      for every d-cart AO, plus RHF/6-31G* water density integral
      recovering n_electrons on a 0.25-Bohr Cartesian grid.
- [x] **KS/d-shell direct SCF** — root cause: `native::eri` always
      dispatched to `eri_prim_md`, whose hand-unrolled E-coefficient
      tables and Boys-function buffer are sized for l ≤ 1 only (Ltot ≤
      4). Any (dd|..) or (..|dd) call overflowed the F_m stack array
      and aborted via the stack canary. Fix: `native::eri` now routes
      quartets with max-l ≥ 2 through `eri_prim` (the general-l path
      already used by the derivative-raising code), which sizes its F_m
      and E tables from `kMaxL` / `kMDEDim`. PBE/6-31G* on water now
      converges in 12 iters to -76.327 Ha; grid-integrated n_electrons
      recovers 5.0 ± 0.003 per spin on the SG-0-ish grid. Test:
      `test_cpu_scf_basis_eval_dshell::PBE_631Gs_WaterConvergesAndIntegratesNElectrons`.
- [x] **LBFGS line search** — `LBFGSStep::line_search(c1, contract,
      min_alpha)` enables Armijo backtracking across successive driver
      steps: on step k the search_direction phase compares
      `E_k - E_{k-1} <= c1 * alpha * g_{k-1}.d_{k-1}`, rolling back to
      `x_{k-1}` and shrinking `alpha` by `contract` on failure without
      advancing the L-BFGS history. Integrates cleanly with the
      "one-force-per-step" driver loop — no new component or phase
      binding. Off by default; tutorials and benchmarks opt in via
      `step.line_search()`. Test:
      `test_cpu_minimizer_lbfgs::H2RHFArmijoLineSearchReachesTightForce`
      drives H2/STO-3G to the PySCF bond length within 0.01 Å at a
      1 mHa/Bohr force tolerance.
- [x] Expanded geometry set to 12: `benchmarks/qc/molecules/` now ships
      h2o, nh3, ch4, ethanol, glycine, benzene, naphthalene, aspirin,
      caffeine (RDKit / UFF-embedded) plus the open-shell trio
      (o2_triplet, no_doublet, ch3_radical). Regenerate via
      `python benchmarks/qc/molecules/_generate.py`.
- [x] Regenerated references at PBE/def2-SVP / cart=True for ten
      molecules (h2o, nh3, ch4, ethanol, glycine, benzene, naphthalene,
      o2_triplet, no_doublet, ch3_radical). Aspirin and caffeine have
      checked-in xyz but are deferred to the nightly CI path because the
      geomeTRIC optimization takes ~10+ min per molecule on a laptop.
      `generate_references.py` now supports open-shell systems
      (`spin` key in the JSON) and routes automatically to UHF/UKS.
- [x] CI label `benchmark-qc-fast` —
      `.github/workflows/benchmark-qc-fast.yml` runs the H2O/NH3/CH4
      triplet on push and the full-refs set on the nightly schedule.
      `bench_qc_suite` itself was rewritten to parse the JSON references
      at runtime (no embedded geometry lists) and hands off to
      `SCFEvaluator` + `LBFGSStep::line_search()` at the level of
      theory encoded in each reference.
- [x] **RKS energy/density convention bug in `scf_solve`** — discovered
      while driving the bench: the restricted-KS branch reused the HF
      trace formula `0.5·Tr[(H+F)P]` (which silently mixes Vxc into the
      energy for KS) and the UKS-style `P_total = P_alpha + P_beta` for
      J/XC feeding (which double-counts when the restricted convention
      stores `P_alpha` as the full density with `P_beta` aliased). Split
      the branch so RKS uses `Tr[H·P_α] + 0.5·Tr[J·P_α] + E_xc + V_nuc`
      and the Fock builder receives `P_total_J = P_α` plus spin-halved
      `(P_α/2, P_α/2)` for the XC engine. Verified by new test
      `BasisEvalDShell::PBE_def2SVP_RKS_Water_MatchesPyscf`.
- [ ] **KS/d-shell analytic-gradient accuracy** — the bench surfaces
      max-force ~1 eV/Bohr at the PySCF-optimized water geometry, where
      PySCF reports ~1e-4 Ha/Bohr. SCF energy matches PySCF start to
      ~80 μHa on the same geometry, so the issue is in the KS gradient
      assembly for d-shells (not SCF). Armijo line search rescues
      convergence but takes many backtracking iterations. Follow-up:
      validate `gradient.cpp` XC piece against FD at 6-31G* / def2-SVP
      and add a d-shell KS FD test alongside `test_def2svp_gradient`.
- [ ] **Denser SG-1/SG-2 grid** — the bench leaves a ~5 mHa residue
      between ATV and PySCF on def2-SVP because the SG-0-ish (50, 110)
      grid is coarser than PySCF's default. Tightening the radial/
      angular scheme is the prerequisite for the spec's 1 μHa
      acceptance; `bench_qc_suite` threshold is 1e-2 Ha until then.

### 1.3 Standard outputs (components, no workflow classes)
- [x] `cpu::minimizer::TrajectoryWriter` — molrs XYZ appender, stride support
- [x] `cpu::minimizer::OptLog` — TSV: step, energy_eV, max_force, rms_force,
      step_norm (tracks previous positions internally)
- [x] `cpu::minimizer::FinalGeometry` — single-frame XYZ, overwritten each
      end_of_step; after `run()` the file holds the last geometry visited.
      (Cross-component reference to `LBFGSStep::converged()` is impractical
      because the driver owns components by value; the overwrite model
      delivers the same observable outcome without the coupling.)
- [x] `cpu::scf::OrbitalInfoObserver` — HOMO/LUMO/gap TSV from `SCFResult`
- [x] `cpu::scf::DensityCubeWriter` — Gaussian `.cube`, respects RHF vs
      UHF/UKS density convention (`P_alpha` is already double-occupied in
      restricted, summed with `P_beta` only for unrestricted). s/p only.
- [x] `cpu::scf::OrbitalCubeWriter` — MO `.cube` by index + spin channel
- [x] **SCFEvaluator state exposure**: evaluator now caches the latest
      `SCFResult`/`IntegralPack`/`SCFConfig` in `shared_ptr` members and
      publishes a `SCFResultHandle` bundle. Observer components capture the
      bundle once and keep reading from it across the driver's value-copy
      of the evaluator.
- [x] Tests: `test_cpu_minimizer_recorders`,
      `test_cpu_scf_orbital_info_observer`, `test_cpu_recorder_density_cube`
      (all 22 suite tests green).

### 1.4 Locked default level of theory
- [x] `cpu::scf::default_config(integrals)` — restricted KS / PBE with
      MINAO guess + SG-0-ish grid (n_radial=50, n_angular=110).
      PBE/def2-SVP SCF itself still blocked on the d-shell XC-grid
      basis_eval extension tracked in 1.2; default_config runs today
      on STO-3G and 6-31G.
- [x] `cpu::scf::baseline_config(integrals)` — restricted HF / MINAO,
      tighter iteration budget, zero damping. Speed-focused regression
      preset for STO-3G fixtures.
- [x] `build_integrals_sto3g` wrapper + call sites removed across
      src/, tests/, examples/, benchmarks/. All consumers now call
      `build_integrals(frame, cpu::scf::BasisSet::STO3G, ...)` directly.
- [x] STO-3G migrated to the generic `basis_data::sto3g()` path (emitted
      by `tools/gen_basis_data.py`); `sto3g_data.h` and the STO-3G-only
      `build_shell_list_sto3g` branch are both deleted. libcint bridge,
      libcint reference header, and bench_libcint all re-plumbed onto
      the generic table.
- [x] `examples/tutorial/05_geom_opt.cpp` switched to `default_config`.
      Method-specific demos (02/03/04) keep their `rhf_solve`/`uhf_config`/
      `uks_config` plumbing — those tutorials *are* demonstrating the
      non-default level.

## Summary

Close the loop on everything that already works: extend basis-set support (prerequisite for any meaningful default), lock in a production level of theory, add a 10-molecule benchmark harness against PySCF, and wire **component-level** recorders/observers for the minimizer driver so optimization runs produce inspectable outputs. No new algorithms and no god-classes.

## Binding rules (apply here)

- All new "outputs" are `cpu::ComponentBase` subclasses bound to a `minimizer::Phase` or `scf::Phase`. No `Workflow`, no `Runner`, no `Report` facade classes.
- Basis-set data flows through plain value types (`BasisSpec`, `ShellBlockBasis`). No runtime string lookups inside hot loops.
- PySCF is ground truth for integrals, SCF energies, orbital energies, and geometries.

## 1.1 Basis-set foundation

`build_integrals_sto3g` hardcodes STO-3G. Required work:

**New API (replaces STO-3G-only entry point):**
```cpp
namespace cpu::scf {
  /// @brief Identifier for a basis set (compile-time enumerated).
  enum class BasisSet { STO3G, Pople631G, Pople631Gs, Def2SVP, Def2TZVP };

  /// @brief Immutable descriptor resolved from a BasisSet at build time.
  struct BasisSpec {
    std::string_view name;
    std::span<const ShellPrimitive> primitives;
    std::span<const ShellContraction> contractions;
  };

  [[nodiscard]] BasisSpec resolve_basis(BasisSet id);

  [[nodiscard]] IntegralPack build_integrals(
      const Frame& frame, BasisSet basis, int charge = 0);
}
```

- Data generator `tools/gen_basis_data.py` emits headers under `include/atomiverse/cpu/scf/basis_data/` (same pattern as MINAO). One header per basis set.
- `build_integrals_sto3g` becomes a thin inline wrapper `{ return build_integrals(frame, BasisSet::STO3G, charge); }` — **kept for backward compatibility only until Phase 1.4**, then removed.
- Rys ERI engine extended to l=2 (d-shells) — def2-SVP and 6-31G* both need d-polarization. FD self-consistency test on `(d,s|s,s)`, `(d,p|s,s)`, `(d,d|d,d)` quartets and derivatives.

### Acceptance (1.1)
- **Correctness:** max `|E_atv - E_pyscf|` ≤ 1e-7 Ha for H2O/CH4/benzene/aspirin × {STO3G, 6-31G, 6-31G*, def2-SVP, def2-TZVP}, RHF and UKS-PBE.
- **Correctness:** max `|∂E/∂R|` divergence from libcint ≤ 1e-5 Ha/Bohr on the same set.
- **Performance:** aspirin RHF/def2-SVP single SCF iteration ≤ 1.2× PySCF wall time.
- New test: `tests/cpu/scf/test_basis_vs_pyscf.cpp` with checked-in PySCF reference JSON.
- New test: `tests/cpu/scf/test_dshell_rys.cpp` — FD self-consistency for all d-containing quartets.

## 1.2 Benchmark suite (component-driven, not ad-hoc)

**Not a Python script orchestrating C++ internals.** Instead: a C++ benchmark binary that constructs drivers + components exactly as a user would, and emits machine-readable results.

**Molecule set** (checked in as `benchmarks/qc/molecules/*.xyz`):

| Category | Molecules |
|---|---|
| Rigid small | H2O, NH3, CH4 |
| Conjugated | benzene, naphthalene |
| Flexible | ethanol, glycine, aspirin, caffeine |
| Open-shell | O2 (triplet), NO (doublet), CH3• |

**Deliverables:**
- `benchmarks/qc/bench_small_molecule_suite.cpp` — loops molecules, builds `cpu::scf::SinglePoint` driver (Phase 2 dep — use `scf_solve` free-function in Phase 1 interim) + `cpu::minimizer::LBFGS`, records timings and energies into JSON.
- `benchmarks/qc/reference/*.json` — PySCF reference energies/geometries, produced by a checked-in `benchmarks/qc/reference/generate_references.py`.
- CI label `benchmark-qc-fast`: runs H2O+NH3+CH4 on every push; full suite nightly.

### Acceptance (1.2)
- **All 12 molecules converge** at PBE/def2-SVP under LBFGS with ≤ 1.5× PySCF step count.
- **Final geometry** RMS deviation from PySCF-optimized ≤ 1e-3 Å per molecule.
- **Final energy** within 1e-6 Ha of PySCF single-point on PySCF-optimized geometry.
- **Wall time** per full opt ≤ 1.5× PySCF on reference build (AMD/Intel, 10 threads).

## 1.3 Standard outputs (pure components, no workflow classes)

Every item here is a `ComponentBase` that attaches to an existing driver phase. Users compose with `opt.on(phase).add(component)` — no new driver types.

**New components (header + impl, Doxygen required):**

| Component | Driver | Bound phase | Role |
|---|---|---|---|
| `cpu::minimizer::TrajectoryWriter` | `LBFGS` | `minimizer::Phase::end_of_step` | Appends current geometry to `*.xyz` |
| `cpu::minimizer::OptLog` | `LBFGS` | `minimizer::Phase::end_of_step` | TSV: step, E[eV], max_F, rms_F, step_norm |
| `cpu::minimizer::FinalGeometry` | `LBFGS` | `minimizer::Phase::end_of_step` (runs only on converged flag) | Writes `optimized.xyz` |
| `cpu::scf::OrbitalInfoObserver` | any SCF driver | `scf::Phase::end_of_step` | Pulls HOMO/LUMO/gap from `SCFResult`, logs them |
| `cpu::scf::DensityCubeWriter` | any SCF driver | `scf::Phase::end_of_step` | Gaussian `.cube` on user-supplied grid |
| `cpu::scf::OrbitalCubeWriter` | any SCF driver | `scf::Phase::end_of_step` | Same for a chosen MO index |

**Files (new):**
- `include/atomiverse/cpu/minimizer/trajectory_writer.h`
- `include/atomiverse/cpu/minimizer/opt_log.h`
- `include/atomiverse/cpu/minimizer/final_geometry.h`
- `include/atomiverse/cpu/scf/orbital_info_observer.h`
- `include/atomiverse/cpu/recorder/density_cube_writer.h`
- `include/atomiverse/cpu/recorder/orbital_cube_writer.h`

Composition example (this is how users consume it):
```cpp
cpu::minimizer::LBFGS opt;
opt.add(cpu::minimizer::LBFGSStep{/*...*/});
opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
opt.on(minimizer::Phase::end_of_step)
   .add(cpu::minimizer::TrajectoryWriter{"aspirin_traj.xyz", /*stride=*/1})
   .add(cpu::minimizer::OptLog{"aspirin_opt.tsv"})
   .add(cpu::minimizer::FinalGeometry{"aspirin_final.xyz"});
opt.build();
opt.run(frame, 200);
```

### Acceptance (1.3)
- **Correctness:** HOMO/LUMO match PySCF to ≤ 1e-4 eV on all 12 benchmark molecules.
- **Correctness:** density cube integral ≈ N_electrons to ≤ 1e-3 absolute (grid-resolution dependent; fix grid in test).
- **Correctness:** trajectory `.xyz` loads in ASE without error; atom count + symbols match input.
- **Performance:** recorder overhead ≤ 1% of total opt wall time (verified with/without recorders).
- New test: `tests/cpu/minimizer/test_trajectory_writer.cpp`, `tests/cpu/scf/test_orbital_info_observer.cpp`, `tests/cpu/recorder/test_density_cube.cpp`.

## 1.4 Locked default level of theory

**Free function, not a class:**
```cpp
namespace cpu::scf {
  /// @brief Production default (PBE/def2-SVP, RKS closed-shell, dual convergence).
  [[nodiscard]] SCFConfig default_config(const IntegralPack& integrals);

  /// @brief Regression baseline (RHF/STO-3G, speed-focused).
  [[nodiscard]] SCFConfig baseline_config(const IntegralPack& integrals);
}
```

- Remove the interim `build_integrals_sto3g` wrapper — all examples migrate to `build_integrals(frame, BasisSet::Def2SVP)` or `default_config`.
- Every tutorial (`examples/tutorial/0*.cpp`) updates to the default unless explicitly demonstrating another.

### Acceptance (1.4)
- `default_config` + benchmark suite reproduces PySCF PBE/def2-SVP to 1e-6 Ha on all 12 molecules.
- No call site in the repo references `build_integrals_sto3g` outside deprecation shim tests.

## Files to create / modify

| File | Action |
|---|---|
| `tools/gen_basis_data.py` | new |
| `include/atomiverse/cpu/scf/basis_data/*.h` | new (per basis) |
| `include/atomiverse/cpu/scf/integrals.h` | add `BasisSet`, `BasisSpec`, `build_integrals` |
| `src/cpu/scf/integrals.cpp` | dispatch by `BasisSet` |
| `src/cpu/scf/native_integrals.cpp` | extend to l=2 |
| `src/cpu/scf/rys_eri.cpp` | audit `kRysMaxRoots` for d-quartet derivatives |
| `include/atomiverse/cpu/minimizer/{trajectory_writer,opt_log,final_geometry}.h` | new |
| `include/atomiverse/cpu/scf/orbital_info_observer.h` | new |
| `include/atomiverse/cpu/recorder/{density_cube_writer,orbital_cube_writer}.h` | new |
| `src/cpu/**/*.cpp` | impls for above |
| `benchmarks/qc/bench_small_molecule_suite.cpp` | new |
| `benchmarks/qc/molecules/*.xyz` | new (12 files) |
| `benchmarks/qc/reference/generate_references.py` | new |
| `benchmarks/qc/reference/*.json` | new (checked-in PySCF refs) |
| `tests/cpu/scf/test_basis_vs_pyscf.cpp` | new |
| `tests/cpu/scf/test_dshell_rys.cpp` | new |
| `tests/cpu/minimizer/test_trajectory_writer.cpp` | new |
| `tests/cpu/scf/test_orbital_info_observer.cpp` | new |
| `tests/cpu/recorder/test_density_cube.cpp` | new |
| `examples/tutorial/0*.cpp` | migrate to defaults + recorders |

## Dependencies

None upstream (this is the foundation). Downstream: Phases 2–5 all require `BasisSet` enum and default_config to exist.

## Risks

- **d-shell Rys extension** is the highest-risk line item. Gate: FD self-consistency test on all d-quartets must pass before any other 1.x work.
- **PySCF reference generation** must be version-pinned (`pyscf==2.x.y`) and checked in — re-running it on a new PySCF minor release will cause spurious test regressions.

## Exit gate

All of 1.1–1.4 acceptance criteria green, 12-molecule benchmark suite in CI for 7 consecutive days without failure, user sign-off.

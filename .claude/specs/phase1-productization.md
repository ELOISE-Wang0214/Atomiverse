---
title: Phase 1 — Productize the CPU SCF + gradient + opt stack
status: in progress
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
- [ ] **KS/d-shell XC grid** — `basis_eval.h` only handles l ∈ {0, 1}.
      Blocks PBE/def2-SVP SCF itself (KS grid evaluation). Separate
      from the gradient issue; surfaces as SCF hang before any force
      is computed.
- [ ] **LBFGS line search** — current fixed-step LBFGSStep converges
      to ~10–30 mHa drift from the PySCF minimum because dtol triggers
      on small-step iterations that haven't actually minimized. Full
      1.2 acceptance (1e-6 Ha energy parity) needs a line search.
- [ ] Expand geometry set to 12 (benzene, naphthalene, ethanol,
      glycine, aspirin, caffeine, O2 triplet, NO doublet, CH3•)
- [ ] Regenerate references for the full set at PBE/def2-SVP (d-shell
      gradient unblocked; still needs KS/d-shell XC grid below)
- [ ] CI label `benchmark-qc-fast` (H2O+NH3+CH4 per push) + nightly
      full suite

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
- [ ] `cpu::scf::default_config(integrals)` — PBE/def2-SVP production
- [ ] `cpu::scf::baseline_config(integrals)` — RHF/STO-3G speed path
- [ ] Remove `build_integrals_sto3g` from all non-deprecation call
      sites (migrate `integral_derivatives.cpp`, tutorials, tests)
- [ ] Migrate STO-3G itself to the generic `basis_data::` path, drop
      `sto3g_data.h`
- [ ] All `examples/tutorial/0*.cpp` use `default_config` unless
      demonstrating another level

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

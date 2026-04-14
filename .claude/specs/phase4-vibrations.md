---
title: Phase 4 — Vibrational analysis (Hessian driver + thermochemistry)
status: draft
created: 2026-04-14
---

# Phase 4 — Vibrations

## Summary

Introduce a new driver type `cpu::vib::Hessian` with its own phase enum `ATV::vib::Phase`, composed of small single-responsibility components: geometry displacement, gradient evaluation (reuses `SCFEvaluator`), Hessian accumulation, mass-weighting / diagonalization, thermochemistry. Numerical-Hessian first; analytic Hessian deferred to a future spec.

## Binding rules (apply here)

- **New driver type** because vibrations have a distinct phase sequence that neither MD, minimizer, nor SCF captures (displace → force × 6N → assemble → diagonalize → thermochem).
- All physics contributions are components bound to phases of this new driver. No `HessianCalculator` god-class.
- Force evaluation reuses the existing `SCFEvaluator` component (bound to `vib::Phase::force`), unchanged — the whole point of the driver/component model.
- PySCF `pyscf.hessian.RHF(mf).kernel()` + `pyscf.hessian.thermo.harmonic_analysis(mol, hess)` is ground truth for frequencies, ZPE, and thermodynamic corrections.

## New phase enum

```cpp
// include/atomiverse/facade/stage.h
namespace ATV::vib {

enum class Phase {
    setup,              ///< Allocate Hessian buffer, cache reference geometry
    displace,           ///< Advance to next ±h displacement
    force,              ///< Evaluate gradient at displaced geometry (SCFEvaluator reused)
    accumulate,         ///< Add ∓g / (2h) row into Hessian buffer
    diagonalize,        ///< Mass-weight + symmetrize + dsyevd
    thermochem,         ///< ZPE, thermal correction, Gibbs free energy
    end_of_step,        ///< Observers / recorders
    COUNT
};

inline std::vector<std::string> phase_names(Phase);

} // namespace
```

The driver's `run()` loop iterates **displacement index** (0 .. 6N−1) rather than "time step". On each iteration, phases `displace`, `force`, `accumulate` execute; after the final iteration, `diagonalize` and `thermochem` execute once.

## Driver design

```cpp
namespace ATV::cpu::vib {

class Hessian : public Driver<ATV::vib::Phase, cpu::device_tag> {
public:
    struct Config {
        double displacement = 1e-3;   ///< Finite-difference step h [Å]
        double temperature  = 298.15; ///< For thermochemistry [K]
        double pressure     = 101325; ///< For Gibbs free energy [Pa]
        bool   use_central_difference = true;  ///< else forward FD
    };

    explicit Hessian(Config config);

private:
    void default_bindings() override;  // DisplaceStep, HessianAccumulator,
                                       // MassWeighting, Diagonalize, Thermochem
    Config config_;
};

} // namespace
```

## Components (all new, single responsibility)

| Component | Phase | Role |
|---|---|---|
| `cpu::vib::ReferenceGeometry` | `setup` | Cache initial `frame.atoms.x/y/z` for later restoration |
| `cpu::vib::DisplaceStep` | `displace` | Move atom `i`, axis `d` by `±h` per FD index; restore on exit |
| `cpu::vib::HessianAccumulator` | `accumulate` | `H[i,d, j,e] += g(j,e) / (2h)` (or `(g_plus - g_minus)/(2h)` for central) |
| `cpu::vib::Symmetrize` | `diagonalize` (order -1) | `H = (H + H^T) / 2` |
| `cpu::vib::MassWeighting` | `diagonalize` | `H_mw = M^(-1/2) H M^(-1/2)` |
| `cpu::vib::EigenSolver` | `diagonalize` (order +1) | `dsyevd` → eigenvalues (√λ → ω), eigenvectors → normal modes |
| `cpu::vib::Thermochem` | `thermochem` | ZPE, U(T), H(T), S(T), G(T) from harmonic-oscillator + rigid-rotor + translation partition functions |
| `cpu::vib::FrequencyReport` | `end_of_step` | Logs frequencies [cm⁻¹], imaginary count, ZPE, G(T) |
| `cpu::vib::NormalModeWriter` | `end_of_step` | Optional XYZ animation of each mode for visualization |

`cpu::scf::SCFEvaluator` is bound to `vib::Phase::force` by the user exactly like MD/minimizer:

```cpp
cpu::vib::Hessian hess{{/*displacement=*/1e-3, /*T=*/298.15}};
hess.on(vib::Phase::force).add(
    cpu::scf::SCFEvaluator{cpu::scf::b3lyp_config(integrals)});
hess.on(vib::Phase::end_of_step).add(
    cpu::vib::FrequencyReport{"frequencies.tsv"});
hess.build();
hess.run(frame, /*n_displacements=*/6 * n_atoms);
```

Mixed workflow (opt + freq) is **two drivers in user code**:
```cpp
// 1) Optimize geometry
cpu::minimizer::LBFGS opt;
opt.add(cpu::minimizer::LBFGSStep{/*...*/});
opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
opt.build(); opt.run(frame, 200);

// 2) Frequencies at the optimized geometry
cpu::vib::Hessian hess{{1e-3, 298.15}};
hess.on(vib::Phase::force).add(cpu::scf::SCFEvaluator{config});
hess.build(); hess.run(frame, 6 * frame.atoms.size());
```

## Frame additions

New `frame.hessian` sub-object, allocated in `setup`:
```cpp
namespace ATV::cpu {
struct FrameHessian {
    xt::xtensor<double, 2> H;         ///< (3N, 3N), Hartree/Bohr²
    xt::xtensor<double, 1> omega;     ///< 3N frequencies [cm⁻¹]
    xt::xtensor<double, 2> modes;     ///< (3N, 3N) eigenvectors
    double zpe = 0.0;                 ///< eV
    double g_thermal = 0.0;           ///< eV (G(T) - E_elec)
    int    n_imaginary = 0;
};
}
```

Attached to `Frame` as `std::optional<FrameHessian> hess`; allocated only when a Hessian driver runs.

## Acceptance (numerical Hessian only)

- **Correctness:** frequencies match PySCF `hessian.thermo.harmonic_analysis` to ≤ 5 cm⁻¹ on {H2O, NH3, CH4, benzene, formaldehyde}.
- **Correctness:** ZPE matches PySCF to ≤ 0.1 kJ/mol (≈ 0.001 eV).
- **Correctness:** G(298.15 K, 1 atm) matches PySCF to ≤ 0.5 kJ/mol.
- **Correctness:** translational (3) + rotational (3 for nonlinear, 2 for linear) modes have |ω| < 10 cm⁻¹ at an optimized geometry.
- **Correctness:** transition-state check — NH3 inversion saddle produces exactly one imaginary frequency.
- **Performance:** H2O frequency calculation ≤ 1.2× PySCF `hessian.RHF`; aspirin ≤ 2× PySCF (PySCF uses analytic Hessian, we use numerical — 2× is a reasonable cap for a numerical implementation).
- **Wall-time parity:** aspirin numerical Hessian at PBE/def2-SVP ≤ 10 min on reference 10-thread build (126 gradient evaluations × ~3 s each expected).
- Test: `tests/cpu/vib/test_hessian_vs_pyscf.cpp`.

## Analytic Hessian — deferred

Analytic Hessian requires CPHF solver + second derivatives of ERIs. Scope-separate; will receive its own `phase4b-analytic-hessian.md` spec only if numerical Hessian proves prohibitively slow for routine use. Not part of this phase.

## Files to create

| File | Action |
|---|---|
| `include/atomiverse/facade/stage.h` | add `vib::Phase` enum |
| `include/atomiverse/cpu/vib/hessian.h` | new driver |
| `include/atomiverse/cpu/vib/components/*.h` | new (DisplaceStep, HessianAccumulator, MassWeighting, etc.) |
| `src/cpu/vib/hessian.cpp` | new |
| `src/cpu/vib/components/*.cpp` | new |
| `include/atomiverse/cpu/frame.h` | optional `hess` field |
| `tests/cpu/vib/test_hessian_vs_pyscf.cpp` | new |
| `tests/cpu/vib/test_displace_step.cpp` | new |
| `tests/cpu/vib/test_thermochem_vs_pyscf.cpp` | new |
| `benchmarks/qc/reference/generate_references.py` | add PySCF hessian refs |
| `examples/tutorial/08_opt_and_freq.cpp` | new (composition example) |

## Dependencies

- **Phase 3 required** so frequencies can be requested at B3LYP/6-31G* (the stated target level).
- **Phase 1 required** for basis-set support.
- **Phase 2 required** because `SCFEvaluator` internally uses `SinglePoint`.

## Risks

- **Finite-difference noise:** h=1e-3 Å may introduce ~1 cm⁻¹ error in low modes. Mitigation: central difference, h tunable in config, regression test pinned to PySCF with 5 cm⁻¹ slack.
- **Symmetry / rotational-frame contamination:** at non-stationary geometries, FD Hessian is asymmetric. Mitigation: symmetrize before diagonalization; document that Hessian must be computed at a converged optimum.
- **Performance:** 6N gradients. For caffeine (24 atoms, 144 gradients × ~5 s each = 12 min). Document as expected; users aware that frequency analysis is expensive.

## Exit gate

All acceptance criteria green + opt+freq example produces physically sensible output on aspirin (no imaginary modes at the optimum) + user sign-off.

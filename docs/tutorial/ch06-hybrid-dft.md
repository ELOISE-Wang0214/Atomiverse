# Chapter 6: Hybrid Functionals and Analytic Gradients

| | |
|---|---|
| **Examples** | [`04_uks.cpp`](https://github.com/MolCrafts/Atomiverse/blob/dev/examples/04_uks.cpp) (LSDA baseline), [`05_geom_opt.cpp`](https://github.com/MolCrafts/Atomiverse/blob/dev/examples/05_geom_opt.cpp) (PBE forces + L-BFGS) |
| **API** | `cpu::scf::uks_config` / hybrid `XCType::{B3LYP,PBE0}`, `scf_gradient` |

Hybrids mix a fraction of exact Hartree–Fock exchange into the KS functional.
Analytic nuclear gradients (`scf_gradient`, or `SCF{}.with_gradient(true)`) are
what turn a single-point energy into geometry optimization and AIMD
([Chapter 12](ch12-geometry-optimization.md), [Chapter 15](ch15-cross-method.md)).

```bash
# UKS-LSDA scan (energy components) — results in examples/results/04_uks.json
./build/examples/04_uks
# Geometry opt with DFT forces (aspirin) — longer run
./build/examples/05_geom_opt
```

Validate gradients with finite differences (policy from
[Chapter 1](ch01-foundations.md)): analytic and FD forces must agree within a
documented tolerance before trusting an optimizer step.


[Chapter 5](ch05-kohn-sham-dft.md) introduced Kohn-Sham DFT with the pure functionals
LSDA and PBE, where the exchange-correlation energy is a semi-local integral
over the density on a grid. This chapter adds the two most widely used
**global hybrid** functionals — B3LYP and PBE0 — which mix a fixed fraction
of *exact* Hartree-Fock exchange into that semi-local recipe, and then shows
how to get **analytic nuclear gradients** (forces) for them. Hybrids are the
workhorse of molecular DFT precisely because that exchange admixture
corrects the self-interaction error that plagues pure functionals; analytic
gradients are what make geometry optimization and vibrational analysis on top
of them affordable.

The six questions:

1. **The concept.** Why mix exact exchange in, and what "global hybrid"
   means.
2. **The Atomiverse representation.** `XCType::B3LYP` / `PBE0`, the
   exact-exchange fraction, and how the Fock builder splits the work.
3. **The C++ API.** `RKS` with `XCType::B3LYP` / `PBE0`, and
   `scf_gradient`.
4. **The expected result.** A converged hybrid energy and forces that match
   PySCF's `grad.RKS`.
5. **The validation.** Analytic force vs PySCF *and* vs an independent
   finite-difference oracle, to 1e-5 Ha/Bohr.
6. **The recurring pattern.** The gradient plugs into the same
   `SCFEvaluator` that drives AIMD and geometry optimization.

This chapter assumes you are comfortable with the SCF loop and the
`SCFConfig` / `SCFResult` types from Tutorials 2 and 4.

## The Concept: A Fraction of Exact Exchange

Pure semi-local functionals (LDA, GGA) describe exchange approximately, from
the local density and its gradient. That approximation carries a
**self-interaction error**: an electron feels a spurious repulsion from its
own charge density, which pure functionals do not fully cancel. Hartree-Fock
exchange, being computed from the orbitals themselves, cancels
self-interaction *exactly* — but HF has its own failing (no correlation, and
it overshoots exchange for delocalized systems).

A **global hybrid** takes the pragmatic middle: replace a fixed fraction
$\alpha_{\text{HF}}$ of the functional's semi-local exchange with exact HF
exchange, and keep the rest semi-local. The KS exchange-correlation energy
becomes

$$
E_{\text{xc}}^{\text{hybrid}}
  = \alpha_{\text{HF}}\, E_{\text{x}}^{\text{HF}}
  + (1 - \alpha_{\text{HF}})\, E_{\text{x}}^{\text{semi-local}}
  + E_{\text{c}}^{\text{semi-local}}.
$$

The two functionals in this chapter fix $\alpha_{\text{HF}}$ at values
established by their authors:

- **B3LYP** — $\alpha_{\text{HF}} = 0.20$, over the mix
  `0.08 Slater-X + 0.72 B88-X + 0.19 VWN_RPA-C + 0.81 LYP-C`
  (libxc `XC_HYB_GGA_XC_B3LYP`).
- **PBE0** — $\alpha_{\text{HF}} = 0.25$, over `0.75 PBE-X + 1.00 PBE-C`
  (libxc `XC_HYB_GGA_XC_PBEH`).

The word "global" means the fraction is a single constant everywhere in
space (as opposed to range-separated hybrids, where it varies with
interelectronic distance — out of scope here).

## The Representation: The Fock Builder Splits the Work

Atomiverse selects a functional with the `XCType` enum
([`xc_functional.h`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/scf/xc_functional.h)):

```cpp
enum class XCType { LSDA, PBE, B3LYP, PBE0 };
```

`LSDA` and `PBE` are pure; `B3LYP` and `PBE0` are hybrids. The split between
semi-local and exact exchange is handled cleanly: the grid engine's
`evaluate_xc` returns only the **semi-local** part, and the Fock builder adds
the HF-exchange admixture separately, weighted by the functional's fraction:

```cpp
double a = cpu::scf::hybrid_exchange_fraction(cpu::scf::XCType::B3LYP); // 0.20
bool   g = cpu::scf::is_gga(cpu::scf::XCType::B3LYP);                   // true
```

`hybrid_exchange_fraction` returns `0.0` for pure functionals — the signal
"no exchange-matrix build needed" — and the nonzero fraction for hybrids.
`is_gga` is `true` for PBE, B3LYP, and PBE0 (all carry a density-gradient
term), which tells the grid engine it must build basis-function gradients.
The practical consequence is cost: a hybrid pays for *both* a KS grid *and* a
Hartree-Fock exchange (K) matrix build every SCF iteration, so it is
meaningfully more expensive than a pure GGA. That expense is the price of the
self-interaction correction.

## The C++ API: Two Config Builders, Then Gradients

The convenience builders in
[`scf_types.h`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/scf/scf_types.h) preset a
production-grade closed-shell hybrid configuration — restricted KS, the
PySCF-matching `GridLevel::MatchPySCF` integration grid, MINAO initial guess,
`tol = 1e-9`, DIIS on — and differ only in the functional:

```cpp
using namespace ATV;

cpu::Frame frame = /* H2O, positions in Å, type = Z */;

/* integrals private inside RHF/UHF/RKS/UKS */
cpu::scf::RKS rks(cpu::scf::XCType::B3LYP);  // or XCType::PBE0
rks.basis = cpu::scf::BasisSet::Pople631Gs;
cfg.max_iter = 300;

const auto& res = rks.run(frame);

const auto& result = rhf.run(frame);                    // SCFResult
// result.energy (Hartree), result.converged, result.xc_energy, ...
```

`b3lyp_config` and `pbe0_config` both clone `default_config` and set only the
`xc_type`, so everything else — grid, guess, convergence — is the same locked
production default. The hybrids use a GGA grid, so a polarized basis such as
6-31G* (`Pople631Gs`) or def2-SVP is the natural choice; see
[Chapter 4](ch04-basis-sets.md).

Forces come from a *separate* call on the converged result. From
[`gradient.h`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/scf/gradient.h):

```cpp
cpu::scf::GradientResult grad =
    cpu::scf::scf_gradient(frame, integrals, result, cfg);

// grad.forces : xtensor [n_atoms, 3] in eV/Å   (FORCE = -dE/dR)
// grad.energy : total energy in eV
```

`scf_gradient` is the *unified* analytic gradient for every mean-field
method — it dispatches on the config, computing the HF Coulomb+exchange
two-electron terms for HF configs and Coulomb + XC-grid contributions for KS,
with the hybrid case adding the exchange gradient weighted by
$\alpha_{\text{HF}}$. The same function serves RHF, UHF, RKS, UKS, and both
hybrids; you never pick a method-specific gradient routine. Note the units
and sign: `scf_gradient` returns **forces** (the negative energy gradient) in
eV/Å, consistent with what the frame's `fx/fy/fz` arrays expect, whereas
PySCF's `grad` returns $\partial E/\partial R$ in Hartree/Bohr — convert with
`core/units.h` before comparing.

Putting the pieces together, a complete single-point-energy-and-forces
routine for a hybrid is:

```cpp
auto run_hybrid = [](cpu::Frame& frame, cpu::scf::XCType xc) {
    /* integrals private inside RHF/UHF/RKS/UKS */
    cpu::scf::SCFConfig cfg = (xc == cpu::scf::XCType::B3LYP)
                                  ? cpu::scf::XCType::B3LYP
                                  : cpu::scf::XCType::PBE0;
    cpu::scf::RKS rks(xc);
    cfg.max_iter = 300;
    const auto& res = rks.run(frame);
    solver.build(ints);
    const auto& result = solver.run(ints);
    return cpu::scf::scf_gradient(frame, ints, result, cfg);
};
```

## The Expected Result: Converged Energy, PySCF-Grade Forces

For a closed-shell molecule at a sensible geometry, a hybrid SCF should
converge in a modest number of DIIS iterations (the MINAO guess and
`GridLevel::MatchPySCF` grid are chosen for reliable convergence), and the
analytic forces should reproduce PySCF's `grad.RKS` component by component.
As a concrete anchor, B3LYP/6-31G* on water yields forces on the order of a
few times $10^{-3}$ Ha/Bohr with the correct symmetry (the two O–H forces
mirror across the C₂ axis), and the sum of forces over all atoms vanishes to
$\sim 10^{-8}$ Ha/Bohr — translational invariance, a free structural check
that no net force can push the isolated molecule.

The optimized geometry follows the forces: minimizing the B3LYP/6-31G*
energy of water lands at an O–H bond of ≈ 0.969 Å and an H–O–H angle of
≈ 103.7°, close to experiment and matching a PySCF-optimized reference to
sub-milliångström RMS. This is the pay-off of the self-interaction
correction — pure LDA over-binds and gives too-short bonds; the hybrid gets
the structure right.

## The Validation: Analytic vs PySCF *and* vs Finite Difference

Hybrid gradients are validated two independent ways, and the redundancy is
deliberate:

1. **Against PySCF (`grad.RKS`).** The Atomiverse analytic force is compared,
   component by component, against a genuine PySCF 2.x `nuc_grad_method`
   reference at the same functional / basis / grid, to **≤ 1e-5 Ha/Bohr** on
   H₂O, CH₄, and NH₃. This certifies the gradient against a trusted external
   engine.
2. **Against a moving-grid finite difference.** Because a shared bug could in
   principle affect *both* Atomiverse and PySCF's analytic derivations, the
   test *also* compares the analytic force against
   `scf_finite_difference_gradient` — central differences of the hybrid
   energy, with the XC grid regenerated at each displaced geometry — again to
   **≤ 1e-5 Ha/Bohr**. This oracle trusts nothing but the energy itself (the
   pattern from [Chapter 1](ch01-foundations.md)):

```cpp
const double h_ang = 1e-3 / ATV::units::kAngstromInBohr;   // 1e-3 Bohr in Å
auto fd = cpu::scf::scf_finite_difference_gradient(
    frame, cfg, /*charge=*/0, h_ang, cpu::scf::BasisSet::Pople631Gs);
// compare fd.forces(A, d) against scf_gradient(...).forces(A, d)
```

The moving-grid FD check matters for hybrids specifically: the XC grid moves
with the nuclei, so the gradient has a **grid-response** term (the
weight-derivative and grid-coordinate contributions) that a naïve fixed-grid
implementation would miss. Agreement between the analytic gradient and a
genuinely moving-grid finite difference is what proves that term is present
and correct.

Underneath, the semi-local part of each hybrid is itself pinned: the native
XC functionals match libxc to **1e-12**, so the only remaining freedom is the
exact-exchange admixture and its gradient, which the two checks above nail
down. Each result carries its class: the analytic-vs-FD agreement is a
**numerical-tolerance** claim (same quantity, two computations); the
optimized-geometry-vs-experiment gap is a **physical-approximation** claim
(the functional, not the code).

## The Recurring Pattern: The Gradient Is a Force Provider

Everything in this chapter produced a `GradientResult` by hand, which is the
right way to *validate* a gradient. In a real workflow you rarely call
`scf_gradient` directly — you let the **`SCFEvaluator`** call it for you.
`SCFEvaluator` wraps exactly this "build integrals → run SCF → compute
analytic gradient → write forces into the frame" sequence into a `force`-phase
component, and it is functional-agnostic: hand it a `b3lyp_config` and it
provides B3LYP forces to whichever driver it is plugged into.

```cpp
// B3LYP geometry optimization: hybrid forces into the Minimize driver
cpu::minimizer::Minimize opt;
opt.add(cpu::minimizer::LBFGS{/* memory, alpha0, max_ls */});
opt.add_handler(cpu::ForceConverged{/* the five tolerances */});
opt.on(minimizer::Phase::force)
   .add(cpu::scf::RKS(cpu::scf::XCType::B3LYP) /* basis via member */,
                               cpu::scf::BasisSet::Pople631Gs});
opt.build();
opt.run(frame, /*max_iter=*/100);
```

The same `RKS{XCType::B3LYP}` drops into `cpu::md::MD` for
hybrid-functional AIMD ([Chapter 15](ch15-cross-method.md)) and into `cpu::vib::Hessian`
for hybrid vibrational frequencies ([Chapter 14](ch14-vibrations.md)). One
gradient routine, one evaluator component, three drivers — the "swap the
evaluator, keep the driver" architecture that
[Chapter 15](ch15-cross-method.md) makes explicit. The functional you chose
in this chapter travels with the evaluator wherever forces are needed.

## Validation Reference

Hybrid analytic gradients are validated in
[`tests/cpu/scf/test_hybrid_gradient.cpp`](https://github.com/MolCrafts/Atomiverse/blob/dev/tests/cpu/scf/test_hybrid_gradient.cpp):
B3LYP and PBE0 forces vs PySCF `grad.RKS` (**≤ 1e-5 Ha/Bohr** on H₂O/CH₄/NH₃),
vs moving-grid finite difference (**≤ 1e-5 Ha/Bohr**, an oracle independent of
PySCF), translational invariance (`|Σ_A F_A| ≤ 1e-8 Ha/Bohr`), and a
B3LYP/6-31G* geometry-optimization endpoint vs a PySCF-optimized reference.
Hybrid single-point energies are pinned in
[`tests/cpu/scf/test_hybrid_singlepoint.cpp`](https://github.com/MolCrafts/Atomiverse/blob/dev/tests/cpu/scf/test_hybrid_singlepoint.cpp),
and the semi-local XC parts vs libxc (1e-12) in `test_xc_vs_libxc`.

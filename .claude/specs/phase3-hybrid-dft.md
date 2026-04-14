---
title: Phase 3 — Hybrid DFT (B3LYP, PBE0)
status: draft
created: 2026-04-14
---

# Phase 3 — Hybrid DFT

## Summary

Add hybrid functionals (B3LYP, PBE0) to the existing SCF + gradient components. No new drivers, no new wrappers — only extending the `HamiltonianBuilder` and gradient accumulator to scale the HF-exchange contribution by a functional-specific `alpha_hf`. Composition with `SinglePoint` and `LBFGS` drivers is unchanged from user perspective.

## Binding rules (apply here)

- **No new driver.** Hybrid DFT uses the same `SinglePoint` and `LBFGS` drivers composed with a different `SCFConfig`. Zero new classes with "Hybrid" in the name.
- Hybrid blending is a **parameter** on existing components, not a new component. A hybrid SCF run differs from a pure DFT run only by `alpha_hf != 0` in config.
- PySCF `pyscf.dft.RKS(mol, xc='b3lyp').kernel()` and `pyscf.grad.RKS(mf).kernel()` are ground truth.
- libxc is the ground-truth reference for XC energy density values; HF-exchange mixing is under our control and verified against PySCF end-to-end.

## Landing order (do not overlap)

Three sequential sub-stages, each gated independently:

### 3.1 Hybrid single-point

**Data changes:**
```cpp
namespace cpu::scf {

enum class XCFunctional {
    LDA_VWN5, PBE, B3LYP, PBE0   // extend existing XCType
};

struct XCSpec {
    XCFunctional functional;
    double alpha_hf;       ///< HF-exchange fraction (0 for pure, 0.20 for B3LYP, 0.25 for PBE0).
    // (additional blend coefficients for LYP/VWN only populated internally)
};

[[nodiscard]] XCSpec make_xc_spec(XCFunctional f);

// SCFConfig gains:
struct SCFConfig {
    // ... existing fields ...
    XCSpec xc{make_xc_spec(XCFunctional::PBE)};
};

[[nodiscard]] SCFConfig b3lyp_config(const IntegralPack& integrals);
[[nodiscard]] SCFConfig pbe0_config(const IntegralPack& integrals);
}
```

**Behavior:** `HamiltonianBuilder` (Phase 2 component) already builds `H = h_core + J + V_xc`. Extend it to include `-alpha_hf * K` when `alpha_hf > 0`. The K matrix is already computed in the RHF/UHF path — no new integral code.

**Ground-truth coefficients** — hard-coded from PySCF `pyscf/dft/libxc.py` XC-code tables, cross-checked against libxc headers:
- **B3LYP** (Becke 1993 VWN-III variant used by PySCF default): 0.20 HF + 0.80 (0.08 Slater + 0.72 B88) + (0.19 VWN3 + 0.81 LYP)
- **PBE0**: 0.25 HF + 0.75 PBE_X + 1.00 PBE_C

**Acceptance (3.1):**
- **Correctness:** B3LYP single-point energy matches PySCF `pyscf.dft.RKS(mol, xc='b3lyp').kernel()` to ≤ 1e-7 Ha on {H2O, CH4, NH3, benzene, glycine} × 6-31G*.
- **Correctness:** PBE0 energy matches PySCF `xc='pbe0'` to ≤ 1e-7 Ha on same set.
- **Correctness:** orbital energies (HOMO, LUMO) match PySCF to ≤ 1e-5 eV.
- **Performance:** B3LYP/6-31G* SCF wall time on aspirin ≤ 1.5× PySCF (hybrid is K-heavy; target is looser than pure-functional baseline).
- Test: `tests/cpu/scf/test_hybrid_scf_vs_pyscf.cpp`.

### 3.2 Hybrid gradient

The analytic gradient already splits `F = h + J - alpha*K + V_xc` into contributions per shell-pair derivative. Extend the 2e-gradient contraction from Phase 2's `integral_derivatives.cpp` so the exchange-channel multiplier uses `alpha_hf` from config rather than the hard-coded 1.0 (HF) / 0.0 (pure DFT) dichotomy.

**No change to `rys_eri_block_deriv`.** Only the density-weighted accumulation changes:
```cpp
const double factor_J = 0.5 * mult * P(i,j) * P(k,l);
const double factor_K = -0.25 * mult * config.xc.alpha_hf * P(i,k) * P(j,l);
```

**Acceptance (3.2):**
- **Correctness:** B3LYP force vs PySCF `pyscf.grad.RKS(mf).kernel()` max |Δ| ≤ 1e-5 Ha/Bohr on same 5 molecules.
- **Correctness:** PBE0 same tolerance.
- **Correctness:** force translational invariance `Σ F_atom = 0` to ≤ 1e-10 Ha/Bohr.
- **Performance:** B3LYP/6-31G* aspirin analytic-gradient wall time ≤ 2× PySCF (hybrid K-derivative is dominant, no density fitting planned).
- Test: `tests/cpu/scf/test_hybrid_gradient_vs_pyscf.cpp`.

### 3.3 Hybrid geometry optimization (composition-only)

**No code change** beyond handing a `b3lyp_config` to the existing LBFGS composition:

```cpp
cpu::minimizer::LBFGS opt;
opt.add(cpu::minimizer::LBFGSStep{/*...*/});
opt.on(minimizer::Phase::force).add(
    cpu::scf::SCFEvaluator{cpu::scf::b3lyp_config(integrals)});
opt.on(minimizer::Phase::end_of_step).add(
    cpu::minimizer::TrajectoryWriter{"aspirin_b3lyp.xyz", 1});
opt.build();
opt.run(frame, 200);
```

**Acceptance (3.3):**
- **Correctness:** B3LYP/6-31G* optimized geometry RMS deviation from PySCF ≤ 1e-3 Å on {aspirin, caffeine, glycine, benzene}.
- **Correctness:** final B3LYP energy matches PySCF to ≤ 1e-6 Ha on same molecules.
- **Performance:** full opt wall time ≤ 2× PySCF `geometric` optimizer on aspirin.
- Test: `tests/cpu/scf/test_hybrid_opt_vs_pyscf.cpp` (may be tagged `heavy` — runs nightly, not on every push).

## Files to modify

| File | Action |
|---|---|
| `include/atomiverse/cpu/scf/xc_functional.h` | add `XCFunctional`, `XCSpec`, `make_xc_spec` |
| `src/cpu/scf/xc_functional.cpp` | implement hybrid XC energy-density blending |
| `include/atomiverse/cpu/scf/scf_solve.h` | `SCFConfig.xc` replaces scalar `XCType`; add `b3lyp_config` / `pbe0_config` |
| `src/cpu/scf/components/hamiltonian_builder.cpp` | use `alpha_hf` when composing F |
| `src/cpu/scf/integral_derivatives.cpp` | use `config.xc.alpha_hf` in K-channel factor |
| `tests/cpu/scf/test_hybrid_scf_vs_pyscf.cpp` | new |
| `tests/cpu/scf/test_hybrid_gradient_vs_pyscf.cpp` | new |
| `tests/cpu/scf/test_hybrid_opt_vs_pyscf.cpp` | new |
| `benchmarks/qc/reference/generate_references.py` | add B3LYP/PBE0 reference JSONs |
| `examples/tutorial/07_b3lyp_opt.cpp` | new (composition example) |

## Dependencies

- **Phase 2 required:** `SinglePoint` driver + `HamiltonianBuilder` component as the integration point.
- **Phase 1 required:** 6-31G* basis (for target comparison with PySCF default B3LYP papers) and def2-SVP (for cross-checks).

## Risks

- **B3LYP convention ambiguity:** Gaussian's B3LYP uses VWN5 correlation; PySCF default is VWN3 (a.k.a. VWN-III). Ours will match **PySCF**. Any user comparing against Gaussian output will see differences at 1e-4 Ha level — document prominently.
- **Hybrid K-derivative performance:** exact exchange derivatives are the dominant cost (no density fitting in scope). Expect 2× PySCF wall time; flagged in the performance acceptance. If acceptance fails on aspirin, escalate — either optimize or document as known.
- **α=0 regression:** setting `alpha_hf = 0` on a hybrid code path must reproduce the pure-DFT result bit-for-bit. Add as explicit regression test.

## Exit gate

All 3.1–3.3 acceptance criteria green + pure-DFT regression bit-exact + user sign-off.

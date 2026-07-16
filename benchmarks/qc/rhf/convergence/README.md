# qc/rhf/convergence — RHF energy vs basis set

Convergence study for the manual "basis sets & convergence" chapter
(manual-system §3 chapter 4, listed there as a **new** basis-convergence
benchmark). A third study type alongside the case's `validation/` (vs PySCF)
and `performance/` (vs libcint) parts.

## What it checks

`bench_qc_basis_convergence` runs a closed-shell **RHF single point on a fixed
H2O geometry** across the systematically-improving basis ladder

```
STO-3G → 6-31G → 6-31G* → def2-SVP → def2-TZVP
```

and asserts the **physical-approximation** result class (manual-system §5): not
agreement with one reference number, but the *trend*:

1. **Variational monotonicity** — a larger basis lowers (or, within a small
   cross-family slack, holds) the HF energy. The Rayleigh–Ritz variational
   principle guarantees `E(larger ⊇ smaller) ≤ E(smaller)` for a strict
   superset; across *families* (6-31G* vs def2-SVP) the sets are not nested, so
   a small positive `--slack` (default 2e-3 Ha) absorbs legitimate
   near-degeneracy while still catching a basis that returns nonsense.
2. **Convergence** — the successive energy lowering shrinks toward the
   complete-basis-set (CBS) limit (reported as `dE_vs_prev`, informational).

Exit code 0 = monotone-convergent trend; 1 = a rung raised the energy past the
slack or an SCF failed to converge.

## Why a skeleton

The physics is real and runs today (all five bases are supported —
`cpu::scf::BasisSet`), but this benchmark's *added* value over the
PySCF-validated `validation/` suite is the cross-basis trend. It is deliberately
minimal. Extend as the chapter needs:

- more molecules (CH4, NH3, benzene) to show the trend is systematic;
- a two-point CBS extrapolation (e.g. def2-SVP/def2-TZVP → estimated limit);
- a companion **DFT** ladder (`SCFSolver` / PBE) once chapter 5 needs it.

## Expected numbers (orientation, not a gate)

RHF/H2O total energy is monotone-decreasing down the ladder; STO-3G is ~0.5–0.7
Ha above the def2-TZVP value, and def2-SVP → def2-TZVP typically lowers the
energy by only a few mHa (the trend is flattening — the hallmark of
convergence). The exact numbers are whatever the current integrals + SCF
produce; the `validation/` suite is what pins absolute energies to PySCF.

## Running

```bash
cmake -B build -DATV_BUILD_RHF=ON -DATV_BUILD_BENCHMARKS=ON ...
cmake --build build --target bench_qc_basis_convergence
./build/benchmarks/bench_qc_basis_convergence            # default slack 2e-3 Ha
./build/benchmarks/bench_qc_basis_convergence --slack 5e-3
```

The geometry is embedded in the source (no data-directory dependency), so the
binary runs from any working directory.

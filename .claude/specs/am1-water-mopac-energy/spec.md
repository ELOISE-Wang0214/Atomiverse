# Spec — `am1-water-mopac-energy`

Wave 4, `tb-method-loop`, Phase B (NDDO product depth, oracle = **OpenMOPAC**).

## Summary

The AM1/RM1 NDDO surface on `dev` ships an energy path but **no hard OpenMOPAC
energy pin anywhere**: `reg_am1_water_force` gates a *self*-consistency identity
(analytic force vs its own re-SCF finite difference) plus a ±8 kcal/mol window
around ΔHf, and `reg_rm1_water_energy` is method-discriminating only. Nothing in
the tree would catch a systematic AM1 parameter or integral defect that moved
the absolute energy by, say, 0.5 eV.

Deliver the first **absolute OpenMOPAC pin** on the NDDO surface: H₂O at a fully
specified geometry, run through the public `cpu::semiempirical::AM1` /
`AM1Evaluator` API, gated against a published MOPAC AM1 output — total energy
[eV], heat of formation [kcal/mol], and Mulliken charges — plus the geometry
convention that makes the pin defensible.

Also seed `tools/tb_oracle/inputs/mopac/`, which does not exist on `dev`.

## Domain basis

**Geometry convention (the load-bearing decision).** A published AM1 ΔHf is
almost always quoted **at the AM1-optimized geometry**, so it cannot be pinned
at an arbitrary structure. Two options:

1. Pin at the experimental geometry → needs an oracle run at *that* geometry;
   no published MOPAC output for it, and `mopac` is not installed on this host.
2. Pin at the **AM1-optimized geometry, taken from the same published output**
   as the energy.

Option 2 is chosen. The archived MOPAC 7.00 output
(`AM1 PRECISE EF T=48H`, Zipse teaching archive, LMU München) prints the
optimized internal coordinates to 5–6 decimals *and* the converged energies in
the same transcript, so geometry and energy come from one self-consistent run.

Two things make this a *sharper* gate rather than a weaker one:

- The published structure is a **stationary point of MOPAC's AM1 surface**
  (`RMS GRADIENT = 0.00141` kcal/mol/Å ≈ 6×10⁻⁵ eV/Å). If ATV's AM1 is the same
  functional, that geometry must be a stationary point of **ATV's** surface too.
  So the regression additionally asserts ‖F‖ ≈ 0 there — a gradient-level
  agreement check the current AM1 regressions cannot make, and one that is
  *independent* of the energy pin.
- Stationarity also makes the pin robust to the truncation in the published
  coordinates: near a minimum the energy is quadratic, so the ≤5×10⁻⁶ Å / ≤10⁻⁵°
  rounding in the transcript perturbs the energy at second order (≪ 10⁻⁶ eV) and
  cannot be confused with a parameter defect.

**ΔHf convention.** MOPAC's ΔHf is
`ΔHf = (E_elec + E_core − Σ_A EISOL_A)·23.061 + Σ_A ΔHf_atom(A)` [kcal/mol]
(openmopac.net/Manual/calc_of_hof.html) — already implemented in
`src/cpu/semiempirical/nddo_solver.cpp`. AM1 `EISOL(H) = −11.396427` eV,
`EISOL(O) = −316.09952` eV, `ΔHf_atom(H) = 52.102`, `ΔHf_atom(O) = 59.559`
kcal/mol — all already in `am1_parameters.cpp` and matching MOPAC.
Pinning ΔHf therefore tests the SCF energy **and** the thermochemistry
constants, and pinning the total energy separately isolates which one broke.

**Oracle family lock:** NDDO → OpenMOPAC only. DFTB+ is a forbidden pin source
for this row.

## Design

One out-of-library regression + one unit test file. **No production change is
planned** — the pin is a gate on existing behaviour. If ATV disagrees beyond a
physically defensible tolerance, that is a finding to report, not a tolerance to
loosen (iron law).

Public API used: `cpu::Frame`, `cpu::Evaluate`,
`cpu::semiempirical::AM1Evaluator`, `cpu::semiempirical::AM1`,
`frame.meta.get(kEnergy)`, `AM1Result::{total_energy_ev,
heat_of_formation_kcal_mol, mulliken_charges, forces_ev_a}`. No `*Config`, no
factory, Frame in / Frame out.

## Files

| Path | Action |
| ---- | ------ |
| `regressions/semiempirical/am1-water-mopac-energy/reg_am1-water-mopac-energy.cpp` | new |
| `regressions/semiempirical/am1-water-mopac-energy/REFERENCE.md` | new |
| `regressions/CMakeLists.txt` | register `atv_add_regression(semiempirical am1-water-mopac-energy)` |
| `tools/tb_oracle/inputs/mopac/am1-water/mopac_am1_water.mop` | new — input deck twin |
| `tools/tb_oracle/inputs/mopac/am1-water/mopac_am1_water.out.txt` | new — archived transcript |
| `tools/tb_oracle/inputs/mopac/am1-water/README.md` | new — offline-twin note |
| `tools/tb_oracle/README.md` | edit — record the mopac tree |
| `tests/cpu/semiempirical/test_am1_thermochemistry.cpp` | new — unit gates on the ΔHf constants |
| `tests/CMakeLists.txt` | register the new unit test |
| `docs/developer/tb-method-queue.md` | move `am1-water-mopac-energy` pending → done |

## Tasks

1. Transcribe the MOPAC deck + output into `tools/tb_oracle/inputs/mopac/am1-water/`.
2. Unit gates: `am1_atomic_heat_kcal_mol` / `eisol` / `kAm1EvToKcalMol` equal the
   MOPAC constants (single-function, milliseconds).
3. Regression: build the frame from the published internal coordinates, run the
   public Evaluate + `AM1Evaluator` path, assert E, ΔHf, Mulliken, stationarity.
4. Measure ATV's values; set tolerances from the *transcript precision*, not
   from the observed residual.
5. Register in `regressions/CMakeLists.txt` + `tests/CMakeLists.txt`; write
   `REFERENCE.md` with `oracle: mopac`, `oracle_version`, geometry, pins,
   tolerance rationale, and `Gaps`.

## Testing

- Unit: `ctest -R '^test_cpu_am1'` (seconds).
- Regression: `cmake --build build --target run_reg_am1-water-mopac-energy`.
- Domain sanity: `ctest -R '^test_cpu_(am1|rm1|dftb0|dftb2|dftb3|xtb_gfn0|nddo_)'`.

## Out of scope

- Closed-form NDDO multipole forces (`nddo-closed-form-forces`).
- RM1 mopac pin (`rm1-water-mopac-energy`), PM3/PM6/PM7 (no parameter tables).
- Geometry optimization in-case: `ATV_BUILD_MINIMIZER` is **OFF** in the SE-only
  build, so the case pins at the published optimum rather than re-optimizing.
- Koopmans IP / eigenvalue pin: `NDDOResult` exposes no orbital energies, and
  adding a public surface for one pin is out of scope for this row.

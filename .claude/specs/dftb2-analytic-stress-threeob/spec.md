# spec — `dftb2-analytic-stress-threeob`

Closed-form cell stress on `DftbSkSet::ThreeOb`; delete the multi-SCC affine
strain-FD product path for that SK set.

## Summary

`DFTB2::with_stress` today gives a closed-form pair virial **only** on the
restricted MioSketch path (`analytic_forces_restricted`,
`src/cpu/semiempirical/dftb2.cpp:1364`). On `DftbSkSet::ThreeOb`,
`launch_3ob` (`src/cpu/semiempirical/dftb2.cpp:2364`) falls back to **twelve**
full re-SCC energy evaluations (6 Voigt × ±h) of an affine strain of the
converged geometry.

ThreeOb nuclear forces are already closed-form on this branch
(`analytic_forces_3ob`, `src/cpu/semiempirical/dftb2.cpp:2184`, landed as
`dftb2-analytic-forces-threeob`), so every virial term is already computed:
Hellmann–Feynman ∂H⁰/∂R, Pulay −∂S/∂R contracted with the energy-weighted
density matrix W, the fixed-Δq SCC γ pair gradient, and the 3ob repulsive
spline derivative.

Deliver the closed-form ThreeOb virial, then **delete** the strain-FD product
path (iron law — no dead slow path left behind); an FD routine survives only
inside the tests, as the oracle.

## Domain basis

SCC-DFTB2 (Elstner et al., Phys. Rev. B **58**, 7260 (1998)). On the ThreeOb
path the total is assembled as

    E = Tr(P H⁰) + ½ Σ_AB γ_AB Δq_A Δq_B + E_rep

(`scc_3ob`, all in Ha / Bohr; the public surface is eV / Å).

The **stress theorem** (Nielsen & Martin, Phys. Rev. B **32**, 3780 (1985)) in
the mechanical convention used by DFTB+ (Hourahine et al., J. Chem. Phys.
**152**, 124101 (2020)) and by ATV's PW / GPW `kStress`:

    σ_αβ = −(1/Ω) ∂E/∂ε_αβ

Under an affine strain of nuclei **and** cell, every bond vector transforms as
R → (1+ε)R, and the ThreeOb energy depends on geometry **only** through the
pair bond vectors (`scc_3ob` is a raw-distance cluster sum — no lattice images,
no min-image). Therefore

    ∂E/∂ε_αβ = Σ_pairs (∂E/∂R_α) R_β

i.e. exactly the pair virial of the same `dE_dbond` vector that
`analytic_forces_3ob` already forms per (a,b) pair, contracted with that pair's
bond vector. This is FD-exact, not an approximation: the FD path it replaces
differentiates the *same* functional.

Sign / Voigt convention is taken **from the restricted-Mio path already in
tree** (`analytic_forces_restricted`, `dftb2.cpp:1473-1491`), not assumed:

    dE_de[0] += dEdb_x R_x                    (xx)
    dE_de[1] += dEdb_y R_y                    (yy)
    dE_de[2] += dEdb_z R_z                    (zz)
    dE_de[3] += dEdb_x R_y + dEdb_y R_x       (xy)
    dE_de[4] += dEdb_x R_z + dEdb_z R_x       (xz)
    dE_de[5] += dEdb_y R_z + dEdb_z R_y       (yz)
    σ[v]      = −dE_de[v] / Ω

(The off-diagonal doubling is because `apply_affine_strain` builds a *symmetric*
strain ε_αβ = ε_βα = h for Voigt 3/4/5.)

Units: `analytic_forces_3ob` works in Ha/Bohr with bond vectors in Bohr, so
`dE_de` accumulates in **Ha** (strain is dimensionless); the public value is
`σ = −dE_de·kHa2eV / Ω_Å³` [eV/Å³].

**Isotropic identity used by the unit gate.** For isotropic scaling
r → (1+ε)r, Ω(ε) = (1+ε)³Ω₀ so dΩ/dε|₀ = 3Ω₀ and
dE/dε = Σ_α ∂E/∂ε_αα = −Ω₀ (σ_xx+σ_yy+σ_zz). Hence

    tr(σ)/3 = −dE/dΩ = −(1/3Ω₀) dE/dε|₀

which the test checks against a central difference of `E` over isotropic
affine scaling.

## Design

Extend `analytic_forces_3ob` (file-local, anonymous namespace) exactly the way
`analytic_forces_restricted` is already extended — no new function, no new
public symbol, no `*Config`:

```cpp
void analytic_forces_3ob(..., AccReal* fx, AccReal* fy, AccReal* fz,
                         bool want_forces,
                         double* stress_voigt /* nullable, length 6 */);
```

Inside the existing `a < b` pair loop, after `dE_dbond[]` is complete (SK
derivative + γ gradient + repulsive spline derivative), accumulate the six
Voigt products above; after the loop convert with `kHa2eV / Ω`.

`launch_3ob`:

- call `analytic_forces_3ob` when `with_forces || with_stress`;
- write `out_stress` and accumulate into `f.meta.get(kStress)` as the Mio path
  does;
- **delete** the strain-FD block and its only helper `energy_3ob_ha`
  (`dftb2.cpp:2171`), which exists solely to feed it.

`apply_affine_strain` stays — the MioSketch **open-shell** path still uses it
(collinear open-shell stress is a separate queued slug,
`dftb2-spin-analytic-stress`).

Restricted-Mio code is **not touched**; a unit test pins its stress bitwise.

## Files

| Path | Change |
| ---- | ------ |
| `src/cpu/semiempirical/dftb2.cpp` | `analytic_forces_3ob` gains `want_forces` + `stress_voigt`; virial accumulation; `launch_3ob` rewired; strain-FD block + `energy_3ob_ha` deleted |
| `include/atomiverse/cpu/semiempirical/dftb2.h` | class/member Doxygen: ThreeOb stress is closed-form |
| `tests/cpu/semiempirical/test_dftb2_3ob.cpp` | new unit gates (below) |
| `tests/cpu/semiempirical/test_dftb2.cpp` | Mio-path stress fossil guard (≤1e-12) |
| `regressions/semiempirical/dftb2-analytic-stress-threeob/` | new `reg_*.cpp` + `REFERENCE.md` |
| `regressions/CMakeLists.txt` | `atv_add_regression(semiempirical dftb2-analytic-stress-threeob)` |
| `tools/tb_oracle/inputs/dftbplus/dftb2_analytic_stress_threeob_note.md` | offline pin note |
| `docs/roadmap/semiempirical.md`, `docs/developer/tb-method-queue.md` | ThreeOb stress ○ → ✓ |
| `docs/tutorial/methods/dftb2-3ob-organic.md`, `dftb2-analytic-forces.md` | stale "ThreeOb forces remain re-SCC FD" claims (already false since `dftb2-analytic-forces-threeob` landed) |

## Tasks

1. Baseline: build + record the restricted-Mio stress at full precision (fossil guard value).
2. `analytic_forces_3ob`: add `want_forces` + `stress_voigt`, accumulate the pair virial.
3. `launch_3ob`: route `with_stress` through it; delete the FD block + `energy_3ob_ha`.
4. Header Doxygen.
5. Unit gates in `tests/cpu/semiempirical/test_dftb2_3ob.cpp`.
6. Mio fossil guard in `tests/cpu/semiempirical/test_dftb2.cpp`.
7. Regression dir + `REFERENCE.md` + CMake registration + oracle note.
8. Docs / roadmap / queue.

## Testing

**Unit** (`tests/cpu/semiempirical/`, seconds, no e2e):

- `AnalyticStressMatchesStrainFdCh4` — **distorted** CH4 (one H displaced so
  all six Voigt components are non-zero; symmetric CH4 has identically zero
  shear) in a 20 Å cubic box; all six analytic σ vs affine strain-FD of the
  converged ThreeOb energy, ≤ **1e-8 eV/Å³**.
- `AnalyticStressMatchesStrainFdWater` — off-eq H₂O (the polar Δq ≠ 0 case that
  caught the γ-sign bug) — same gate.
- `AnalyticStressIsotropicScalingIdentity` — tr(σ)/3 vs −dE/dΩ from an
  isotropic-scaling central difference.
- `test_dftb2.cpp` `MioAnalyticStressFossilGuard` — restricted MioSketch H₂O
  stress vs hard-coded pre-change values, ≤ **1e-12 eV/Å³**.

**Regression** `regressions/semiempirical/dftb2-analytic-stress-threeob/` —
distorted CH₄ with 3ob-3-1 through the public `DFTB2` API (out-of-library
consumer, public headers only, ≤ 60 s): finiteness, analytic-vs-strain-FD
identity ≤ 1e-8 eV/Å³, isotropic trace identity. `REFERENCE.md` with
`oracle: dftb+`, `oracle_version: 24.1` (family oracle; **formula-level** — no
DFTB+ binary on this host), geometry + cell, pinned values, citations,
tolerance rationale, `Gaps`.

## Out of scope

- Collinear **open-shell** analytic stress (queued `dftb2-spin-analytic-stress`).
- Multi-k forces / stress (queued `dftb2-multi-k-forces`).
- D3 dispersion virial (open gap on `dftb2-analytic-stress` REFERENCE).
- Periodic (lattice-image) ThreeOb: `scc_3ob` is a molecular cluster sum; the
  stress inherits exactly that limitation, unchanged by this spec.
- Absolute DFTB+ σ pin (no oracle binary on this host).

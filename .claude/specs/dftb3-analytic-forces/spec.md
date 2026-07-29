# dftb3-analytic-forces

status: code_complete
branch: `wf/dftb3-analytic-forces`
loop: tb-method-loop, wave 4
oracle family: **DFTB+** (no binary on host → formula / internal-consistency gates)

## Summary

`cpu::semiempirical::DFTB3::launch` computes nuclear forces by **re-SCC central
finite difference** of the total energy: `6N` extra full self-consistency solves
per force evaluation (`src/cpu/semiempirical/dftb3.cpp:360-386`). Replace that
product path with the **closed-form gradient at the converged SCC density**, the
same assembly DFTB2 already ships (`analytic_forces_restricted`,
`src/cpu/semiempirical/dftb2.cpp:1364`), extended by the third-order term.

FD stays available as the **test oracle only** — the public energy path
(`with_forces = false` → `last_energy()` / `frame.meta[kEnergy]`) re-SCCs at every
displaced geometry, which is exactly what the SCC-consistency gate needs.

## Domain basis

DFTB3 total energy (Gaus, Cui, Elstner, *J. Chem. Theory Comput.* **7**, 931
(2011), DOI 10.1021/ct100684s; second-order backbone Elstner et al.,
*Phys. Rev. B* **58**, 7260 (1998)):

```
E = Tr(ρ H⁰) + E⁽²⁾ + E⁽³⁾ + E_rep
E⁽²⁾ = ½ Σ_AB γ_AB Δq_A Δq_B
E⁽³⁾ = ⅓ Σ_AB Δq_A² Δq_B Γ_AB ,   Γ_AB = ∂γ_AB/∂q_A  (Γ_AB ≠ Γ_BA)
```

The SCC Hamiltonian is `H = H⁰ + ½ S (V_A + V_B)` with
`V = ∂(E⁽²⁾+E⁽³⁾)/∂Δq` — i.e. `H_μν = ∂E/∂P_μν` exactly (Mulliken
`∂Δq_A/∂P_μν = ½S_μν(δ_{μ∈A}+δ_{ν∈A})`). The Euler-correct `V⁽³⁾` (both
contractions of the non-symmetric Γ) is already on this branch in
`internal/dftb3_scc.h`. Because `H` **is** the derivative of the energy w.r.t.
the density, the converged `ρ` is stationary and the nuclear gradient is the
ordinary Hellmann–Feynman + Pulay expression — no coupled-perturbed charge
response:

```
∂E/∂R = Σ_μν P_μν ∂H⁰_μν/∂R
      + Σ_μν P_μν ½(V_A+V_B) ∂S_μν/∂R
      − Σ_μν W_μν ∂S_μν/∂R                  W_μν = Σ_o n_o ε_o C_μo C_νo
      + ∂(E⁽²⁾+E⁽³⁾)/∂R |_Δq
      + ∂E_rep/∂R
```

with the fixed-Δq charge-kernel term, per unordered pair A≠B and bond
`R = R_B − R_A`:

```
∂E_scc^{AB}/∂R = Δq_A Δq_B ∂γ_AB/∂R
               + ⅓ ( Δq_A² Δq_B ∂Γ_AB/∂R + Δq_B² Δq_A ∂Γ_BA/∂R )
```

For the in-tree Ohno–Klopman-like kernel `γ_AB = (R² + τ²)^{-1/2}`,
`τ = ½(1/U_A + 1/U_B)`:

```
∂γ_AB/∂R = −R γ_AB³
Γ_AB     = τ (R²+τ²)^{-3/2} (½/U_A²) U^d_A
∂Γ_AB/∂R = −3 R τ (R²+τ²)^{-5/2} (½/U_A²) U^d_A
```

On-site `γ_AA = U_A`, `Γ_AA = U^d_A` are geometry-independent ⇒ zero gradient.

## Design

**No duplication across methods** (CLAUDE.md): the ~85-line MioSketch SK block
value + derivative kernel currently lives in `dftb2.cpp`'s anonymous namespace
(`fill_sk_block_directed`, `fill_sk_block_deriv`) and is *also* re-implemented as
a value-only copy in `dftb3.cpp` (`fill_sk_block`). Copying the derivative into
`dftb3.cpp` would be a third copy of the subtle direction-cosine algebra. Extract
it once instead.

1. **New private workspace type `internal::MioSketchPair`**
   (`internal/mio_sk_block.h` + `src/cpu/semiempirical/mio_sk_block.cpp`),
   mirroring the existing `internal::MioHHTable` reader pattern: binds one
   ordered element pair + its static sketch SK row and evaluates, for a bond
   vector `R = R_B − R_A`,
   - `add_block(...)` — accumulate the directed A→B `H⁰`/`S` block,
   - `set_block_derivative(...)` — assign the directed A→B `∂H⁰/∂R`, `∂S/∂R`
     blocks.
   Both carry the mio-1-1 H–H special case (shared `MioHHTable`). Value and
   derivative now live in one type, so they cannot drift apart.
   Internal only (`internal/`), never on a public include path.

2. **`dftb2.cpp`** drops its two local copies and calls `MioSketchPair`.
   Bodies move verbatim → energies and forces are bit-identical (guarded by the
   existing DFTB2 regressions).

3. **`dftb3.cpp`** drops `fill_sk_block` (its H⁰/S build calls
   `MioSketchPair::add_block` on the `a<b` min-image walk, then mirrors the
   strict upper triangle — the only entries that walk writes) and gains
   `analytic_forces(...)`, a closed-form gradient over the same min-image pair
   walk the energy uses: SK-derivative terms inside `pair_rcut_A`, γ/Γ and
   repulsive terms with no cutoff (exactly as `solve_scc` assembles them).
   `with_forces` now selects it; the FD loop is deleted.

4. **`internal::Dftb3Scc`** gains the geometry derivatives of its own kernel,
   next to the energies they differentiate:
   - `gamma_distance_derivative(a,b)` → `∂γ_AB/∂R` [eV/Å]
   - `gamma_charge_distance_derivative(a,b)` → `∂²γ_AB/∂q_A∂R` [eV/(e·Å)]
   - `pair_energy_distance_derivative(a,b,dq)` → `∂(E⁽²⁾+E⁽³⁾)_AB/∂R` [eV/Å],
     owning the ⅓ and the two Γ orderings so the gradient cannot disagree with
     `third_order_energy` about the convention.

5. `SccSolution` carries the converged MO coefficients / eigenvalues and the
   frozen charge kernel so the gradient consumes exactly the solve that
   converged (no second SCC, no re-derived γ).

**No new public surface**: `DFTB3`'s public members are unchanged
(`with_forces` keeps its name and meaning; only its documented implementation
changes). No `*Config` bag, no factory.

## Files

| Path | Change |
| ---- | ------ |
| `include/atomiverse/cpu/semiempirical/internal/mio_sk_block.h` | **new** — `MioSketchPair` |
| `src/cpu/semiempirical/mio_sk_block.cpp` | **new** — its impl (moved verbatim from `dftb2.cpp`) |
| `include/atomiverse/cpu/semiempirical/internal/dftb3_scc.h` | +3 kernel-derivative methods |
| `src/cpu/semiempirical/dftb3_scc.cpp` | their impl |
| `src/cpu/semiempirical/dftb2.cpp` | delete 2 local SK helpers → call `MioSketchPair` |
| `src/cpu/semiempirical/dftb3.cpp` | delete FD force loop + local `fill_sk_block`; add `analytic_forces` |
| `include/atomiverse/cpu/semiempirical/dftb3.h` | docs: forces are closed-form; cite Gaus 2011 |
| `CMakeLists.txt` | register `mio_sk_block.cpp` |
| `tests/cpu/semiempirical/test_dftb3.cpp` | + analytic-vs-FD, ΣF, h-refinement, DFTB2-reduction force gates |
| `tests/cpu/semiempirical/test_dftb3_scc.cpp` | + FD gates on the 3 new kernel derivatives |
| `regressions/semiempirical/dftb3-analytic-forces/` | **new** case + `REFERENCE.md` |
| `regressions/CMakeLists.txt` | `atv_add_regression(semiempirical dftb3-analytic-forces)` |

## Tasks

1. `MioSketchPair` extraction; `dftb2.cpp` + `dftb3.cpp` rewired; DFTB2 suite +
   `reg_dftb3-water-mio` energy pin unchanged (bit-identical).
2. `Dftb3Scc` geometry derivatives + their FD unit gates.
3. `analytic_forces` in `dftb3.cpp`; FD product path deleted.
4. Unit gates in `test_dftb3.cpp`.
5. Regression `dftb3-analytic-forces` + `REFERENCE.md` + CMake registration.

## Testing

**Unit** (`tests/cpu/semiempirical/`, seconds, one method per test):

- `Dftb3Scc.GammaDistanceDerivativeMatchesFd` / `…ChargeDistanceDerivative…` /
  `…PairEnergyDistanceDerivative…` — central FD over two kernels built at `R±h`.
- `DFTB3.AnalyticForceMatchesFdH2O` — all 9 components of off-equilibrium water
  vs central FD (`h = 1e-4` Å) of the **re-SCC** converged energy; ≤ 1e-6 eV/Å.
- `DFTB3.AnalyticForceSumsToZeroH2O` — ΣF ≤ 1e-10 eV/Å per Cartesian
  (translational invariance; exact by pairwise assembly).
- `DFTB3.AnalyticForceFdResidualShrinksWithStep` — the SCC-consistency gate: the
  analytic-vs-FD residual must fall like `h²` (`h = 8e-4` → `2e-4`). A missing or
  mis-signed term leaves an `h`-independent offset and fails here even when a
  single loose-tolerance comparison would pass.
- `DFTB3.AnalyticForceReducesToDftb2ForHomonuclearH2` — `Δq ≡ 0` ⇒
  `E⁽²⁾=E⁽³⁾=0` ⇒ the DFTB3 gradient must equal the DFTB2 gradient to
  1e-12 eV/Å.
- `DFTB3.AnalyticForceMatchesFdCOWithPBlocksOnBothCentres` — water pairs O(s,p)
  only with H(s), so it never reaches the pp-σ/π SK derivative or the off-axis
  direction-cosine terms; off-axis C–O does (6 components, ≤ 1e-6 eV/Å).

**Regression** `regressions/semiempirical/dftb3-analytic-forces/` — out-of-library
consumer (`cpu::Evaluate` + `DFTB3`), off-equilibrium H2O, ≤ 60 s: finite E,
`|F|² > 1e-8`, `|ΣF| ≤ 1e-6`, analytic vs re-SCC FD on `F_O,x` / `F_H1,y` within
1e-5 eV/Å, third-order force actually differs from DFTB2, and a provisional
ATV-twin force pin. `REFERENCE.md` states `oracle: dftb+`, DFTB+ 24.1 as the
family oracle at formula level, and records in `Gaps` that **no DFTB+ binary
exists on this host** so no absolute oracle force is claimed.

## Out of scope

- Periodic (Γ / lattice-image) or multi-k DFTB3 forces — DFTB3 is min-image
  molecular today; the periodic walk is DFTB2's (`dftb2-periodic-forces`).
- DFTB3 cell stress, open shell (spin), 3ob tables for DFTB3, D3 dispersion.
- Absolute DFTB+ energy/force identity (needs an offline DFTB+ pin under
  `tools/tb_oracle/inputs/dftbplus/`).
- Any change to the DFTB3 **energy** — `reg_dftb3-water-mio`'s pin must not move.

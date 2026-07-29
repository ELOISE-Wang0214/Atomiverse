# DFTB MioSketch γ_AB is dimensionally inconsistent (found 2026-07-29)

**Status:** open. Found by the `td-dftb-h2` worker in TB wave 4, independently
re-derived at merge. **Not fixed** — the fix is two lines but re-pins every
polar MioSketch case, so it is its own change (`/mol:fix`), not a drive-by.

## The defect

`src/cpu/semiempirical/dftb2.cpp` `gamma_ab` and the identical
`src/cpu/semiempirical/dftb3_scc.cpp` `Dftb3Scc::gamma`:

```cpp
double gamma_ab(double r, double ua, double ub) {
    const double tau = 0.5 * (1.0 / ua + 1.0 / ub);
    return 1.0 / std::sqrt(r * r + tau * tau);
}
```

`hubbard_u` is **[eV]** and `r` is **[Å]**, so `tau` is [eV⁻¹] and the sum adds
**Å² to eV⁻²**. The Coulomb constant is missing: `units::coulomb = 14.399645`
eV·Å already exists in `include/atomiverse/core/units.h`.

The consistent form (Ohno–Klopman written in a.u., where k_e = 1, then carried
to eV/Å) is

```cpp
γ = k_e / std::sqrt(r*r + (k_e*tau)*(k_e*tau));
```

## Why it survived every pin

| limit | shipped | consistent | note |
| ----- | ------- | ---------- | ---- |
| R → 0 | 1/τ = harmonic mean of U | k_e/(k_e τ) = 1/τ | **identical** — γ_AA = U_A is accidentally right |
| R → ∞ | 1/R [Å⁻¹ read as eV] | k_e/R | **14.4× too small** |
| crossover | τ ≈ 0.088 Å | k_e τ ≈ 1.26 Å | wrong scale |

H–H at 0.74 Å: **1.342 eV** shipped vs **9.856 eV** consistent (7.3×); at 3 Å,
13× low.

Invisible because the only real-mio-1-1 system in tree is H₂, where Δq ≡ 0 by
symmetry so γ_AB never enters the energy — and every polar case is declared a
"sketch" whose pins are ATV-twin drift guards, not oracle numbers. The one
absolute DFTB+ pin (`dftb2-dftbplus-energy-pin`, H₂) is Δq = 0 and therefore
unaffected by the fix.

## Blast radius when fixed

Moves every MioSketch DFTB2/DFTB3 energy, force, stress and charge for any
Δq ≠ 0 system. Re-pins: `h2o_dftb2_mio`, `dftb2-analytic-forces`,
`dftb2-analytic-stress`, `dftb2-periodic-gamma`, `dftb2-periodic-forces`,
`dftb2-kpoints`, `dftb2-multi-k-forces`, `dftb2-dispersion-d3`,
`dftb2-dispersion-d3-periodic`, `spin-dftb2-o2`, `dftb2-spin-analytic-forces`,
`dftb2-spin-analytic-stress`, `dftb3-water-mio`, `dftb3-analytic-forces`,
`td-dftb-h2`.

**All of those gates are invariant under the fix in *form*** — they are
analytic-vs-FD identities, symmetry limits and extensivity identities, not
absolute values. Only the quoted ATV-twin numbers move. That is what makes
this a mechanical (if wide) re-pin rather than a redesign.

Does **not** touch: `dftb0-h2` (non-SCC), the ThreeOb path (3ob has its own
γ), NDDO/AM1/RM1, GFN0.

## Route

`/mol:fix` as one change: correct both call sites, re-derive the drift pins in
a single pass, confirm `dftb2-dftbplus-energy-pin` is bit-unchanged, and check
whether the MioSketch U values were tuned against the broken γ (they are sketch
values, so a re-tune may be owed too — that part is *not* mechanical).

Related: [[gpw-periodic-lattice-sum]], [[e2e-parity-goes-in-benchmarks]].

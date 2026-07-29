# DFTB MioSketch γ_AB is dimensionally inconsistent (found 2026-07-29)

**Status:** **superseded 2026-07-29** by the DFTB+ **Elstner** γ kernel
(`internal/dftb_gamma.h`). The eV·Å Ohno–Klopman form below fixed the
dimensionless bug but is **not** the DFTB+ short-range model — it left
H₂O E² ≈ 0.65 eV vs DFTB+ 0.496 eV at the same Δq. Absolute H₂ pin never
moved (Δq ≡ 0). After Elstner: water \|ΔE\|≈3e-5 eV, 3ob CH₄ \|ΔE\|≈4e-7 eV.
See [[dftb-vs-dftbplus-cleanup]].

## The defect (historical)

`src/cpu/semiempirical/dftb2.cpp` `gamma_ab` and the identical
`src/cpu/semiempirical/dftb3_scc.cpp` `Dftb3Scc::gamma` previously:

```cpp
double gamma_ab(double r, double ua, double ub) {
    const double tau = 0.5 * (1.0 / ua + 1.0 / ub);
    return 1.0 / std::sqrt(r * r + tau * tau);
}
```

`hubbard_u` is **[eV]** and `r` is **[Å]**, so `tau` is [eV⁻¹] and the sum adds
**Å² to eV⁻²**. The Coulomb constant is missing: `units::coulomb = 14.399645`
eV·Å already exists in `include/atomiverse/core/units.h`.

## Fixed form

```cpp
const double ke = units::coulomb;
const double tau = 0.5 * (1.0 / ua + 1.0 / ub);
const double u = ke * tau;
return ke / std::sqrt(r * r + u * u);
// dγ/dR = -r * γ³ / ke²
```

ThreeOb path (`gamma_ab_ha` in Ha/Bohr) was already consistent (k_e = 1 in a.u.)
and was **not** changed.

## Why it survived every pin

| limit | shipped | consistent | note |
| ----- | ------- | ---------- | ---- |
| R → 0 | 1/τ = harmonic mean of U | k_e/(k_e τ) = 1/τ | **identical** — γ_AA = U_A is accidentally right |
| R → ∞ | 1/R [Å⁻¹ read as eV] | k_e/R | **14.4× too small** |

Invisible because the only real-mio-1-1 absolute pin is H₂, where Δq ≡ 0.

Related: [[tb-blocked-capabilities]], [[oracle-pin-provenance]].

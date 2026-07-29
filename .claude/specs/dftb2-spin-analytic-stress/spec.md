# dftb2-spin-analytic-stress

Closed-form collinear open-shell cell stress for `cpu::semiempirical::DFTB2`
(MioSketch), replacing the multi-SCC affine strain-FD product path.

## Summary

`DFTB2` with `multiplicity > 1` already has a closed-form nuclear gradient
(`analytic_forces_unrestricted`, wave-3). Its `with_stress` path still runs
**12 extra SCC solves** (6 Voigt × ±h) of affine strain-FD
(`src/cpu/semiempirical/dftb2.cpp:1660`). Deliver the spin-resolved pair
virial, delete the open-shell FD product path, and keep FD only as the test
oracle.

## Domain basis

Collinear spin SCC-DFTB2 (Köhler, Seifert, Frauenheim; Elstner et al. PRB 58,
7260 (1998); DFTB+ spin/stress module, Hourahine et al. JCP 152, 124101 (2020)):

E = Σ_σ Tr(P^σ H⁰) + ½ Σ_AB γ_AB Δq_A Δq_B + ½ Σ_A W_A m_A² + E_rep

with Δq_A, m_A the Mulliken charge / spin populations and the spin-resolved
Fock H^σ_μν = H⁰_μν + ½ S_μν (ε^σ_{A(μ)} + ε^σ_{A(ν)}),
ε^σ_A = Σ_B γ_AB Δq_B + σ W_A m_A  (σ = +1 α, −1 β).

**The functional is variational in P^α, P^β** — ∂E/∂P^σ = H^σ exactly, using
the symmetrized Mulliken derivative ∂Δq_A/∂P_μν = ½(δ_{μ∈A}+δ_{ν∈A}) S_μν.
Therefore, for any geometry parameter λ (nuclear coordinate **or** affine
strain ε), with C^{σT} S C^σ = I giving Tr(∂P^σ/∂λ · H^σ) = −Tr(W^σ ∂S/∂λ):

dE/dλ = Σ_σ [ Tr(P^σ ∂H⁰/∂λ)
             + Σ_μν P^σ_μν ½(ε^σ_{A(μ)}+ε^σ_{A(ν)}) ∂S_μν/∂λ
             − Tr(W^σ ∂S/∂λ) ]
        + ½ Σ_AB (∂γ_AB/∂λ) Δq_A Δq_B
        + ∂E_rep/∂λ ,          W^σ_μν = Σ_{o∈occ σ} ε^σ_o C^σ_μo C^σ_νo.

**Spin virial — it does not vanish.** E_spin = ½ Σ_A W_A m_A² carries no
*explicit* λ dependence (W_A is a geometry-independent atomic constant and m_A
is held at the SCC value), but m_A = m_A(λ) through S. The implicit part is
**not** dropped: it is exactly the ±W_A m_A half of ε^σ_A multiplying
∂S_μν/∂λ in line 2 above — that term is what the variational cancellation
converts the ∂m/∂λ response into. Omitting it (or replacing ε^σ by the
charge-only shift) breaks F = −∇E and σ = −(1/Ω)∂E/∂ε for any system whose
spin shift is non-uniform.

Stress: bond vectors R = R_B + T − R_A are the only λ-dependent quantities of
every pair term, and affine strain maps R → (1+ε)R, so

∂E/∂ε_αβ = Σ_pairs (∂E/∂R_α) R_β ,  σ = −(1/Ω) ∂E/∂ε  [eV/Å³],

Voigt (xx,yy,zz,xy,xz,yz) with the off-diagonal entries carrying both
ε_αβ = ε_βα = h halves — bit-identical assembly to the shipped restricted
path (`analytic_forces_restricted`).

## Design

**One routine, not two.** `analytic_forces_restricted` and
`analytic_forces_unrestricted` are already ~95 % duplicated; bolting a third
copy of the virial accumulation onto the open-shell one is what let the
wrong-sign γ gradient survive in three functions (see
`regressions/semiempirical/dftb2-analytic-forces/REFERENCE.md`). They collapse
into one spin-resolved `analytic_forces_stress(...)`:

- inputs: (C^α, ε^α, n_α), (C^β, ε^β, n_β), Δq, m, `want_forces`,
  `stress_voigt` (nullable);
- restricted call passes C^α = C^β = C, ε^α = ε^β = ε, n_α = n_β = n_occ,
  m ≡ 0 — the degenerate limit, algebraically identical to the shipped
  restricted assembly (P^α = P^β = P/2, W^α = W^β = W/2).

**Delete** the open-shell branch of the strain-FD block in `launch_mio`
(`dftb2.cpp:1660-1709`). The ThreeOb FD stress block stays (different SK
path, out of scope).

## Files

- `src/cpu/semiempirical/dftb2.cpp` — unify the two gradient assemblies into
  `analytic_forces_stress`; add the virial to it; drop the open-shell FD
  stress fallback in `launch_mio`.
- `include/atomiverse/cpu/semiempirical/dftb2.h` — `with_stress` /
  class-doc update (open shell is now closed-form).
- `tests/cpu/semiempirical/test_dftb2.cpp` — unit gates.
- `regressions/semiempirical/dftb2-spin-analytic-stress/` +
  `regressions/CMakeLists.txt` — new e2e case.

## Tasks

1. Unify the two assemblies into `analytic_forces_stress` with the spin-resolved
   pair virial; route both `launch_mio` branches through it.
2. Delete the open-shell multi-SCC strain-FD stress path.
3. Unit gates in `tests/cpu/semiempirical/test_dftb2.cpp`.
4. Regression `dftb2-spin-analytic-stress` + `REFERENCE.md` + CMake registration.
5. Header Doxygen update.

## Testing

Unit (`tests/cpu/semiempirical/test_dftb2.cpp`):

- O₂ triplet (`multiplicity = 3`), tilted bond, 20 Å cubic box: analytic σ vs
  affine strain-FD of the converged energy, **all six** Voigt components,
  ≤ 1e-8 eV/Å³; ΣF ≈ 0 unchanged.
- `multiplicity = 1` on the same geometry reproduces the shipped restricted
  analytic σ (baseline measured on the parent commit `448dca40`) to
  ≤ 1e-10 eV/Å³ — the equivalence guard for the unification.
- **NO doublet** (`multiplicity = 2`, Z = 8/7): the discriminating case.
  O₂ is homonuclear, so Δq ≡ 0 **and** the spin shift ±W_A m_A is uniform
  over atoms — a uniform cS shift leaves the generalized eigenvectors
  untouched, so the O₂ triplet gradient/virial is numerically identical to
  the restricted one and **cannot** see a dropped spin or γ virial. NO has
  Δq ≠ 0 and W_N m_N ≠ W_O m_O.

Regression `regressions/semiempirical/dftb2-spin-analytic-stress/`: O₂ triplet
through `cpu::Evaluate` + public `DFTB2` (mandated headline case) plus the NO
doublet discriminator, both analytic-vs-strain-FD; ≤ 60 s; `REFERENCE.md`
cites Köhler/Elstner/Hourahine with `oracle: dftb+`, formula-level
(no DFTB+ binary on this host).

## Out of scope

- ThreeOb stress (still strain-FD) — separate queue slug.
- Multi-k forces/stress (still throws).
- D3 dispersion virial (`kStress` excludes it on every path — owed).
- Self-image (a = a, T ≠ 0) virial under dense Γ-PBC: the shipped restricted
  path omits it too; documented gap, not introduced here.

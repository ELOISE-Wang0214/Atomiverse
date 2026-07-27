# Pin notes — pim-npt-polarizable-il

Slug: `pim-npt-polarizable-il` (metalwalls-method-loop wave 5).
Regression: `regressions/md/pim-npt-polarizable-il/`.

## Oracle

`oracle=metalwalls` (queue row). No MetalWalls binary is installed on this
host, so the pin follows the same offline-analytic path as
`pim-induced-dipole-bulk`: an **independent numpy implementation** of the
published bulk-PIM Thole multipole form —

- MetalWalls, *JOSS* **5**, 2373 (2020) — bulk PIM induced-dipole SCF form
  (charge–dipole + Thole-damped dipole–dipole, Hellmann–Feynman forces),
- Thole, *Chem. Phys.* **59**, 341 (1981) — exponential screening
  λ3/λ5/λ7 with u = a·r/(αiαj)^{1/6},

evaluated on the deterministic regression start configuration. LAMMPS is
**not** a fair twin (its polarizable paths are Drude / QEq, a different
algorithm), hence `perf_exempt=true`.

## Script

`tools/metalwalls_oracle/inputs/pim_npt_polarizable_il_pin.py`
(run offline: `~/work/.x86_64/bin/python3 pim_npt_polarizable_il_pin.py`,
numpy 2.5.1). Never executed at build/ctest/regression time.

## Configuration (must match the regression exactly)

- 4×4×4 simple-cubic rock-salt toy polarizable IL, d = 2.35 Å,
  cubic box L = 9.4 Å (V = 830.584 ų), index order ix→iy→iz.
- Charges (−1)^(ix+iy+iz): +1 cation (type 0, α = 0),
  −1 anion (type 1, α = 3.0 ų).
- Deterministic displacement per atom index i, axis a ∈ {0,1,2}:
  `disp = 0.15 · sin(1.7·(i+1) + 2.3·a)` Å (breaks the centro-symmetry
  that would make E_pol ≡ 0 on the perfect lattice).
- PIM: cutoff 4.0 Å, thole_a = 0.39, tt_b = 0, Jacobi SCC to
  max|Δμ| < 1e-13 e·Å (25 iterations).
- Coulomb prefactor k = 14.399645 eV·Å/e² (`units::coulomb`).

## Pinned values (hard-coded in `reg_pim-npt-polarizable-il.cpp`)

| Observable | Value | Units |
| --- | --- | --- |
| E_pol | −8.695749870555906 | eV |
| W_xx | −10.66184373345355 | eV |
| W_yy | −12.125914211727187 | eV |
| W_zz | −9.657852793291474 | eV |
| W_xy | −0.2325271422922529 | eV |
| W_xz | −2.2469237837238 | eV |
| W_yz | +1.375774440282161 | eV |
| P_pol = tr(W)/(3V) | −0.013021203851134545 | eV/ų |
| F0 | (+1.3910743837872235, −0.5347228671088384, +0.3580748823060393) | eV/Å |

Tolerances in the regression: |ΔE| ≤ 1e-8 eV, |ΔW| ≤ 1e-7 eV,
|ΔP_pol| ≤ 1e-10 eV/ų (covers std::sin vs np.sin position ulps, summation
order, and the 1e-13 SCC tail).

## NPT product leg (invariants, not pins)

Same lattice **stretched** to d = 2.9 Å (start tension ≈ −0.12 eV/ų; an
offline scan showed the compressed pin density d = 2.35 Å carries
P₀ ≈ +0.88 eV/ų, which would drive a ballistic MTK expansion through the
soft attractive side, while a stretched start contracts onto the stiff
r⁻¹² wall and self-limits) + `PairLJCut{0.01 eV, 3.0 Å, 4.0 Å}` +
`PairCoulCut{4.0 Å}` (virial-correct pair styles) + `PimInducedDipole`
under `NoseHooverNPT{T=800 K, τ_T=100 fs, P=0, τ_P=500 fs, dt=1 fs}`,
250 steps. Checks: SCC converged every step (component throws otherwise),
finite MTK conserved quantity, contraction happened (V_min < 0.97·V₀),
bounded (0.5·V₀ < V_min, V_max < 1.2·V₀), PIM virial live in the barostat
trace (leg-2 additivity).

Note (found gap, out of slug scope): `PairBorn` / `PairBuck` do not write
`kStress`, so a Fumi–Tosi NPT twin would barostat blind to the Born term —
which is why the product leg uses the LJ + coul/cut virial-correct stack.
Reported upstream in the worker summary for a queue row.

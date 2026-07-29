# spec: rm1-water-mopac-energy

**status:** code_complete
**branch:** `wf/rm1-water-mopac-energy`
**loop:** tb-method-loop, wave 4
**oracle:** OpenMOPAC (NDDO family lock — DFTB+ is forbidden here)

## Summary

`regressions/semiempirical/rm1_water_energy/` only checks that RM1 differs from
AM1 by more than 0.05 eV; there is no absolute OpenMOPAC pin anywhere on the
RM1 path. Deliver the first absolute RM1 gate for water: H₂O relaxed on the RM1
surface, ΔHf compared against a published OpenMOPAC/RM1 quantity, plus the
offline MOPAC input twin under `tools/tb_oracle/inputs/mopac/rm1-water/`.

## Domain basis

MOPAC's heat of formation is

    ΔHf = (E_elec + E_nuc − Σ_A EISOL(A)) · (eV→kcal/mol) + Σ_A ΔHf_exp(A)

`EISOL(A)` is the **free-atom** NDDO energy in the atom's ground-state
configuration. MOPAC does not tabulate it — `src/models/calpar.F90` **derives**
it from the active parameter set:

```fortran
gssc = max(ios-1,0);  k = iop;  gspc = ios*k;  l = min(k, 6-k)
gp2c = (k*(k-1))/2 + 0.5*(l*(l-1))/2
gppc = -0.5*(l*(l-1))/2
hspc = -k*ios*0.5
hpp  = max(0.1, 0.5*(gpp-gp2));  hsp = max(1e-7, hsp)
eisol = uss*ios + upp*iop + udd*iod
      + gss*gssc + gpp*gppc + gsp*gspc + gp2*gp2c + hsp*hspc
```
and in the same loop the multipole additives

```fortran
dd = (2n+1)(4 zs zp)^(n+1/2) / (zs+zp)^(2n+2) / sqrt(3)
qq = sqrt((4n²+6n+2)/20) / zp
am = gss/ev
ad : secant solve of  hsp/ev = 0.5 d − 0.5/sqrt(4 dd² + 1/d²)
aq : secant solve of  hpp/ev = 0.25 q − 0.5/sqrt(4 qq² + 1/q²)
                              + 0.25/sqrt(8 qq² + 1/q²)     (start (16 hpp/(48 ev qq⁴))^0.2)
```
with `ev = 27.211386245988` (MOPAC 2018-CODATA `fpcref(1,4)`).

**Every one of these is a function of the method's OWN parameters.** RM1 changes
`gss/gsp/gpp/gp2/hsp` as well as `uss/upp`, so an RM1 `eisol`/`aq` cannot be
obtained by perturbing the AM1 value.

## Finding (iron law — bug in the surface being pinned)

`src/cpu/semiempirical/rm1_parameters.cpp` builds `eisol` as *"the AM1 free-atom
energy shifted by the one-electron parameter change (uss/upp)"* — the
two-electron part is left at AM1 values. Verified numerically: the AM1 table
reproduces `calpar` for **all** its elements to 1e-9, while RM1 is off by

| Z | eisol stored | eisol per calpar | error [eV] | ΔHf error per atom [kcal/mol] |
|---|---|---|---|---|
| 6 C | −119.7956994 | −117.8673442 | +1.928 | −44.5 |
| 7 N | −202.8194020 | −205.0876425 | −2.268 | +52.3 |
| 8 O | −312.8526806 | −312.0403541 | +0.812 | −18.7 |
| 15 P | −125.0659220 | −123.1797789 | +1.886 | −43.5 |
| 16 S | −179.9286042 | −185.3861381 | −5.458 | +125.9 |

H is right by construction (`ios=1`, `iop=0` ⇒ no two-electron term).

Second defect in the same table: `aq` does not invert to RM1's own
`hpp = (gpp−gp2)/2` (implied/own ratio 0.79…0.95), while AM1's does to 1e-4.
`am`/`ad` additionally used `ev = 27.21` (MOPAC's `OLDCON` value) where the AM1
table and `nddo_two_electron.cpp` use 27.211386245988.

Consequence: ATV RM1 ΔHf(H₂O) = **−38.82** kcal/mol against an experimental
−57.80 — a 19 kcal/mol error, impossible for a model whose published AUE is
5.0 kcal/mol. AM1 on the same code path gives −59.26 (literature −59.24), so
the NDDO engine is sound and only the RM1 table is wrong.

## Design

1. **Fix `rm1_parameters.cpp`** — replace the derived columns (`am`, `ad`, `aq`,
   `eisol`) of every RM1 row with the `calpar.F90` values computed from that
   row's own primary parameters at `ev = 27.211386245988`. `dd`/`qq` were
   already exact and are unchanged. Primary parameters (uss, upp, betas, betap,
   zs, zp, alp, gss, gsp, gpp, gp2, hsp, Gaussian core terms) are untouched —
   they are the published RM1 values.
2. **Unit gate** `tests/cpu/semiempirical/test_nddo_calpar_derived.cpp` — re-derive
   `eisol`/`am`/`ad`/`aq`/`dd`/`qq` from `calpar.F90` for every AM1 **and** RM1
   row and require agreement. This is the invariant that was silently violated;
   it now fails the build if any future table row is hand-shifted again.
3. **Regression** `regressions/semiempirical/rm1-water-mopac-energy/` — H₂O
   relaxed on the RM1 surface in the two symmetry-unique internal coordinates
   (r_OH, ∠HOH) by a plain golden-section/pattern search over `RM1::run`
   (the regression is an out-of-library consumer; MINIMIZER is not built in the
   SE-only configuration and is not needed for a 2-parameter search), then ΔHf
   compared against the published OpenMOPAC/RM1 quantity. Geometry convention
   is stated explicitly in `REFERENCE.md`: a MOPAC ΔHf is quoted **at the
   method's own optimized geometry**, so the case must relax before comparing.
4. **Offline twin** `tools/tb_oracle/inputs/mopac/rm1-water/` — the exact MOPAC
   deck (`rm1-water.mop`) + `README.md` marking it a never-executed twin, and
   one additive row in `tools/tb_oracle/README.md`.

## Files

- `src/cpu/semiempirical/rm1_parameters.cpp` (fix + provenance comment)
- `include/atomiverse/cpu/semiempirical/am1_parameters.h` (eV→kcal/mol constant)
- `tests/cpu/semiempirical/test_nddo_calpar_derived.cpp` (new unit gate)
- `tests/CMakeLists.txt` (register)
- `regressions/semiempirical/rm1-water-mopac-energy/reg_rm1-water-mopac-energy.cpp`
- `regressions/semiempirical/rm1-water-mopac-energy/REFERENCE.md`
- `regressions/CMakeLists.txt` (`atv_add_regression(semiempirical rm1-water-mopac-energy)`)
- `tools/tb_oracle/inputs/mopac/rm1-water/{rm1-water.mop,README.md}`
- `tools/tb_oracle/README.md` (one additive row)

## Tasks

1. Recompute the RM1 derived columns from `calpar.F90`; patch the table.
2. Add the AM1+RM1 `calpar` invariant unit test; confirm it goes RED on the old
   table and GREEN after the fix.
3. Add the regression + REFERENCE.md with the citation, tolerance and Gaps.
4. Seed the offline MOPAC twin.
5. Re-run the semiempirical unit suite and the affected regressions.

## Testing

- `ctest -R '^test_cpu_(am1|rm1|dftb0|dftb2|dftb3|xtb_gfn0|nddo_)'`
- `cmake --build build --target run_reg_rm1_water_energy run_reg_rm1-water-mopac-energy`
- `run_reg_am1_water_force` (AM1 must be untouched by an RM1-table change)

## Outcome

Pin found and verified (both sources fetched and read directly, not taken on
trust):

| Quantity | Published RM1 | ATV | Δ |
|---|---|---|---|
| ΔHf(H₂O) | −57.86 kcal/mol (OpenMOPAC PARAM `RM1.heats` L546 @ 9859905) | −57.8553 | 0.0047 |
| ΔHf(H₂O) | −242.04 kJ/mol = −57.8489 (Scientia Plena 12(6), 069903 (2016), Tabela 1) | −57.8553 | 0.0064 |
| r(OH) | 0.956 Å (`RM1.geos` L165) | 0.95637 | 0.00037 |
| ∠HOH | 103.610° (`RM1.geos` L166) | 103.599 | 0.011 |

Calibration against a distributed MOPAC transcript (`openmopac/mopac`
`tests/keywords/RM1.out`, MOPAC2016 21.298M, CO₂): MOPAC −81.09054 / −81.55398
at d(CO) = 1.16 / 1.173109199 Å, ATV −81.09085 / −81.55429 — a constant
3.1e-4 kcal/mol at both geometries.

Extra fix beyond the plan: `kAm1EvToKcalMol` was the legacy `OLDCON` 23.061
instead of MOPAC's default 23.060547830619029 — that constant was the whole
0.0076 kcal/mol residual against the CO₂ transcript before the change.

Also verified: the AM1 water control in the same Scientia Plena table
(−247.90 kJ/mol = −59.2495 kcal/mol) matches ATV's AM1 −59.2507 to 1.2e-3,
and the RM1 CO₂ entry (−341.20 kJ/mol = −81.5487) matches ATV's −81.55435.

## Out of scope

- PM3/PM6/PM7 tables (no enum rows yet).
- Analytic RM1 gradients (still the fixed-density FD path).
- The AM1 table itself — it already satisfies the `calpar` invariant.

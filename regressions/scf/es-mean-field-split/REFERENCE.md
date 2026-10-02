# Reference — mean-field split on H₂/STO-3G and H₂O/6-31G

Out-of-library consumer. SCF `run` / `forces` stay on `RHF`. On H₂/STO-3G,
Seeger–Pople stability, CIS/TDA, and singlet EOM-EE-CCSD are not members of
`RHF` (`Stability<RHF>::run`, `TDA<RHF>::run`, `CCSD<RHF>::eom_ee_ccsd` —
PySCF `scf.stability.rhf_stability`, `tdscf.rhf.TDA`, `CCSD.eomee_ccsd_singlet`).
MP2, the CPHF nuclear Hessian, C-PCM, and static polarizability stay beside
the method as well.

## Oracle

PySCF **2.14.0**, `tools/pyscf_oracle/pin_scf.py` (and `pyscf.hessian.rhf` /
finite-field `dip_moment` for the two observables that script does not
emit). `cart=True`. Not an Atomiverse measurement.

```bash
python tools/pyscf_oracle/pin_scf.py rhf --geom h2 --basis sto-3g --grad --conv-tol 1e-12
python tools/pyscf_oracle/pin_scf.py tda --geom h2 --basis sto-3g --nstates 1 --conv-tol 1e-12
python tools/pyscf_oracle/pin_scf.py eom-ee-ccsd --geom h2 --basis sto-3g --nstates 1 --conv-tol 1e-12
python tools/pyscf_oracle/pin_scf.py mp2 --geom h2 --basis sto-3g --conv-tol 1e-12
python tools/pyscf_oracle/pin_scf.py rhf --geom water --basis 6-31g --grad --conv-tol 1e-12
python tools/pyscf_oracle/pin_scf.py rhf --geom water --basis 6-31g --pcm --conv-tol 1e-10
```

H₂ Hessian: `pyscf.hessian.RHF(mf).kernel()` at the same geometry, converted
with Ha→eV = 27.211386245988 and Bohr→Å = 1/1.8897261245650618
(eV/Å² = Ha/Bohr² × Ha→eV / (Bohr→Å)²).

H₂O polarizability: central difference of `RHF.dip_moment(unit='A.U.')`
under `hcore += F·r` (`int1e_r`), step 2×10⁻⁴ a.u., same geometry. The
STO-3G twin of this finite-field recipe reproduces the analytic CPHF pin
in `static-polariz-h2o` to < 1×10⁻⁵ a.u.

## Geometry

| System | Basis | Coordinates [Å] |
| --- | --- | --- |
| H₂ | STO-3G | H (0,0,0), H (0,0,0.74) |
| H₂O | 6-31G | O (0,0,0.1173), H (0, ±0.7572, −0.4692) |

20 Å cubic box. C-PCM: ε = 78.4, Bondi × 1.2, x = 0, Lebedev order 29.

## Pins (PySCF 2.14.0)

| Observable | Value | Gate (existing) |
| --- | --- | --- |
| H₂ E_RHF | −1.1167593073964255 Ha | 1×10⁻⁸ Ha |
| H₂ ∇E (atom 0, z) | −0.02767960070679465 Ha/Bohr | 1×10⁻⁶ Ha/Bohr |
| H₂ ∇E (atom 1, z) | +0.02767960070679487 Ha/Bohr | 1×10⁻⁶ Ha/Bohr |
| H₂ λ_min ¹(A+B) | +1.1296173364405577 Ha | 1×10⁻⁶ Ha |
| H₂ λ_min ³(A+B) | +0.40477548837977023 Ha | 1×10⁻⁶ Ha |
| H₂ TDA S1 | 25.807465778338745 eV | 1×10⁻⁴ eV |
| H₂ TDA f(S1) | 1.0950117175889624 | 1×10⁻³ |
| H₂ EOM-EE S1 | 26.365966612792256 eV | 1×10⁻⁴ eV |
| H₂ E_CCSD | −1.1372838344971299 Ha | 1×10⁻⁶ Ha |
| H₂ E_MP2 | −1.1298973809859585 Ha | 1×10⁻⁶ Ha |
| H₂ E_MP2 corr | −0.013138073589533008 Ha | 1×10⁻⁷ Ha |
| H₂ H_xx (0,0) | 1.923435534077 eV/Å² | 1×10⁻⁵ eV/Å² |
| H₂ H_zz (0,0) | 47.077486422795 eV/Å² | 1×10⁻⁵ eV/Å² |
| H₂O E_RHF | −75.98397447272193 Ha | 1×10⁻⁸ Ha |
| H₂O ∇E (O, z) | +0.02381868594901171 Ha/Bohr | 1×10⁻⁶ Ha/Bohr |
| H₂O ∇E (H, y) | ∓0.00441388034973933 Ha/Bohr | 1×10⁻⁶ Ha/Bohr |
| H₂O ∇E (H, z) | −0.011909342974506965 Ha/Bohr | 1×10⁻⁶ Ha/Bohr |
| H₂O ΔG_el (C-PCM) | −0.01614471877573211 Ha | 0.1 kcal/mol |
| H₂O α_xx, α_yy, α_zz | 1.394213713033, 6.648784943839, 4.410031211301 a.u. | 1×10⁻⁵ a.u. |

The full H₂ Hessian is the 6×6 atom-major matrix printed by
`hessian.RHF.kernel` (off-diagonal blocks = −diagonal; transverse xx = yy).

H₂ stability eigenvalues are the lowest eigenvalues of ¹(A+B) and ³(A+B)
built in PySCF from `mol.intor('int2e')` at R = 0.74 Å (same formulas as
`h2_rhf_stability`). `mf.stability(external=True, return_status=True)` at
this geometry reports internal-stable and external-stable. PySCF's Davidson
hop values are not the pin — they carry a different scale.

H₂O forces are `grad_ha_bohr` from `pin_scf.py rhf --grad`. The full 3×3
atom-major matrix is in the regression; transverse components are numerical
zeros (~1×10⁻¹⁶).

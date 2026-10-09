# Complex spin-orbital Hartree–Fock

## Systems

- H₂, bond 0.74 Å, STO-3G, charge 0
- H₂O, r(OH) = 0.95 Å, angle 104.5°, STO-3G, charge 0
  (oxygen at the origin, hydrogens in the xz plane)

## Oracle

PySCF 2.14.0. X2C is `scf.GHF(mol).x2c1e()` with
`with_x2c.xuncontract = False` and `approx='1E'` (the default). The
non-relativistic row is `scf.GHF(mol)` (real for these closed shells;
it matches `scf.RHF`).

```python
from pyscf import gto, scf
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g', verbose=0)
mf = scf.GHF(mol).x2c1e()
mf.with_x2c.xuncontract = False
e = mf.run(conv_tol=1e-12).e_tot
```

H₂O uses the same geometry as `h2o_rhf_sto3g`.

| System | Method | E [Ha] |
|---|---|---|
| H₂ | GHF x2c1e | −1.11677137468639 |
| H₂O | GHF x2c1e | −75.00458277099555 |
| H₂O | GHF (non-relativistic) | −74.96175405132141 |

## Acceptance

`|E − pin| ≤ 1e-8` Ha and the SCF reports converged.

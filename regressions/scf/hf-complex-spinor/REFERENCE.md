# Complex spinor Hartree–Fock

## Systems

Dirac–Hartree–Fock (PySCF `scf.DHF` docstring geometry):

- H₂, bond 1.0 Å (`H 0 0 0; H 0 0 1`), STO-3G and cc-pVDZ, charge 0

Complex two-component GHF (kept):

- H₂, bond 0.74 Å, STO-3G
- H₂O, r(OH) = 0.95 Å, angle 104.5°, STO-3G

## Oracle

PySCF 2.14.0, offline `tools/pyscf_oracle/pin_scf.py`. DHF is
`scf.DHF` with `conv_tol=1e-12`. X2C is `scf.GHF(mol).x2c1e()` with
`with_x2c.xuncontract = False`.

```bash
python tools/pyscf_oracle/pin_scf.py dhf --geom h2-1A --basis sto-3g
python tools/pyscf_oracle/pin_scf.py dhf --geom h2-1A --basis cc-pvdz
```

| System | Method | E [Ha] |
|---|---|---|
| H₂ 1.0 Å STO-3G | DHF | −1.0661226588590456 |
| H₂ 1.0 Å cc-pVDZ | DHF | −1.1001626190803968 |
| H₂ 0.74 Å STO-3G | GHF x2c1e | −1.11677137468639 |
| H₂O STO-3G | GHF x2c1e | −75.00458277099555 |
| H₂O STO-3G | GHF | −74.96175405132141 |

## Acceptance

`|E − pin| ≤ 1e-8` Ha and the SCF reports converged.

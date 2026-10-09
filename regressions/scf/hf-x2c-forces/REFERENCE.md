# H2 / STO-3G and H2O / 6-31G RHF x2c-1e analytic gradient

Classic molecular systems with a nonzero nuclear force. An isolated atom
(Ar / def2-TZVP) is an energy pin only: symmetry makes the force zero.

Oracle (PySCF 2.14.0), from the offline tool:

```
python tools/pyscf_oracle/pin_scf.py rhf --geom h2 --basis sto-3g --grad --x2c --conv-tol 1e-12
python tools/pyscf_oracle/pin_scf.py rhf --geom water --basis 6-31g --grad --x2c --conv-tol 1e-12
```

That is `scf.RHF(mol).x2c1e()` with `with_x2c.xuncontract = False` and
`nuc_grad_method().kernel()` (`pyscf.x2c.sfx2c1e_grad.hcore_grad_generator`,
`gen_sf_hfw`). Cartesian bases (`cart=True`). Geometries are the oracle
presets (Å):

- H2: H(0, 0, 0), H(0, 0, 0.74)
- H2O: O(0, 0, 0.1173), H(0, ±0.7572, −0.4692)

| System | Quantity | Pin |
| --- | --- | --- |
| H2 STO-3G | E_x2c | −1.1167713747002186 Ha |
| H2 STO-3G | g(H1)_z | −0.02767894667735915 Ha/Bohr |
| H2 STO-3G | g(H2)_z | +0.027678946677359373 Ha/Bohr |
| H2O 6-31G | E_x2c | −76.03879282908161 Ha |
| H2O 6-31G | g(O)_z | +0.02372121511333436 Ha/Bohr |
| H2O 6-31G | g(H+y)_y, g(H+y)_z | −0.004279885030643804, −0.011860607556359426 Ha/Bohr |
| H2O 6-31G | g(H−y)_y, g(H−y)_z | +0.004279885031505337, −0.011860607556974934 Ha/Bohr |

Transverse components are numerical noise (≤ 1e-12) and are pinned as 0.
Both molecules have |g| ≈ 0.02 Ha/Bohr, so a vanished force fails the gate.

The H2 x2c−nonrelativistic force shift is ~6.5e-7 Ha/Bohr, below the force
tolerance; the energy shift is −1.206e-5 Ha, so a non-relativistic
Hamiltonian fails the H2 energy gate. On H2O the g(O)_z shift is ~9.7e-5
Ha/Bohr, so a non-relativistic gradient fails the force gate.

Acceptance, same gates as the HCl x2c gradient: |E − pin| ≤ 1e-6 Ha;
element-wise |g − pin| ≤ 1e-6 Ha/Bohr; |Σ_A g_A| ≤ 1e-8 Ha/Bohr.

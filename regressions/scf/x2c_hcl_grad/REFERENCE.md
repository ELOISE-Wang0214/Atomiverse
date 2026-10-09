# HCl / 6-31G RHF x2c-1e analytic gradient

Geometry (Å): Cl at the origin, H at (0, 0, 1.35) — stretched past the
equilibrium bond so the force is nonzero.

Oracle (PySCF 2.14.0), from the offline tool:

```
python tools/pyscf_oracle/pin_scf.py rhf --geom hcl --basis 6-31g --grad --x2c --conv-tol 1e-12
```

That is `scf.RHF(mol).x2c1e()` with `with_x2c.xuncontract = False` and
`nuc_grad_method().kernel()` (`pyscf.x2c.sfx2c1e_grad.hcore_grad_generator`).
Cartesian 6-31G (`cart=True`).

| Quantity | Pin |
| --- | --- |
| E_x2c | −461.49265438490323 Ha |
| g(Cl)_z | −0.028498915583206497 Ha/Bohr |
| g(H)_z | +0.028498915583214046 Ha/Bohr |

Transverse components are numerical noise (~1e-12) because the molecule
lies on z. The x2c shift of g_z relative to the non-relativistic gradient
is about 3.5e-4 Ha/Bohr, so a non-relativistic force fails the gate.

Acceptance: |E − pin| ≤ 1e-6 Ha; element-wise |g − pin| ≤ 1e-6 Ha/Bohr;
|Σ_A g_A| ≤ 1e-8 Ha/Bohr.

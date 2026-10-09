# UHF / STO-3G — H₂⁺, H₂, H₂O⁺

Out-of-library consumer of `cpu::scf::UHF` and, when CUDA is on,
`cuda::scf::UHF`. Same geometry and the same PySCF numbers on both.

## Oracle

PySCF 2.14.0, offline:

```bash
python tools/pyscf_oracle/pin_scf.py uhf --geom h2 --basis sto-3g --charge 1 --spin 1 --grad --conv-tol 1e-12
python tools/pyscf_oracle/pin_scf.py uhf --geom h2 --basis sto-3g --charge 0 --spin 0 --grad --conv-tol 1e-12
python tools/pyscf_oracle/pin_scf.py uhf --geom water --basis sto-3g --charge 1 --spin 1 --grad --conv-tol 1e-12
```

`oracle: pyscf`. Cartesian STO-3G. Energies are `e_tot_ha`. Gradients are
`grad_ha_bohr` from `nuc_grad_method().kernel()` (∂E/∂R).

| System | Electrons (α, β) | Charge | E [Ha] |
| --- | --- | --- | --- |
| H₂⁺ at 0.74 Å | (1, 0) | +1 | −0.538205447565 |
| H₂ at 0.74 Å | (1, 1) | 0 | −1.116759307396 |
| H₂O⁺ (oracle water geometry) | (5, 4) | +1 | −74.655898713621 |

H₂⁺ gradient [Ha/Bohr], atoms along z: (0, 0, +0.191355066881) and
(0, 0, −0.191355066881).

H₂ gradient: (0, 0, −0.027679600707) and (0, 0, +0.027679600707).

H₂O⁺ gradient, atoms O, H(+y), H(−y):

| Atom | x | y | z |
| --- | --- | --- | --- |
| O | 0 | 0 | −0.095196667501 |
| H | 0 | −0.084875484950 | +0.047598333750 |
| H | 0 | +0.084875484950 | +0.047598333750 |

Gates: |ΔE| < 1e-8 Ha, element-wise |Δg| < 1e-6 Ha/Bohr, |Σ g| < 1e-8 Ha/Bohr.

The CUDA path is device s/p integrals and `launch_uhf_forces_sp`. It throws
if the basis is outside STO-3G l ≤ 1.

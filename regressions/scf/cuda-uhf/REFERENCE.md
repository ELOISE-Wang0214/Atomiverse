# CUDA UHF / STO-3G — O₂ and OH

Out-of-library consumer of `cpu::scf::UHF` and, when CUDA is on,
`cuda::scf::UHF`. Same geometry and the same PySCF numbers on both.
The device leg must report `last_scf_on_device` and `last_force_on_device`.

## Oracle

PySCF 2.14.0, offline (`oracle: pyscf`):

```bash
python tools/pyscf_oracle/pin_scf.py uhf --geom o2 --basis sto-3g --spin 2 --grad --conv-tol 1e-12
python tools/pyscf_oracle/pin_scf.py uhf --geom oh --basis sto-3g --spin 1 --grad --conv-tol 1e-12
```

Cartesian STO-3G. Energies are `e_tot_ha`. Gradients are `grad_ha_bohr`
from `nuc_grad_method().kernel()` (∂E/∂R). PySCF's default initial guess
is MINAO (`UHF.init_guess_by_minao`: spin-averaged, because `_break_dm_spin_symm`
does not run when `mol.spin != 0`). The core-Hamiltonian guess is a different
O₂ minimum (−147.378559141664 Ha) and is not the pin.

OH uses r(O–H) = 0.96966 Å along z (the experimental equilibrium distance).
O₂ uses the repository geometry r(O–O) = 1.2075 Å.

| System | Electrons (α, β) | Charge | Spin (2S) | E [Ha] |
| --- | --- | --- | --- | --- |
| O₂ triplet | (9, 7) | 0 | 2 | −147.633946785450 |
| OH doublet | (5, 4) | 0 | 1 | −74.362633276842 |

O₂ gradient [Ha/Bohr], atoms along z: (0, 0, +0.024814046022) and
(0, 0, −0.024814046022).

OH gradient [Ha/Bohr]: (0, 0, +0.056145481913) and (0, 0, −0.056145481913).

Gates: |ΔE| < 1e-8 Ha, element-wise |Δg| < 1e-6 Ha/Bohr, |Σ g| < 1e-8 Ha/Bohr.
Same thresholds as the STO-3G UHF gradient cases.

The CPU call sets `init_guess = MINAO`, no damping, and CDIIS from iteration 1
(PySCF `diis_start_cycle = 1`, `damp = 0`). The `UHF` constructor's historical
Hcore guess plus 0.7 damping is the other O₂ basin. The CUDA path packs the
same spin-averaged MINAO density once at setup; the SCF and the analytic
gradient stay on device (`launch_eri_sp`, `launch_uhf_forces_sp`).

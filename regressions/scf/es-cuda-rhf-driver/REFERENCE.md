# Reference — staged CUDA RHF driver, H₂/STO-3G and H₂O/6-31G

## Method

Restricted Hartree–Fock total energy through the public Frame API.

- CPU: `cpu::scf::RHF{basis, charge, tol, max_iter}.run(frame)`.
- CUDA: `cuda::scf::RHF` is `cuda::Driver<scf::Phase>` (fock, mix,
  diagonalize, density, check). The H₂/STO-3G leg must report
  `last_scf_on_device()` and `last_force_on_device()`. A host ERI fill
  under a CUDA banner is not a pass.

**One regression, multi-device.** `reg_es-cuda-rhf-driver --device=cpu|cuda|both`
(default `both` when built with CUDA). The `.cpp` is the CUDA-off path.

## Oracle

Offline PySCF 2.14.0, `tools/pyscf_oracle/pin_scf.py`, `conv_tol=1e-12`,
Cartesian (`cart=True`; neither basis has a d shell, so pure vs Cartesian
does not change the energy).

```
python tools/pyscf_oracle/pin_scf.py rhf --geom h2 --basis sto-3g
python tools/pyscf_oracle/pin_scf.py rhf --geom water --basis 6-31g
```

| System | Basis | Geometry [Å] | E_RHF [Ha] |
| --- | --- | --- | --- |
| H₂ | STO-3G | H (0,0,0), H (0,0,0.74) | −1.1167593073964255 |
| H₂O | 6-31G | O (0,0,0.1173), H (0, ±0.7572, −0.4692) | −75.98397447272204 |

`oracle: pyscf` 2.14.0 (`pyscf.scf.RHF.kernel`). Not an Atomiverse measurement.

## Tolerance

Same bands as `scf/h2o_rhf_sto3g_grad` for the RHF total energy:

| Device | \|E − E_PySCF\| |
| --- | --- |
| cpu | 1e-8 Ha |
| cuda (H₂/STO-3G only) | 1e-6 Ha |

H₂O/6-31G is outside the device s/p STO-3G integral kernel (`device_sp_sto3g_ok`).
The CUDA leg does **not** run it. Grading the transitional host-ERI fill as a
device result is forbidden. The 6-31G pin is the CPU leg only.

## What is checked

Both systems converge and `energy()` [Hartree] and `kEnergy` [eV] / Ha agree
with the pin inside the device band. The CUDA H₂ leg additionally requires
the device SCF and analytic-force flags. Non-convergence, a host-integral
CUDA path, or an out-of-tolerance energy exits non-zero.

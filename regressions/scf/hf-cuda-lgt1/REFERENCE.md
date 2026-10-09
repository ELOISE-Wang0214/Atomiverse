# Reference — H₂O RHF, Cartesian d (6-31G*) and Cartesian f (def2-TZVP)

**One case, multi-device:** `reg_hf-cuda-lgt1 --device=cpu|cuda|both`
(default `both` when built with CUDA). Same PySCF pins on every device.
Public `RHF` + Frame only.

## Systems

Geometry [Å], the `water` preset of `tools/pyscf_oracle/pin_scf.py`:

| Atom | x | y | z |
|---|---|---|---|
| O | 0 | 0 | 0.1173 |
| H | 0 | 0.7572 | −0.4692 |
| H | 0 | −0.7572 | −0.4692 |

20 Å cubic box (gas phase; the cell is not used by the isolated-molecule SCF).

| System | Basis | What l > 1 it exercises | PySCF `cart` | nao | E_RHF [Ha] |
|---|---|---|---|---|---|
| H₂O | 6-31G* | Cartesian d on O (6 functions) | true | 19 | −76.01050498825705 |
| H₂O | def2-TZVP | Cartesian f on O (10 functions), then the basis's spherical sandwich (7f / 5d) | false | 43 | −76.05900476057711 |

def2-TZVP is a spherical basis. The device (and the CPU) still **evaluates**
the oxygen f shell as ten Cartesian GTOs (`fill_cart_comps`, libcint order)
and the oxygen d shells as six, then applies the same `CartToSph` sandwich
as `IntegralWorkspace` (`integrals.cpp`). A Cartesian-only PySCF run
(`--cart`, nao 48) prints **−76.05941121596425 Ha** and is **not** this pin:
that number is the unsandwiched basis, which the public `BasisSet::Def2TZVP`
does not publish.

H₂ / STO-3G and H₂O / 6-31G never put l > 1 on an atom, so they are not
this case.

## Oracle

PySCF 2.14.0, `conv_tol = 1e-12`, offline:

```
python tools/pyscf_oracle/pin_scf.py rhf --geom water --basis 6-31g*
python tools/pyscf_oracle/pin_scf.py rhf --geom water --basis def2-tzvp
```

`cart` follows the script default (Pople Cartesian, def2 spherical).

CPU tolerance **1e-8 Ha**, the same window as H₂O RHF/def2-SVP
(`regressions/scf/h2o_rhf_def2svp_sph`). The CPU leg matches the pin to
~1e-12 Ha.

CUDA tolerance **1e-6 Ha**, the same window as the device H₂/STO-3G energy
pin. The device McMurchie–Davidson ERI is not bit-identical to the CPU Rys
engine; on this geometry the device is 3e-8 Ha (6-31G*) and 4e-7 Ha
(def2-TZVP) from the pin.

## Device leg

`cuda::scf::RHF` with `with_gradient(false)` (analytic forces outside
STO-3G l ≤ 1 still throw `launch_rhf_forces_lgt1`). The SCF must report
`last_scf_on_device()` and must not report `last_force_on_device()`.
There is no host ERI fill. l ≥ 4 throws `launch_eri_g`.

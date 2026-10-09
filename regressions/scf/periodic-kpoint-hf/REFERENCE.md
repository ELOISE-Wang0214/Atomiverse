# periodic-kpoint-hf

H₂ / STO-3G in a cubic 4 Bohr cell. The cell edge is smaller than the
Gaussian pair reach (`cell.rcut` ≈ 16.51 Bohr at `precision = 1e-8`), so
lattice images of the basis overlap the home cell.

Monkhorst–Pack meshes use PySCF `make_kpts` (`wrap_around=False`).
H–H lies on c, so `k_mesh = {2,1,1}` samples the bond-perpendicular
direction and `{1,1,2}` samples along the bond.

## Oracle

PySCF 2.14.0 `pyscf.pbc.scf`, default `FFTDF`, `exxdiv='ewald'`.

```
python tools/pyscf_oracle/pin_scf.py krhf
python tools/pyscf_oracle/pin_scf.py krhf --kmesh 2,1,1
python tools/pyscf_oracle/pin_scf.py kuhf --kmesh 2,1,1 --spin 2
```

| Case | Call | E [Ha] |
| --- | --- | --- |
| KRHF `{1,1,2}` | `KRHF` | −1.3541438798330374 |
| KRHF `{2,1,1}` | `KRHF` | −1.4471748215213558 |
| triplet KUHF `{2,1,1}` | `KUHF`, `cell.spin = 2` | −0.7658470471886634 |

`KUHF.nelec` (`pyscf/pbc/scf/kuhf.py`) is the mesh total
`(tot_electrons(nk) + spin)//2`. `tot_electrons(2) = 4` and `spin = 2`
give Nα = 3, Nβ = 1. Spin is not multiplied by the number of k-points.
The same call prints `enuc = -0.6281049380287604` Ha
(`pyscf.pbc.gto.cell.ewald`).

## Gate

`|E - E_pin| ≤ 1e-6` Ha. The lattice-sum cutoff on S and T
(PySCF `cell.rcut`) is the leading discrepancy versus libcint
`int1e_*`; J and K on a matched AO grid agree much more tightly.

## Not this case

Analytic nuclear forces (`krhf_nuclear_forces`, `kuhf_nuclear_forces`).
Shells with l > 0 (`krhf_ao_lgt0`). Device KUHF.

# periodic-kpoint-hf

H₂ / STO-3G restricted Hartree–Fock on a cubic 4 Bohr cell with a
Monkhorst–Pack mesh `k_mesh = {1,1,2}` (Γ and the zone-boundary point
along c; PySCF `make_kpts` default `wrap_around=False`).

## Oracle

PySCF 2.14.0 `pyscf.pbc.scf.KRHF`, default `FFTDF`, `exxdiv='ewald'`.

```
python tools/pyscf_oracle/pin_scf.py krhf
```

```
E_tot = -1.3541438798330374 Ha
```

The same call prints `enuc = -0.6281049380287604` Ha
(`pyscf.pbc.gto.cell.ewald`).

## Gate

`|E - E_pin| ≤ 1e-6` Ha. The lattice-sum cutoff on S and T
(PySCF `cell.rcut`) is the leading discrepancy versus libcint
`int1e_*`; J and K on a matched AO grid agree much more tightly.

## Not this case

Analytic nuclear forces (`krhf_nuclear_forces`). k-point UHF.
Shells with l > 0 (`krhf_ao_lgt0`).

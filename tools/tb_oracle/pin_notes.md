
## 2026-07-27 — dftb2-dftbplus-energy-pin (H2)

- Binary: DFTB+ **24.1** x86_64-linux release tarball (not built from 25.1; GF 14
  build of 25.1 hit a Fortran field mismatch on this host).
- Parameters: **mio-1-1** `H-H.skf` (dftbparams/mio).
- Geometry: H–H = 0.74 Å, cluster, SCC on, tol 1e-8.
- Total energy: −0.6749509324 Ha → −18.36635052 eV (×27.211386245988).
- ATV matches via transcribed static tables + DFTB+ “new” SK interpolation +
  spline repulsive (see `basis_data/dftb_mio_hh_sk.h`).
- Wall (DFTB+ cold): ~0.15 s; ATV measured in REFERENCE after build.

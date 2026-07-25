# periodic-gw-kpoints — Multi-k spectral G0W0 (k-diagonal product)

## Goal

Close the **multi-k residual** after Γ-point spectral G0W0
(`periodic-gw-si` ✓). `GpwDft::g0w0` must accept a converged multi-k
closed-shell RKS reference and return quasiparticle energies without
throwing.

## Science (product)

**k-diagonal spectral G0W0@KS** on the same RI-GPW / complex-Bloch
infrastructure as multi-k dRPA (`periodic-rpa-kpoints`):

- Per sampled \(k\): Coulomb-only Casida spectrum from **same-\(k\)**
  (vertical, \(q=0\)) OV pair densities → Furche reduction with \(X+Y\)
  amplitudes.
- Correlation self-energy: frequency-exact spectral sum over those poles
  (molecular / Γ `RHF::g0w0` / `GpwDft::g0w0` form with \(|t_{smp}|^2\)).
- Exact-exchange diagonals: truncated (MT) Coulomb on complex pair
  densities (G=0-safe foothold; periodic HF still blocked).
- Local \(v_{\mathrm{xc}}\): same pure GGA/LDA as the SCF on the
  BZ-weighted density; MO diagonals per \(k\).
- Virtual space for screening: `mp2_n_virt` (0 = all); QP levels for the
  full MO spectrum at each \(k\).

**Not in product this slug:** full inter-\(k\) / finite-\(q\) dielectric
\(W(\mathbf{q})\); frequency-integrated RI-GW (CP2K default); UKS/UGW;
EXX total-energy assembly.

## Public API

```cpp
cpu::gpw::GpwDft m;
m.k_mesh[0] = 2; m.k_mesh[1] = 1; m.k_mesh[2] = 1;
m.k_time_reversal = false;
m.mp2_ri = true;
// build(); run(atoms);
const cpu::gpw::G0W0Result& gw = m.g0w0(atoms);  // multi-k, no throw
// gw.n_kpts == 2; gw.gap_ev = min_k direct QP gap [eV]
```

`G0W0Result` gains `n_kpts` (0 = real Γ path; multi-k count otherwise).
`gap_ev` / `mf_gap_ev` are the **minimum direct** HOMO→LUMO gaps over the
sampled mesh. `qp_energy_ha` / `mf_energy_ha` / `sigma_x_ha` / `vxc_ha`
hold the spectrum at the \(k\) that realises the min direct QP gap
(diagnostic; same layout as Γ).

## Acceptance

See `acceptance.md`. Offline CP2K pin input under
`tools/cp2k_oracle/inputs/periodic-gw-kpoints.inp`. Regression
`regressions/gpw/periodic-gw-kpoints/` vs first-green band (spectral ≠
CP2K frequency-integrated RI-GW → no invented absolute CP2K QP digits).

## Deps

- `periodic-gw-si` ✓ (Γ spectral product)
- `periodic-rpa-kpoints` ✓ (multi-k k-diagonal dRPA infrastructure)

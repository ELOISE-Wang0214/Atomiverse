# Spec: `periodic-gw-si`

## Goal

Ship **Γ-point spectral G0W0** quasiparticle energies on a converged
**periodic GPW PBE** reference for artificial sc-Si₂ — the periodic GPW
counterpart of molecular `RHF::g0w0` (dRPA-screened spectral self-energy +
Newton QP solve), reporting HOMO/LUMO **QP gaps**.

Roadmap §10 **Periodic GW** ○ → product foothold. Multi-k GW / frequency-
quadrature RI-GW remain residuals. **Not** molecular G0W0 (PySCF loop).

## Public API (no Config bag, no factory)

```cpp
cpu::gpw::GpwDft m;
// ... members; build(); run(atoms);
const cpu::gpw::G0W0Result& gw = m.g0w0(atoms);
// gw.qp_energy_ha, gw.mf_energy_ha [Ha]; gw.gap_ev [eV]; gw.converged
```

Parameters reuse `mp2_n_virt` (virtual cap for the dRPA screening space)
and `mp2_ri` (RI-GPW vs Poisson foothold for neutral OV pair densities).
No new public Config type.

## Science

Closed-shell **spectral G0W0@KS** (PBE) at Γ:

1. Build `(ia|jb)` via product RI-GPW (default) or Poisson foothold —
   same kernels as MP2/RPA (`internal/mp2_kernels.h`).
2. Singlet Coulomb-only dRPA Casida (`α_hf = 0`): recover positive
   poles \(\Omega_s\) and amplitudes \(X,Y\) (Furche reduction —
   same as molecular `OrbitalHessian::direct_rpa`).
3. Transition densities
   \(\rho_{s,pq}=\sum_{ia}\sqrt{2}\,(X+Y)_{ia}^{s}\,(ia|pq)\)
   (neutral \(\rho_{ia}\); periodic Coulomb OK).
4. Correlation self-energy (spectral, real-axis, small \(\eta\)):
   \(\mathrm{Re}\,\Sigma_{c,p}(\omega)\) from RPA poles (Berkelbach /
   PySCF `GWExact` form used by `RHF::g0w0`).
5. Exact-exchange diagonal \(\Sigma_{x,p}=-\sum_i(pi|pi)\) via
   **truncated (MT) Coulomb** (G=0-safe Γ foothold; periodic HF/X
   remains blocked as a product path).
6. Local KS \(v_{\mathrm{xc},p}=\langle\psi_p|V_{\mathrm{xc}}|\psi_p\rangle\)
   from the same PBE functional / density as the SCF.
7. Newton solve of
   \(\omega=\varepsilon_p+\mathrm{Re}\,\Sigma_c(\omega)+\Sigma_{x,p}-v_{\mathrm{xc},p}\)
   for every MO; gap from HOMO/LUMO QP levels [eV on the public result].

References: Hybertsen & Louie, PRB **34**, 5390 (1986); Furche, JCP
**114**, 5982 (2001); molecular twin `reg_g0w0-h2-sto3g`; GPW dRPA
`reg_periodic-rpa-si`.

## System (oracle deck)

| Item | Value |
|------|--------|
| Crystal | sc-Si₂ (same as `periodic-rpa-si` / `periodic-ri-mp2-si`) |
| Cell | cubic a = 5.43 Å, periodic, Γ only |
| Basis / PP | DZVP-MOLOPT-GTH / GTH-PBE-q4 |
| XC | PBE |
| Cutoff | 30 Ry |
| SCF | tol 1e-6, mix 0.3, level_shift 0.3 Ha |

## Out of scope

- Multi-k / frequency-quadrature RI-GW (CP2K default machinery)
- BSE@GW, self-consistent GW, vertex corrections
- Periodic HF / EXX product path (MT exchange is Γ foothold only)
- CUDA
- Molecular Gaussian G0W0 (PySCF loop owns that)

## Acceptance

See `acceptance.md`.

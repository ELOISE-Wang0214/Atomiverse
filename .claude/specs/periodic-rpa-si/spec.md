# Spec: `periodic-rpa-si`

## Goal

Ship **Γ-point direct RPA (dRPA)** correlation energy on a converged
**periodic GPW PBE** reference for artificial sc-Si₂ — the periodic
counterpart of molecular `RHF::direct_rpa` (Furche plasmon formula), using
the same occupied–virtual MO repulsion tensor path as
`GpwDft::mp2` / `periodic-ri-mp2-si`.

Roadmap §10 **Periodic RPA** ○ → product foothold. Multi-k RPA and
periodic GW remain residuals (`periodic-gw-si` held).

## Public API (no Config bag, no factory)

```cpp
cpu::gpw::GpwDft m;
// ... members; build(); run(atoms);
const cpu::gpw::RpaResult& r = m.direct_rpa(atoms);
// r.e_corr, r.e_scf, r.e_total [Ha]; r.omega_ha; r.ri
```

Parameters reuse `mp2_n_virt` (virtual cap) and `mp2_ri` (RI-GPW vs
Poisson foothold for `(ia|jb)`). No new public Config type.

## Science

Closed-shell **direct RPA** on KS (PBE) orbitals:

1. Build `(ia|jb)` via product RI-GPW (default) or Poisson foothold —
   same kernels as MP2 (`internal/mp2_kernels.h`).
2. Singlet Coulomb-only Casida matrices (`α_hf = 0`):
   \(A_{ia,jb}=\delta_{ij}\delta_{ab}(\varepsilon_a-\varepsilon_i)+2(ia|jb)\),
   \(B_{ia,jb}=2(ia|jb)\) (real orbitals).
3. Hermitian reduction \((A-B)^{1/2}(A+B)(A-B)^{1/2}\); positive roots
   \(\Omega_s\).
4. Furche plasmon:
   \(E_c=\tfrac12\bigl(\sum_s\Omega_s-\mathrm{Tr}\,A\bigr)\).
5. \(E_{\mathrm{total}}=E_{\mathrm{SCF}}+E_c\) (dRPA@PBE — no EXX
   replacement; periodic HF remains blocked).

References: Furche, JCP **114**, 5982 (2001); Ren et al., NJP **14**,
053020 (2012); molecular twin `reg_direct-rpa-h2-sto3g`.

## System (oracle deck)

| Item | Value |
|------|--------|
| Crystal | sc-Si₂ (same as `periodic-ri-mp2-si`) |
| Cell | cubic a = 5.43 Å, periodic, Γ only |
| Basis / PP | DZVP-MOLOPT-GTH / GTH-PBE-q4 |
| XC | PBE |
| Cutoff | 30 Ry |
| SCF | tol 1e-6, mix 0.3, level_shift 0.3 Ha |

## Out of scope

- Multi-k RPA, GW, BSE@RPA
- EXX+RPA total energies (needs periodic HF)
- CUDA
- Frequency-quadrature / Gaussian-RI_AUX path (CP2K default machinery);
  pin is plasmon E_c on GPW integrals; CP2K input retained for offline
  correlation re-pin notes

## Acceptance

See `acceptance.md`.

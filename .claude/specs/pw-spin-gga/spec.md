# Spec — `pw-spin-gga`

## Goal

Productize **true collinear spin-GGA** on the plane-wave FFT grid, closing the
residual after:

- closed-shell PW true GGA PBE ✓ (`reg_pw-pbe-gga-si`)
- collinear spin-polarized PW (σ=0 LSDA-limit PBE per channel) ✓
  (`reg_pw-spin-polarized-o2`, `reg_pw-uks-fe-smearing`, multi-k UKS)
- GPW collinear spin-GGA ✓ (`xc_grid_gga_spin` / UKS-PBE)

Today `XCBuilder::step_spin` calls `evaluate_xc(PBE, ρ↑, ρ↓)` with default
`σ=0` — the **LDA limit of PBE**, not true spin-GGA. Density gradients never
enter the spin XC energy or potential.

## Science

Collinear UKS XC energy

\[
E_{\mathrm{xc}}[\rho_\uparrow,\rho_\downarrow]
= \int e\big(\rho_\uparrow,\rho_\downarrow,
  \sigma_{\uparrow\uparrow},\sigma_{\uparrow\downarrow},\sigma_{\downarrow\downarrow}\big)\,d^3r
\]

with spectral gradients on the density FFT grid

\[
\nabla\rho_\sigma(\mathbf{r})
= \mathcal{F}^{-1}\big[i\mathbf{G}\,\rho_\sigma(\mathbf{G})\big],
\quad
\sigma_{\sigma\sigma'}=(\nabla\rho_\sigma)\cdot(\nabla\rho_{\sigma'}).
\]

Spin potentials (libxc / GPW convention):

\[
v_\sigma
= v_{\rho_\sigma}
- \nabla\cdot\big(
  2\,v_{\sigma_{\sigma\sigma}}\,\nabla\rho_\sigma
  + v_{\sigma_{\uparrow\downarrow}}\,\nabla\rho_{\bar\sigma}
\big).
\]

Same collapse as GPW `xc_grid_gga_spin`. Closed-shell recovery: when
\(\rho_\uparrow=\rho_\downarrow=\rho/2\), \(v_\uparrow=v_\downarrow\) must match
`XCBuilder::step` (restricted true GGA).

`XCType::LSDA` keeps the pure local spin path (no gradient). GGA selectors
(`PBE`, …) use the true spin-GGA path. No public `*Config` / factory.

## System (regression)

| Item | Value |
|------|--------|
| Molecule | O₂ triplet (M = 3 → N↑=7, N↓=5 for 12 valence e⁻) |
| Geometry | R(O–O) = 1.2075 Å |
| Cell | cubic L = 8.0 Å |
| PP | SG15 ONCV PBE O `tests/data/pw/O_ONCV_PBE-1.2.upf` |
| E_cut | 12 Ry |
| k-mesh | Γ only |
| XC | `XCType::PBE` true collinear spin-GGA |
| SCF | tol 5×10⁻⁴, max_iter 60, mix_beta 0.3 |

## API

```cpp
cpu::pw::PW pw;
pw.upf_path = ".../O_ONCV_PBE-1.2.upf";
pw.ecut_ry = 12.0;
pw.k_mesh = 1;
pw.spin_polarized = true;
pw.multiplicity = 3;
pw.xc = cpu::scf::XCType::PBE; // true spin-GGA (not σ=0)
pw.build(frame);
const auto& r = pw.run(frame);
// r.magnetization == 2; r.energy_ha finite & pinned
```

Public Frame energy remains eV via evaluators; `PWResult` reports Ha/Ry as today.

## Out of scope

- Noncollinear / spin-orbit
- CUDA PW spin-GGA
- Spin-GGA forces/stress productization (follow-on if needed)
- Live CP2K/QE in ctest (offline pins only)
- Molecular post-HF, MetalWalls, GAPW/PAW

## Deps

`pw-pbe-gga-si` ✓, `pw-spin-polarized-o2` ✓, GPW `xc_grid_gga_spin` ✓

# gpw-kpoint-uks-stress-gga

## Goal

Expand multi-k UKS analytic cell stress (`gpw_stress` / `energy_fixed_P_k`)
from `XcKind` ∈ {PADE, PBE, PBEsol} to pure spin-GGA
**RevPBE / RPBE / BLYP** (same spin XC kernels as multi-k UKS SCF and
Γ-UKS stress). Foothold `gpw-kpoint-uks-stress` (PADE) is **not** the
GGA product — this slug closes the residual.

## Scope

**In**

- `energy_fixed_P_k` UKS branch: accept RevPBE / RPBE / BLYP / OptB88VdW;
  dispatch `xc_grid_*_spin` like Γ-UKS `energy_fixed_P`.
- Unit: dual-gate multi-k UKS analytic vs multi-SCF FD for each GGA.
- Regression `regressions/gpw/gpw-kpoint-uks-stress-gga/`: high-spin H₂,
  MP 1×1×2, dual-gate for RevPBE / RPBE / BLYP + activation vs PBE.
- Offline CP2K documentation twin under `tools/cp2k_oracle/inputs/…`
  (no live CP2K in ctest).

**Out**

- Meta-GGA, hybrid, GAPW/PAW.
- Full OptB88 product regression (kernel allowed; product pin is GGA trio).
- Molecular post-HF / MetalWalls / CUDA.

## Public API

No new symbols. Callers set `GpwDft::xc = XcKind::RevPBE|RPBE|BLYP` with
`n_alpha`/`n_beta` + `k_mesh` and call `gpw_stress` / `GPWEvaluator`
`compute_stress`. No `*Config`, no `make_*` / `build_*`.

## Science

Hellmann–Feynman fixed-\(\{P_k^\sigma,W_k^\sigma\}\) Lagrangian under
orthorhombic diagonal strain (Nielsen–Martin / Quickstep). Spin XC on
collocated \(\rho_\alpha,\rho_\beta\) via pure GGA spin kernels. Units:
internal Ha / Ha·Bohr⁻³; public Frame energy eV, stress eV/Å³ when exposed.

## Deps

`gpw-kpoint-uks-stress` ✓, `gpw-uks-revpbe` ✓, `gpw-uks-rpbe` ✓,
`gpw-uks-blyp` ✓, `gpw-uks-stress` ✓ (Γ GGA UKS stress).

## Acceptance

See `acceptance.md`.

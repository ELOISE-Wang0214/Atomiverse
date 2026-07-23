# Spec: `gpw-pbesol-si`

**Slug:** `gpw-pbesol-si`  
**Loop:** CP2K capability (GPW / periodic)  
**Branch:** `wf/gpw-pbesol-si`  
**Chapter:** gpw  
**Status:** open  

## Goal

Ship **PBEsol GGA** on the GPW path for a textbook silicon system
(Γ-point closed-shell, DZVP-MOLOPT-GTH + GTH-PBE-q4), with offline CP2K-style
pinning of \(E_{\mathrm{tot}}\) [Ha]. Extends the GGA ladder after PBE pins
(`gpw-water-pbe-mt`, `pw-pbe-gga-si`).

## Public API (no Config bags, no factories)

| Surface | Shape |
| --- | --- |
| `cpu::scf::XCType::PBEsol` | Pure GGA; `is_gga` true; hybrid fraction 0 |
| `cpu::gpw::XcKind::PBEsol` | GPW grid XC kind (RKS) |
| `cpu::gpw::GPWEvaluator` / `GpwDft` | `xc = XcKind::PBEsol` member / ctor arg |
| `GaussianBasis(Si, "DZVP-MOLOPT-GTH")` | lookup constructor (table row) |
| `GthPotential(Si, "GTH-PBE-q4")` | lookup constructor (table row) |

Parameters remain constructor args / public members. No `*Config`, no
`make_*` / `build_*`.

## Science

- **PBEsol** (Perdew et al., *Phys. Rev. Lett.* **100**, 136406 (2008)):
  same analytic form as PBE with solid-state parameters
  \(\mu = 10/81\), \(\kappa = 0.804\) (exchange) and
  \(\beta = 0.046\) (correlation, libxc / CP2K `XC_FUNCTIONAL PBEsol`
  convention); PW92 LDA correlation base unchanged.
- Shared pointwise engine: `cpu::scf::evaluate_xc(XCType::PBEsol, …)`.
- GPW grid path reuses the closed-shell GGA assembly (spectral \(\nabla\rho\),
  divergence of \(v_\sigma\) flux) already used for PBE.

## System (regression)

- Si₂ diamond fragment in cubic orthorhombic cell \(a = 5.43\) Å
  (GPW is orthorhombic-only; Γ-point foothold — not full 8-atom conventional
  bulk, which exceeds the ≤ 60 s e2e budget with DZVP-MOLOPT).
- Positions [Å]: \((0,0,0)\), \((a/4,a/4,a/4)\).
- Basis: **DZVP-MOLOPT-GTH**; PP: **GTH-PBE-q4** (zion 4); \(n_{\mathrm{occ}}=4\).
- Cutoff: 300 Ry; mix 0.4; tol \(10^{-6}\); max_iter 80; periodic.
- Observable: \(E_{\mathrm{tot}}\) [Ha] from `last_energy_hartree()` / `GpwResult::e_total`.

## Tasks

1. Add `XCType::PBEsol` + parameterized PBE-family exchange/correlation
   (\(\mu,\beta\)) in `xc_functional.{h,cpp}`; metadata (`is_gga`, hybrids).
2. Add `XcKind::PBEsol`; route GPW SCF + analytic-force XC grid through
   `evaluate_xc(PBEsol, …)` (RKS only; UKS GGA remains out of scope).
3. Transcribe Si **GTH-PBE-q4** + **DZVP-MOLOPT-GTH** into
   `basis_data/gpw_tables.h` (CP2K data files, published sources).
4. Unit tests: PBEsol@\(\sigma=0\) ≡ PBE@\(\sigma=0\); PBEsol≠PBE at finite
   \(\sigma\); Si table lookup; GGA grid LDA-limit for PBEsol.
5. Regression `regressions/gpw/gpw-pbesol-si/` with REFERENCE.md offline pin
   + dual-gate (PBEsol vs PBE energy gap).
6. CP2K twin input under `tools/cp2k_oracle/inputs/` (never launched by ctest).
7. Perf gate: fair wall if CP2K binary present; else `perf_exempt` + reason.

## Out of scope

- UKS + GGA PBEsol; k-point GPW; PW-only productization (PW inherits
  `XCType::PBEsol` via shared `evaluate_xc` if callers set it).
- Molecular post-HF (PySCF loop).
- Live CP2K in ctest.

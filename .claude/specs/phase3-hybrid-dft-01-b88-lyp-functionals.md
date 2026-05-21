---
title: B88 / LYP / VWN-RPA semi-local XC functionals for the native engine
status: approved
created: 2026-05-21
---

# B88 / LYP / VWN-RPA semi-local XC functionals for the native engine

## Summary
The native exchange-correlation engine (`src/cpu/scf/xc_functional.cpp`) implements Slater, VWN5, PW92, PBE-X and PBE-C, but lacks the three semi-local functionals that the upcoming B3LYP hybrid (link 2 of the `phase3-hybrid-dft` chain) is built from: Becke-1988 GGA exchange, Lee-Yang-Parr GGA correlation, and the RPA-parametrized Vosko-Wilk-Nusair LDA correlation. This link adds exactly those three as native, spin-polarized building blocks — each returning energy density plus the potential derivatives (`vrho`, and for the GGA pair `vsigma`) the downstream hybrid Fock/gradient code consumes — and validates every output channel against the libxc oracle to 1e-12. No user-visible behavior changes in this link: the functionals are internal building blocks, not new `XCType` enum values; link 2 wires them into B3LYP/PBE0.

## Domain basis
All formulas use atomic units (Hartree energy, Bohr length). Each functional is spin-polarized: energy density depends on `rho_alpha`, `rho_beta` (and, for GGAs, `sigma_xx = grad rho_x . grad rho_y`); each returns `vrho` and, for GGAs, `vsigma`.

- **B88 exchange** — libxc `XC_GGA_X_B88` (id 106). Becke 1988 GGA exchange. The exchange enhancement is applied per spin channel over the spin-scaled Slater (Dirac) LDA exchange:
  `E_x = sum_sigma integral rho_sigma * eps_x^LDA(rho_sigma) * F(x_sigma) dtau`,
  with reduced gradient `x_sigma = |grad rho_sigma| / rho_sigma^{4/3}` and Becke enhancement
  `F(x) = 1 + (beta/C_x) * x^2 / (1 + 6*beta*x*arcsinh(x))`, `beta = 0.0042`, `C_x` the Dirac exchange constant. Reference: A. D. Becke, *Phys. Rev. A* **38**, 3098 (1988). DOI: 10.1103/PhysRevA.38.3098.
- **LYP correlation** — libxc `XC_GGA_C_LYP` (id 131). Lee-Yang-Parr 1988 GGA correlation. libxc implements the Miehlich-Savin-Stoll-Preuss reformulation, which is the closed analytic form in terms of `rho_alpha`, `rho_beta`, `sigma_aa`, `sigma_ab`, `sigma_bb`. Empirical constants `a = 0.04918`, `b = 0.132`, `c = 0.2533`, `d = 0.349`. References: C. Lee, W. Yang, R. G. Parr, *Phys. Rev. B* **37**, 785 (1988), DOI: 10.1103/PhysRevB.37.785; B. Miehlich, A. Savin, H. Stoll, H. Preuss, *Chem. Phys. Lett.* **157**, 200 (1989), DOI: 10.1016/0009-2614(89)87234-3.
- **VWN-RPA correlation** — libxc `XC_LDA_C_VWN_RPA` (id 8). The RPA-parametrized Vosko-Wilk-Nusair LDA correlation. It uses the **same** Pade-approximant functional form as the existing `vwn5_correlation` (the `vwn5_F` / `vwn5_dFdx` machinery), with a **different** parameter set `(A, x0, b, c)` for the paramagnetic, ferromagnetic, and spin-stiffness channels. Reference: S. H. Vosko, L. Wilk, M. Nusair, *Can. J. Phys.* **58**, 1200 (1980), DOI: 10.1139/p80-159 (RPA columns).
  **Naming trap:** PySCF's alias `VWN3` maps to libxc id 8 (`VWN_RPA`); libxc's own `XC_LDA_C_VWN_3` is id 30, a different functional. B3LYP (libxc id 402) uses id 8. This link implements id 8 and the symbol MUST be named `vwn_rpa_correlation` — never `vwn3`.

## Design
Three new static functions are added to the native XC engine, each mirroring the existing per-functional pattern (`slater_exchange`, `pbe_exchange_spin`, `pbe_correlation`, `vwn5_correlation`): spin-polarized signatures that take densities (and sigmas for GGAs) and return energy density plus potential derivatives by out-parameter.

- `vwn_rpa_correlation(rho_a, rho_b, ec, vc_a, vc_b)` — reuses the existing `VWN5Params` struct, `vwn5_F`, `vwn5_dFdx`, `spin_f`, `spin_df` machinery verbatim; only three new `constexpr VWN5Params` (RPA paramagnetic / ferromagnetic / spin-stiffness) are introduced. This is a re-parametrization, not new math.
- `b88_exchange_spin(rho_s, sigma_ss, ex_s, vrho_s, vsigma_ss)` — per-spin GGA exchange in the same shape as `pbe_exchange_spin`, including the `rho_s < kRhoMin/2` early return that yields exactly zero for an empty spin channel (matching libxc).
- `lyp_correlation(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, ec, vc_a, vc_b, vsigma_aa, vsigma_ab, vsigma_bb)` — full spin-polarized closed-shell + open-shell GGA correlation in the Miehlich form.

All three are closed-shell-correct as the unpolarized limit of the spin-polarized expression (no separate RKS path), consistent with how the engine already handles RKS via `rho_alpha = rho_beta`.

**File-size split (required).** `xc_functional.cpp` is 757 lines; `CLAUDE.md` caps source files at 800. Adding B88 + LYP would cross that cap. Per `coding-style.md` ("many small files"), the two new GGA functionals are placed in a new sibling translation unit `src/cpu/scf/xc_gga.cpp`, with internal declarations in a new private header `include/atomiverse/cpu/scf/xc_internal.h` (not facade-public, consistent with the `cpu/scf/internal/` precedent). `xc_gga.cpp` holds `b88_exchange_spin` and `lyp_correlation`; `vwn_rpa_correlation` stays in `xc_functional.cpp` next to `vwn5_correlation` since it shares that file's VWN machinery and adds < 30 lines. `xc_internal.h` declares the three new functions plus the shared `VWN5Params` / `vwn5_F` helpers needed across the two translation units. The public `xc_functional.h` and the `XCType` enum are **unchanged** in this link.

## Files to create or modify
- `include/atomiverse/cpu/scf/xc_internal.h` (new) — private header declaring `b88_exchange_spin`, `lyp_correlation`, `vwn_rpa_correlation`, and the shared `VWN5Params` / `vwn5_F` / `vwn5_dFdx` helpers.
- `src/cpu/scf/xc_functional.cpp` — add `vwn_rpa_correlation` + RPA parameter sets; move shared VWN helpers behind `xc_internal.h`.
- `src/cpu/scf/xc_gga.cpp` (new) — `b88_exchange_spin` and `lyp_correlation`.
- `tests/cpu/scf/oracle/xc_libxc_oracle.h` — declare single-functional oracle helpers `evaluate_gga_x_b88`, `evaluate_gga_c_lyp`, `evaluate_lda_c_vwn_rpa`.
- `tests/cpu/scf/oracle/xc_libxc_oracle.cpp` — implement those helpers (`XC_GGA_X_B88`, `XC_GGA_C_LYP`, `XC_LDA_C_VWN_RPA`).
- `tests/cpu/scf/test_xc_vs_libxc.cpp` — add B88 / LYP / VWN-RPA parity batteries over the existing density + sigma sweep.
- `tests/CMakeLists.txt` — register `src/cpu/scf/xc_gga.cpp` into the `atomiverse` build sources (or confirm glob coverage) and ensure the new test cases compile under `ATV_BUILD_RHF` + `ATV_BUILD_XC_ORACLE`.

## Tasks
- [ ] Write failing oracle helpers `evaluate_gga_x_b88` / `evaluate_gga_c_lyp` / `evaluate_lda_c_vwn_rpa` and their declarations (`tests/cpu/scf/oracle/xc_libxc_oracle.h`, `tests/cpu/scf/oracle/xc_libxc_oracle.cpp`)
- [ ] Write failing VWN-RPA parity test over the density battery incl. spin-polarized points (`tests/cpu/scf/test_xc_vs_libxc.cpp`)
- [ ] Create private header `include/atomiverse/cpu/scf/xc_internal.h` exposing shared VWN helpers + the three new function declarations
- [ ] Implement `vwn_rpa_correlation` with RPA parameter sets in `src/cpu/scf/xc_functional.cpp`
- [ ] Write failing B88 + LYP parity tests over the density + sigma sweep incl. spin-polarized points (`tests/cpu/scf/test_xc_vs_libxc.cpp`)
- [ ] Implement `b88_exchange_spin` and `lyp_correlation` in new file `src/cpu/scf/xc_gga.cpp` and register it in `tests/CMakeLists.txt` build sources
- [ ] Add Doxygen `@brief` / `@param` (with units) / `@return` per doc.style for all three functionals and the new header
- [ ] Verify all three functionals match libxc (ids 106 / 131 / 8) to 1e-12 on every channel across the full battery
- [ ] Run full check + test suite

## Testing strategy
- **Happy path:** B88 / LYP / VWN-RPA each match their libxc oracle (`XC_GGA_X_B88`, `XC_GGA_C_LYP`, `XC_LDA_C_VWN_RPA`) on `exc`, `vrho_alpha`, `vrho_beta` — and for B88/LYP additionally `vsigma_aa`, `vsigma_ab`, `vsigma_bb` — to the combined 1e-12 abs/rel tolerance, over the existing `kPoints` density battery (unpolarized, partially polarized, highly polarized, low density, high density) and, for the GGAs, the `kSigmaValues` sweep.
- **Edge cases:** spin-polarized / highly-polarized points (e.g. `1.0/0.05`, `zeta ~ 0.91`) exercise the `(1-zeta)` near-singular regime — must match libxc, matching the established battery's treatment; zero-sigma grid points for the GGAs; low-density point near the `kRhoMin` cutoff returns zero from both native and oracle.
- **Domain validation:** the 1e-12 libxc match IS the domain validation — libxc is the reference implementation of all three functional IDs. Additionally a textbook spot-check that the unpolarized closed-shell limit (`rho_alpha = rho_beta`) of each spin-polarized expression is correct (B88 reduces to its restricted form; VWN-RPA `zeta=0` reduces to the paramagnetic Pade term).
- **Regression:** the existing `XCvsLibxc.LsdaBattery` / `PbeBattery` / `DensityCutoff` cases and the full `ctest` suite stay green — no shared code path (Slater, VWN5, PW92, PBE) is altered in behavior by the file split.

## Out of scope
- The hybrid functionals B3LYP and PBE0 themselves — link 2 of this chain.
- Hybrid exact-exchange Fock construction and hybrid analytic gradients — link 3.
- Range-separated / long-range-corrected functionals (CAM-B3LYP, wB97X).
- Exposing the new functionals as `XCType` enum values or via `evaluate_xc` / `evaluate_xc_batch` dispatch — deliberately deferred; link 2 introduces the hybrid `XCType` entries that consume these building blocks.

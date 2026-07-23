# Spec: `gpw-dftu-nio`

**Slug:** `gpw-dftu-nio`  
**Branch:** `wf/gpw-dftu-nio`  
**Queue:** CP2K capability loop — Periodic DFT+U (GPW)  
**Status:** open  
**Created:** 2026-07-23  
**Resume:** incomplete queue row (`docs/developer/cp2k-method-queue.md` `## in_progress`)

## Goal

Ship **periodic Γ-point UKS GPW + Dudarev DFT+U** on textbook **AFM NiO**,
pinned offline against **CP2K Quickstep PLUS_U**. Roadmap `docs/roadmap/dft.md`
marks DFT+U / Periodic DFT+U as ○.

## Scope (in)

1. **Dudarev DFT+U** on `cpu::gpw::GpwDft` / `GPWEvaluator` (Γ, UKS + RKS):
   - Per-kind (element) `U_eff = U − J` with angular momentum `L` (Ni `L=2`).
   - Occupation matrix from **Mulliken** projector: \(n_{mm'} = (PS)_{mm'}\) on
     the AO subspace of atom \(A\) with \(\ell = L\) (symmetrized).
   - Energy (per spin channel σ):  
     \(E_U = \frac{U_{\mathrm{eff}}}{2}\sum_{A\sigma}\big[\mathrm{Tr}\,n^{A\sigma}
     - \mathrm{Tr}(n^{A\sigma}n^{A\sigma})\big]\) [Ha; \(U_{\mathrm{eff}}\) converted
     from eV].
   - Potential added to \(H^\sigma\):  
     \(V^{A\sigma}_{mm'} = U_{\mathrm{eff}}\big(\tfrac12\delta_{mm'} - n^{A\sigma}_{mm'}\big)\)
     on the same AO subspace.
2. **Public API** (no `*Config` bag, no `make_*` / `build_*`):
   - `struct PlusU { molrs::Element element; int l; double u_minus_j_ev; };`
   - `GpwDft::plus_u` / `GPWEvaluator::plus_u` as `std::vector<PlusU>`.
   - Optional AFM seed: `std::vector<double> atom_spin_guess` (per-atom, units-free
     sign/magnitude hint applied only to the initial density guess).
   - Result: `GpwResult::e_dftu` [Ha]; `GpwResult::mulliken_spin` [e] per atom.
3. **Ni data tables** (in-tree static, no parsers):
   - `GthPotential(Ni, "GTH-PADE-q18")` and `GaussianBasis(Ni, "SZV-GTH-PADE-q18")`.
   - `GthPotential(O, "GTH-PADE-q6")` and `GaussianBasis(O, "SZV-GTH-PADE")` if missing.
4. **Unit tests** under `tests/cpu/gpw/` (no e2e): table lookup; DFT+U energy/potential
   on a tiny synthetic occupation; empty `plus_u` leaves legacy path unchanged.
5. **One regression** `regressions/gpw/gpw-dftu-nio/`:
   - AFM NiO (2 formula units), Γ UKS PADE + PLUS_U on Ni \(d\) (`L=2`),
   - Assert \(E_{\mathrm{tot}}\) [Ha] and Ni local Mulliken spin signs/magnitudes
     vs offline CP2K pin in `REFERENCE.md`.
6. **Offline CP2K pin** via `tools/cp2k_oracle/inputs/gpw-dftu-nio/`;
   ctest never launches CP2K.

## Scope (out)

- Liechtenstein full \(U,J\) tensor; orbital-dependent U ramping productization.
- k-point DFT+U; UKS PBE-GGA; analytic DFT+U forces/stress.
- GAPW / PAW projectors; MPI scaling.

## Public API (binding shape)

```cpp
cpu::gpw::Kind Ni{molrs::Element::Ni,
                  cpu::gpw::GaussianBasis(molrs::Element::Ni, "SZV-GTH-PADE-q18"),
                  cpu::gpw::GthPotential(molrs::Element::Ni, "GTH-PADE-q18")};
cpu::gpw::Kind O{molrs::Element::O,
                 cpu::gpw::GaussianBasis(molrs::Element::O, "SZV-GTH-PADE"),
                 cpu::gpw::GthPotential(molrs::Element::O, "GTH-PADE-q6")};

cpu::gpw::GPWEvaluator eval({Ni, O});
eval.cutoff_ry = 300.0;
eval.xc = cpu::gpw::XcKind::PADE;
eval.n_alpha = 24;
eval.n_beta  = 24;
eval.plus_u.push_back({molrs::Element::Ni, /*l=*/2, /*u_minus_j_ev=*/6.0});
// atom_spin_guess size == n_atoms; AFM seed e.g. {+1, 0, -1, 0}
eval.atom_spin_guess = {+1.0, 0.0, -1.0, 0.0};
cpu::single_point(frame, eval);
// eval.last_result().e_total [Ha], .e_dftu [Ha], .mulliken_spin
```

## Tasks

- [x] T1: Spec + CP2K twin input on `wf/gpw-dftu-nio`.
- [x] T2: Offline CP2K pin → `REFERENCE.md` (E_tot Ha, mag, settings, version).
- [x] T3: Ni/O PADE GTH + SZV tables in `gpw_tables.h` + unit lookup tests.
- [x] T4: `internal/dft_plus_u.h` Dudarev E/V + AO subspace helper.
- [x] T5: Wire into UKS SCF; `GpwResult::{e_dftu,mulliken_spin}`; spin guess.
- [x] T6: `GPWEvaluator` members + wiring.
- [x] T7: Unit tests (`test_cpu_gpw_dftu.cpp`).
- [x] T8: Regression `reg_gpw-dftu-nio` + CMake `atv_add_regression(gpw …)`.
- [x] T9: Build `ATV_BUILD_GPW=ON`; regression PASS; perf in REFERENCE.

## Acceptance criteria

| # | Criterion | Verify |
|---|-----------|--------|
| A1 | `plus_u` empty: existing GPW UKS/RKS unit tests still pass | `ctest -R gpw` |
| A2 | Ni/O PADE table lookup constructs without throw | unit |
| A3 | AFM NiO DFT+U SCF converges | regression |
| A4 | \|E_ATV − E_CP2K\| ≤ **5e-2 Ha** (first-green SZV/single-grid vs OT multigrid twin) | regression |
| A5 | Ni Mulliken spin \(m_{\mathrm{Ni}} > 0.3\) and within **0.5 e** of CP2K pin | regression |
| A6 | No public `*Config` / `make_*` / `build_*` for this feature | review / scrub |
| A7 | Offline pin documented: CP2K version, input path, geometry, U, L | REFERENCE.md |
| A8 | Perf: \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} ≤ 10\) **or** `perf_exempt` with reason | REFERENCE.md |

**`acceptance_met=true` only when A1–A8 hold.**

## Dual gate

1. **Correctness:** A3–A5 vs offline CP2K pin (no live CP2K in ctest).
2. **Speed:** fair wall when CP2K binary available; else `perf_exempt=true`.

## Geometry / settings (pin twin)

| Item | Value |
|------|--------|
| System | **NiO molecule** (solid AFM Γ-only does not converge on SZV in CP2K or ATV first green; molecule is the PLUS_U foothold) |
| Bond | R(Ni–O) = **1.668 Å** along *z*, centered in cell |
| Coordinates (Å) | Ni (5, 5, 4.166), O (5, 5, 5.834) |
| Cell | **10 × 10 × 10 Å**, `PERIODIC NONE`, MT Poisson |
| Basis | Ni `SZV-GTH-PADE-q18`, O `SZV-GTH-PADE` |
| PP | Ni `GTH-PADE-q18`, O `GTH-PADE-q6` |
| XC | PADE (LDA LSD) |
| DFT+U | Ni \(L=2\), \(U-J = 6.0\) eV, Mulliken |
| Cutoff | 250 Ry (Atomiverse single-grid; CP2K twin NGRIDS 4 / REL_CUTOFF 40) |
| Spin | UKS multiplicity 3 → **Nα=13, Nβ=11** (24 valence e⁻) |
| CP2K SCF | OT (FULL_SINGLE_INVERSE / DIIS) + OUTER_SCF (converges; diag+Broyden does not on this twin) |

Solid rocksalt AFM NiO remains a follow-on once k-points / better bases land.

## Dependencies

- Existing UKS GPW path (`n_alpha`/`n_beta`, PADE LSD) — `gpw-uks-o2-gamma` ✓.
- FFTW3 / `ATV_BUILD_GPW=ON`.

## Notes

- \(U_{\mathrm{eff}}\) is **eV on the public surface** (CP2K `U_MINUS_J [eV]`);
  converted with `units::kHartreeInEv` internally.
- First green uses SZV + single density grid; tighten tol after multigrid parity.

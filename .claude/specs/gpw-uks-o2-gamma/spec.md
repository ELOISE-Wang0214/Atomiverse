# Spec: `gpw-uks-o2-gamma`

**Slug:** `gpw-uks-o2-gamma`  
**Branch:** `wf/gpw-uks-o2-gamma`  
**Queue:** CP2K capability loop — Γ-point UKS GPW  
**Status:** open  
**Created:** 2026-07-23

## Goal

Ship **Γ-point unrestricted Kohn–Sham (UKS) GPW** for open-shell systems.
Today `GpwDft` is closed-shell only (`n_occ` doubly occupied). Product need:
O₂ triplet energy and ⟨S²⟩ vs offline **CP2K Quickstep LSD** pin.

## Scope (in)

1. **UKS path in `cpu::gpw::GpwDft`** (Γ only):
   - Public members `n_alpha`, `n_beta` (occupation counts per spin).
   - When `n_alpha + n_beta > 0`, run unrestricted SCF; otherwise keep
     closed-shell `n_occ` path bit-identical.
   - Spin densities ρ_α, ρ_β → collocate separately; Hartree from ρ_α+ρ_β;
     spin-resolved PADE LSD XC (`pade_lsd`); separate H_α / H_β diagonalization.
   - Result fields: total energy [Ha], ⟨S²⟩ (Löwdin), converged flag; expose
     spin densities / MO sets as needed for unit tests (α/β C,P).
2. **Public evaluator surface** (`GPWEvaluator`): optional `n_alpha` / `n_beta`
   members; if both zero, closed-shell from valence/2 (unchanged).
3. **XC:** PADE (LDA LSD) required for this pin. PBE GGA UKS is **out of
   scope** (throw or document as unsupported when `xc == PBE` and UKS).
4. **Unit tests** under `tests/cpu/gpw/` (no e2e): UKS H₂ or O₂-tiny converges;
   closed-shell path still passes; S² formula unit-checked if extracted.
5. **One regression** `regressions/gpw/gpw-uks-o2-gamma/`:
   - O₂, R = 1.2075 Å, cubic cell 12 Å, `PERIODIC NONE` + MT Poisson twin,
   - DZVP-MOLOPT-GTH + GTH-PBE-q6, CUTOFF 400 Ry, PADE, multiplicity 3
     (Nα=7, Nβ=5 valence),
   - Assert E_tot [Ha] and ⟨S²⟩ vs offline CP2K pin in `REFERENCE.md`.
6. **Offline CP2K pin** via `tools/cp2k_oracle/inputs/gpw-uks-o2-gamma/`;
   ctest never launches CP2K.

## Scope (out)

- k-point UKS, ROKS GPW, analytic UKS forces, PBE-GGA spin UKS, MPI scaling.
- Molecular all-electron UKS (already in SCF / pyscf loop).

## Public API (binding shape)

```cpp
cpu::gpw::GpwDft gpw;
// … cell, cutoff_ry, kinds via atoms …
gpw.n_alpha = 8;
gpw.n_beta  = 4;
gpw.xc = cpu::gpw::XcKind::PADE;
gpw.nonperiodic = true;
gpw.build(atoms);
auto& r = gpw.run(atoms);
// r.e_total [Ha], r.s2, r.converged

// Evaluator (Frame I/O, energy → eV on frame.meta):
cpu::gpw::GPWEvaluator eval({Kind{Element::O, basis, pot}});
eval.n_alpha = 8;
eval.n_beta = 4;
eval.nonperiodic = true;
eval.xc = cpu::gpw::XcKind::PADE;
// no *Config bag, no make_* factory
```

## Tasks

- [x] T1: Spec + CP2K twin input committed on `wf/gpw-uks-o2-gamma`.
- [x] T2: Offline CP2K 2026.1 pin → `REFERENCE.md` (E_tot Ha, ⟨S²⟩, settings).
- [x] T3: Spin-resolved `xc_grid_pade_lsd(rho_a, rho_b)` (+ wire into solver).
- [x] T4: UKS SCF loop in `GpwDft` (`n_alpha`/`n_beta`, dual H, S²).
- [x] T5: `GPWEvaluator` members + occupation wiring.
- [x] T6: Unit test(s) for UKS path (fast, no full textbook e2e).
- [x] T7: Regression `reg_gpw-uks-o2-gamma` + CMake `atv_add_regression(gpw …)`.
- [x] T8: Build with `ATV_BUILD_GPW=ON`; regression PASS; perf note in REFERENCE.

## Acceptance criteria

| # | Criterion | Verify |
|---|-----------|--------|
| A1 | Closed-shell GPW (`n_occ` only) still converges H₂ unit tests | `ctest -R gpw` |
| A2 | UKS mode: O₂ (8,4) PADE/DZVP-MOLOPT-GTH/MT converges | regression / unit |
| A3 | \|E_ATV − E_CP2K\| ≤ **1e-3 Ha** (first-green; tighten later) | regression |
| A4 | ⟨S²⟩ within **0.15** of CP2K pin (UKS contamination) | regression |
| A5 | No public `*Config` / `make_*` / `build_*` for this feature | review / scrub |
| A6 | Offline pin documented: CP2K version, input path, geometry | REFERENCE.md |
| A7 | Perf: \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} ≤ 10\) **or** `perf_exempt` with reason | REFERENCE.md |
| A8 | Frame energy path (evaluator) writes eV consistent with Ha×conversion | optional unit / reg |

**`acceptance_met=true` only when A1–A8 hold.**

## Dual gate

1. **Correctness:** A3–A4 vs offline CP2K pin (no live CP2K in ctest).
2. **Speed:** fair wall when CP2K binary available; else `perf_exempt=true`.

## Geometry / settings (pin twin)

| Item | Value |
|------|--------|
| Molecule | O₂, R = 1.2075 Å along z |
| Cell | 12 × 12 × 12 Å, `PERIODIC NONE`, MT Poisson |
| Basis | DZVP-MOLOPT-GTH |
| PP | GTH-PBE-q6 |
| XC | PADE (LDA LSD) |
| Cutoff | 400 Ry |
| Spin | LSD, multiplicity 3 → Nα=7, Nβ=5 (12 valence e⁻) |
| Γ-only | yes |

## Dependencies

- In-tree GTH/MOLOPT tables (`GthPotential(O,"GTH-PBE-q6")`,
  `GaussianBasis(O,"DZVP-MOLOPT-GTH")`).
- `cpu::xc::pade_lsd` (already spin-resolved).
- FFTW3 / `ATV_BUILD_GPW=ON`.

## Notes

- Roadmap: `docs/roadmap/dft.md` Γ-point periodic UKS ○ → this is the
  **isolated Γ-UKS foothold** (MT), not k-point.
- Forces/stress UKS deferred; do not extend analytic force path unless free.

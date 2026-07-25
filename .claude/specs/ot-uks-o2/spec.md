# Spec: `ot-uks-o2`

**Slug:** `ot-uks-o2`  
**Branch:** `wf/ot-uks-o2`  
**Queue:** CP2K capability loop — OT for Γ-UKS  
**Status:** acceptance_met (dual gate green)  
**Created:** 2026-07-26

## Goal

Ship the **Orbital Transformation (OT)** SCF minimizer for **Γ-point UKS**
GPW. RKS OT is product (`use_ot`; `reg_ot-minimizer-water`); UKS currently
throws `"OT path is not implemented for unrestricted SCF"`. Close that gap
on O₂ triplet vs offline CP2K OT pin.

## Scope (in)

1. **Γ-UKS OT path** in `cpu::gpw::GpwDft::run` when `is_uks() && use_ot`:
   - Separate occupied manifolds \(C_\alpha\) (\(n_\alpha\)) and \(C_\beta\)
     (\(n_\beta\)), each S-orthonormal.
   - Densities \(P_\sigma = C_\sigma C_\sigma^\top\) (occupation 1).
   - Spin-resolved KS build (same as UKS DIIS: collocate α/β, Hartree from
     total ρ, spin XC, dual \(H_\alpha/H_\beta\)).
   - Projected open-shell gradient per spin:
     \(g_\sigma = (I - S C_\sigma C_\sigma^\top)(2 H_\sigma C_\sigma)\).
   - FULL_SINGLE_INVERSE preconditioner per spin; joint Polak–Ribière CG
     step + Armijo line search on the rebuilt total energy (same pattern
     as closed-shell OT).
   - Convergence: \(\|g\|_\infty < k_{\mathrm{OT}}\) and \(|\Delta E| < \mathrm{tol}\).
   - On success: `store_converged_state_uks` for canonical MOs + Löwdin ⟨S²⟩.
2. **Public API:** no new flags — existing `use_ot` + `n_alpha`/`n_beta`.
3. **Unit test** (fast): UKS OT vs DIIS energy parity on toy H₂ high-spin
   (or closed-shell still OK); no full O₂ e2e in `tests/`.
4. **Regression** `regressions/gpw/ot-uks-o2/`:
   - Same O₂ geometry as `gpw-uks-o2-gamma` (centered, R=1.2075 Å, 12 Å MT).
   - PADE, DZVP-MOLOPT-GTH, GTH-PBE-q6, Nα=7, Nβ=5, cutoff 250–400 Ry
     (budget ≤60 s; document).
   - Assert E vs CP2K OT pin + OT/DIIS parity in-tree.
5. **Offline CP2K pin** `tools/cp2k_oracle/inputs/ot-uks-o2/`; never in ctest.

## Scope (out)

- Multi-k OT, OT + DFT+U, OT forces/stress, ROKS OT, CPMD/XL with UKS OT.
- Changing closed-shell OT or DIIS UKS behaviour beyond removing the throw.

## Public API (binding)

```cpp
cpu::gpw::GPWEvaluator eval({Kind{Element::O, basis, pot}});
eval.n_alpha = 7;
eval.n_beta = 5;
eval.use_ot = true;   // product OT; no *Config, no make_*
eval.xc = cpu::gpw::XcKind::PADE;
eval.nonperiodic = true;
cpu::single_point(frame, eval);  // energy eV on Frame; Ha on last_energy_hartree()
```

## Tasks

- [x] T1: Spec + acceptance on `wf/ot-uks-o2`.
- [x] T2: CP2K OT twin input + offline pin → `REFERENCE.md`.
- [x] T3: Open-shell OT gradient (occ=1) + UKS OT SCF loop.
- [x] T4: Unit test OT UKS vs DIIS parity (fast).
- [x] T5: Regression `reg_ot-uks-o2` + CMake registration.
- [x] T6: Dual gate (correctness + fair wall; ratio ≈ 8.5).

## Dual gate

1. **Correctness:** \|E_ATV,OT − E_CP2K,OT\| ≤ pin tol; \|E_OT − E_DIIS\| ≤ 1e-4 Ha.
2. **Speed:** \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} ≤ 10\) or `perf_exempt` + reason.

## Geometry / settings

| Item | Value |
|------|--------|
| Molecule | O₂ triplet, centered in 12 Å cube |
| R(O–O) | 1.2075 Å |
| Basis / PP | DZVP-MOLOPT-GTH / GTH-PBE-q6 |
| XC | PADE LSD |
| Spin | Nα=7, Nβ=5 (multiplicity 3) |
| Poisson | MT, `PERIODIC NONE` |
| OT (CP2K) | `MINIMIZER DIIS`, `PRECONDITIONER FULL_SINGLE_INVERSE` |

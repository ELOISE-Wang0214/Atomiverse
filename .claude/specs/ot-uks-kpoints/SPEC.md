# Spec — `ot-uks-kpoints`

| Field | Value |
| ----- | ----- |
| Slug | `ot-uks-kpoints` |
| Loop | CP2K capability (GPW / periodic) |
| Branch | `wf/ot-uks-kpoints` |
| Domain | `gpw` |
| Deps | `ot-uks-o2` ✓, `ot-kpoints` ✓, `gpw-kpoints-uks` ✓ |
| Public API | `cpu::gpw::GpwDft::{use_ot, n_alpha, n_beta, k_mesh / k_points}` |

## Intent

Productize **multi-k collinear UKS Orbital Transformation (OT)** on the complex
Bloch GPW path. Today:

- Γ-UKS OT is product (`reg_ot-uks-o2`);
- multi-k **RKS** OT is product (`reg_ot-kpoints`);
- multi-k **UKS** + `use_ot` **throws** at `GpwDft::build()`.

Close the residual: joint per-(k, spin) OT manifold (occ=1 open-shell
gradient), same self-consistent fixed point as multi-k UKS DIIS on identical
grids.

## Science

- OT parameterization: VandeVondele & Hutter, *J. Chem. Phys.* **118**, 4365
  (2003); Kühne et al., *J. Chem. Phys.* **152**, 194103 (2020) §VII C.
- Per spin σ ∈ {α, β}, occupied Bloch \(C_{k\sigma}\) on
  \(C^\dagger S C = I\); density \(P_k^\sigma = C_{k\sigma} C_{k\sigma}^\dagger\)
  (factor 1).
- Projected gradient \(g = (I - S C C^\dagger)(2\cdot\mathrm{occ}\, H C)\) with
  `occ = 1` (Hermitian multi-k kernels in `internal/ot.h`).
- Outer SCF: rebuild dual \(H_k^\alpha, H_k^\beta\) from k-summed spin densities;
  joint Polak–Ribière+ CG + Armijo on total energy over all (k, σ).
- k sampling: explicit list or Monkhorst–Pack; optional time-reversal fold.
- Units: method channel Hartree; Frame energy still eV when routed through an
  evaluator.

## API (no public Config / no factories)

```cpp
cpu::gpw::GpwDft m;
m.n_alpha = 2; m.n_beta = 0;   // or general UKS fillings
m.use_ot = true;
m.k_points = {{0.25, 0.0, 0.0}};
m.k_point_weights = {1.0};
// or m.k_mesh = {2,1,1};
m.build(atoms);
auto& r = m.run(atoms);  // r.e_total [Ha]
```

## Fail-fast (unchanged / residual)

| Combo | Behaviour |
| ----- | --------- |
| multi-k UKS OT | **product** (this slug) |
| multi-k OT + DFT+U | still throws (`ot-dftu-kpoints` residual) |
| OT + metallic smearing | still throws (`ot-smearing` residual) |
| multigrid / nonperiodic / `local_pp_on_grid=false` on k-path | still throws |

## Acceptance

1. `GpwDft::build` **accepts** `use_ot` + multi-k + UKS (integer occupations).
2. Unit: multi-k UKS OT energy matches multi-k UKS DIIS within **1e-5 Ha** on
   compact high-spin H₂ (κ = ¼ or MP 2×1×1).
3. Regression `regressions/gpw/ot-uks-kpoints`:
   - A1: explicit-Γ complex UKS OT ≡ real Γ UKS OT (1e-6 Ha)
   - A2: multi-k UKS OT ≡ multi-k UKS DIIS (κ=¼) (1e-4 Ha)
   - A3: MP 2×1×1 TR OT ≡ DIIS (1e-4 Ha)
   - A4: absolute band (frozen first-green or offline CP2K when available)
4. Dual gate green; `perf_exempt` if no fair CP2K twin wall, else
   \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\).
5. Update prior fail-fast gates (`test_cpu_gpw_ot` MultiKUksOtThrows,
   `reg_ot-kpoints` A4) to product semantics.

## Out of scope

- Multi-k OT + PLUS_U (`ot-dftu-kpoints`)
- OT + smear (`ot-smearing`)
- Forces/stress under OT (use DIIS path or existing force kernels)
- MPI industrial OT scaling

## Files

| Path | Change |
| ---- | ------ |
| `src/cpu/gpw/gpw_solver.cpp` | multi-k UKS OT loop; drop build reject |
| `include/atomiverse/cpu/gpw/gpw_solver.h` | doc: multi-k UKS OT product |
| `tests/cpu/gpw/test_cpu_gpw_ot.cpp` | flip throw test → OT/DIIS parity |
| `regressions/gpw/ot-kpoints/*` | A4 no longer requires throw |
| `regressions/gpw/ot-uks-kpoints/*` | new e2e |
| `tools/cp2k_oracle/inputs/ot-uks-kpoints/` | offline twin input |

# Spec — `pw-variable-cell-si`

## Goal

Productize plane-wave **Nielsen–Martin stress** into **variable-cell
relaxation** of diamond Si (NC-UPF), mirroring `gpw-variable-cell-ice`
(Minimize + cell optimizer + stress criterion + evaluator stress write).

## Scope (IN)

- `cpu::minimizer::Minimize` + isotropic cell step driven by PW intensive
  stress on `kStress` [eV/Å³].
- Diamond Si 2-atom fcc primitive, ONCVPSP-LDA `Si_lda_nc.upf`, tutorial
  ecut/k matching `pw-lda-si` / `pw-analytic-stress-si`.
- Offline pin: relaxed conventional lattice \(a^\*\) and residual
  hydrostatic stress; optional CP2K CELL_OPT GTH wall note (not absolute
  energy comparable).
- Stress convention: **mechanical / CP2K** \(\sigma = -(1/\Omega)\partial E/\partial\varepsilon\)
  so `IsotropicCellRelax` / `AnisotropicCellRelax` expand when compressed
  (same sign as GPW). Free `pw_stress` / FD oracle share this convention
  (dual-gate of analytic vs FD unchanged; pin sign flips vs pre-product
  energy-derivative pin).

## Scope (OUT)

- Joint ion+cell relax (fixed fractional coordinates only).
- Shear / triclinic free relax (KEEP_ANGLES / isotropic lattice).
- Production ecut (≥ 40 Ry) or dense k-mesh (budget ≤ 60 s).
- Molecular post-HF, classical MD, CUDA.

## Public API shape

```cpp
cpu::pw::PW method;
method.upf_path = ...;
method.ecut_ry = 15.0;
method.k_mesh = 2;
// ...
cpu::minimizer::Minimize opt;
cpu::minimizer::IsotropicCellRelax cell(/*P_target*/ 0.0, /*alpha*/ 1.0);
cell.stress_is_virial = false; // PW intensive kStress
cell.bulk_mod = /* measured B at a* */;
opt.set_optimizer_ref(cell);
opt.add_criterion(cpu::StressConverged{0.0, stol});
opt.on(minimizer::Phase::force)
    .add(cpu::pw::PWEvaluator{std::move(method)}.with_stress(true));
opt.build();
opt.run(frame, budget);
```

No public `*Config` bags; no `make_*` factories.

## Implementation plan

1. Unify `pw_stress` / `pw_stress_fd` sign to CP2K/GPW mechanical convention;
   update `pw-analytic-stress-si` pin sign and docs.
2. `Box::scale_isotropic`: true full-matrix scale (required for fcc
   primitive cells whose diagonal is zero).
3. `IsotropicCellRelax::stress_is_virial` (default `true` classical) —
   when false, \(P=\mathrm{tr}(\sigma)/3\) without dividing by \(V\).
4. Unit tests: intensive isotropic step; full-matrix isotropic scale on
   off-diagonal cell.
5. Regression `regressions/pw/pw-variable-cell-si/` vs offline pin.
6. Offline CP2K GTH CELL_OPT input under `tools/cp2k_oracle/inputs/` for
   wall note only.

## Acceptance

See `acceptance.md`.

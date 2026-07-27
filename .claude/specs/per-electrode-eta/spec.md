# per-electrode-eta — per-electrode Gaussian width η_e on ConstPotBO

Status: in_progress
Loop: metalwalls-method-loop (wave 5)
Queue row: `per-electrode-eta` | Per-electrode Gaussian width η_e (v1 shares one η)
| electrochem | multi-wall different η | Q, A diag vs shared-η twin |
oracle=metalwalls (no binary installed → literature pin per
`tools/metalwalls_oracle/README.md` contract, perf_exempt).

## Problem

MetalWalls (map §1 multi-electrode V,η) lets every electrode carry its own
Gaussian charge width η_e. Atomiverse `ConstPotBO` already has mixed-width
internals — `pair_A_(r, η_i, η_j)` implements
η_ij = η_i η_j / √(η_i² + η_j²) — but the public surface exposes exactly one
shared `eta`: the matrix build, RHS, and electrode-force paths all pass the
single member. Multi-wall systems with dissimilar electrode materials
(e.g. LAMMPS-default gold η = 1.979 Å⁻¹ vs Merlet/Salanne graphite
η = 0.955 Å⁻¹) cannot be expressed.

## Design (ARCH LAW compliant)

- No new Component, no Driver change, no Config bag, no factory.
- `ElectrodeGroup` gains `double eta = 0.0` (0 ⇒ inherit the shared member).
- New overload `ConstPotBO::add_electrode(int type, double psi, double eta)`
  (`&` and `&&`, mirroring the existing 2-arg pair); throws
  `std::invalid_argument` unless η > 0. The 2-arg overload is unchanged
  (shared-η v1 path stays byte-identical).
- `collect_electrodes_` resolves one η per electrode atom
  (`eta_vec_`, private workspace).
- Physics (Reed–Lanning–Madden JCP 126, 084704 (2007); Scalfi et al.
  JCP 153, 174704 (2020); MetalWalls JOSS 5(53), 2373 (2020)):
  - A_ii self term: k √(2/π) η_i (unchanged formula, per-atom η_i).
  - A_ij cross term: k erf(η_ij r)/r with η_ij = η_iη_j/√(η_i²+η_j²)
    (existing `pair_A_`, now fed per-atom widths) — including 2D xy images.
  - Same-site 2D images: k erf((η_i/√2) r)/r (Gaussian–Gaussian equal width).
  - RHS b_i: k q_j erf(η_i r)/r (existing `pair_b_`, per-atom η_i).
  - `write_electrode_forces`: force = −Q_i ∂/∂r [k q_j erf(η_i r)/r].
    **Iron-law fix folded in:** the shipped force kernel used α = η/√2
    while the RHS energy kernel uses erf(η r)/r — the force was not the
    gradient of the functional the solve used. Per-atom α = η_i now matches
    `pair_b_` exactly. (No numeric pin existed on the old α; the
    electrode-atom-forces regression only asserts non-zero forces.)
- Out of scope (documented residual, not silent): `ConstPotCG` keeps its
  shared-η v1 surface (independent duplicate of the matrix path — future
  queue row); `ConstPotMaZe` reaches per-electrode η through
  `solver().add_electrode(type, psi, eta)` without any change.

## Acceptance criteria

1. `add_electrode(type, psi, eta)` registers a per-group η; η ≤ 0 throws;
   2-arg overload behavior unchanged.
2. Shared-η twin equality: registering every group with the 3-arg overload at
   the shared value reproduces the 2-arg path charges exactly (≤ 1e-14 e).
3. A-diagonal observable: with per-group η, no images, no TF,
   `matrix_condition()` = η_max/η_min (analytic k√(2/π)η_i diagonal).
4. Vacuum two-site capacitor, η₁ = 1.979 Å⁻¹, η₂ = 0.955 Å⁻¹, d = 4 Å,
   ΔΨ = 1 V, electroneutrality: Q = ΔΨ/(A₁₁+A₂₂−2A₁₂) =
   0.037722100829787 e (literature formulas above, k = 14.399645 eV·Å/e²).
5. RHS per-η: single electrode site (group η ≠ shared member η) + one point
   charge, Ψ = 0, no neutrality ⇒ Q = −erf(η_g r)/r · q / (√(2/π) η_g);
   group η wins over the shared member.
6. Electrode force matches the analytic gradient of the RHS kernel with the
   group η (pins the α fix).
7. Unit tests (single-function, seconds) in
   `tests/cpu/dynamics/test_const_pot_per_eta.cpp`, registered in
   `tests/CMakeLists.txt`.
8. Regression `regressions/md/per-electrode-eta/` — out-of-library consumer,
   REFERENCE.md (oracle: literature + full citations, oracle_version,
   units), hard-coded pins, ≤ 60 s; targets `reg_per-electrode-eta` and
   `run_reg_per-electrode-eta` exist.
9. Oracle notes: `tools/metalwalls_oracle/pin_notes_per-electrode-eta.md` +
   `inputs/metalwalls/per-electrode-eta/` documenting the would-be MW twin.
10. MD-only dual gate green (HF/PW/GPW/CUDA/TORCH all OFF);
    perf_exempt=true (no runnable oracle binary for this physics).

---
title: Numerical Hessian — Response driver category + finite-difference Hessian components
status: approved
created: 2026-05-21
---

# Numerical Hessian — Response driver category + finite-difference Hessian components

## Summary

This spec adds the foundation of vibrational analysis to the Atomiverse CPU backend: a fourth outer driver category, `Response`, and the components that assemble a raw numerical (central-difference) molecular Hessian. After this link a user can register any `Capability::Force` evaluator (an FF evaluator or `cpu::scf::SCFEvaluator`) into `cpu::vib::Hessian`, call `run(frame)`, and receive a `cpu::vib::HessianResult` carrying the symmetrized 3N×3N second-derivative matrix of the energy. The new category is the genuine fourth outer-loop semantic that CLAUDE.md's driver-taxonomy section explicitly sanctions ("Response for finite-difference Hessian"): it runs a fixed `6N` displacement budget with no `outer_status` polling — every displacement must complete for the Hessian to be defined. Frequencies, mass-weighting and normal modes (link 2) and thermochemistry (link 3) build on `HessianResult` and are out of scope here.

## Domain basis

The molecular Hessian is the matrix of second derivatives of the Born–Oppenheimer energy with respect to nuclear Cartesian coordinates:

```
H_{iα,jβ} = ∂²E / ∂R_{iα} ∂R_{jβ}
```

a symmetric 3N×3N matrix, where `i,j` index atoms and `α,β ∈ {x,y,z}`. In Atomiverse internal physics units the natural Hessian unit is Hartree/Bohr²; however `cpu::Frame` stores positions in Ångström and forces in eV/Å, and the registered evaluator writes forces (the negative gradient `g = −∂E/∂R`) into `frame.atoms.fx/fy/fz` in eV/Å. To keep the implementation consistent with the existing force pipeline, this link assembles the Hessian directly in the **frame's working units, eV/Å²**, and records the unit in `HessianResult`. Unit conversion to Hartree/Bohr² is deferred to link 2 (mass-weighting / frequency conversion), where the conversion constants are needed anyway.

The Hessian is built by central-differencing the gradient that the registered evaluator already produces. For each Cartesian coordinate `iα` we displace the atom by ±h, evaluate the force, and form one Hessian row:

```
H[iα, :] = −( f(R + h·e_{iα}) − f(R − h·e_{iα}) ) / (2h)
```

Here `f` is the frame force (`−g`), so `−Δf/(2h) = Δg/(2h)` is the second derivative. A full Hessian costs `6N` force evaluations (two per coordinate). Central differencing has truncation error `O(h²)` with leading coefficient `g'''/6`.

Step size. The draft proposes `h = 1e-3 Å`. Displacement is applied **directly in the frame's Ångström coordinates**: `h = 1.0e-3 Å` exactly. For callers reasoning in Bohr, the equivalent is `h = 1.889726e-3 Bohr` (1 Bohr = 0.52917721 Å) — recorded as a non-binding note; the binding value is `1.0e-3 Å`.

Symmetrization. The raw row-wise Hessian is only approximately symmetric (independent gradient noise per row). The standard treatment, used by Gaussian, PySCF and ASE, is the orthogonal projection onto the symmetric subspace:

```
H ← (H + Hᵀ) / 2
```

This is the least-squares-optimal symmetric approximation of a noisy near-symmetric matrix.

Accuracy requirement. The dominant error source is gradient noise `δg` amplified by the `1/(2h)` factor. Therefore the energy/gradient evaluation at **every displaced geometry** must be converged at least as tightly as at the equilibrium geometry — for SCF evaluators, gradient tolerance ≤ 1e-7 Ha/Bohr. The acceptance contract enforces this by requiring the driver-produced Hessian to match an in-test oracle (same method, same `h`, same convergence) to ≤ 1e-6 in Hessian units.

References:
- Wilson, Decius & Cross, *Molecular Vibrations*, McGraw-Hill (1955) — Cartesian Hessian / GF method foundation.
- Pulay, *Mol. Phys.* 17, 197 (1969), DOI:10.1080/00268976900100941 — analytic forces differenced here.
- Helgaker, Jørgensen & Olsen, *Molecular Electronic-Structure Theory*, Wiley (2000), Ch. 9–12 — molecular gradients and Hessians.
- Neugebauer & Reiher, *J. Comput. Chem.* 25, 587 (2004), DOI:10.1002/jcc.10374 — finite-difference Hessian construction in practice.

## Design

A new outer category `Response` is added as a header-only base `cpu::ResponseDriver` at the cpu backend root, alongside `EvaluateDriver` / `DynamicsDriver` / `OptimizeDriver`. It subclasses `CPUDriver<ATV::vib::Phase>` (which wraps `Driver<vib::Phase, cpu::device_tag>`). Its `run(cpu::Frame&)` derives `N = frame.atoms.size()` internally, runs `setup_pass(frame)`, then `step_once(d)` for `d` in `0 … 6N-1` (the displacement passes), then one final `step_once(6N)` pass (the diagonalize/finalize pass where `Symmetrize` runs), then `teardown_pass()`. `6N` is derived in `run()` — never a `run()` argument and never a Config field. `ResponseDriver` performs no `outer_status` polling: a Hessian must complete all displacements, so `iteration_status.h` plays no role in this category. `ResponseDriver` overrides only `validate_components()`.

The phase enum `ATV::vib::Phase` is added to `facade/stage.h` next to the existing `md` / `minimizer` / `evaluate` / `scf` phase namespaces, with `phase_names()`:

```
enum class Phase { setup, displace, force, accumulate,
                   diagonalize, thermochem, end_of_step, COUNT };
```

`force` is the slot where the user-registered evaluator runs; `thermochem` is reserved for link 3 and carries no component in this link. A `CategoryTraits<ATV::vib::Phase>` specialization is added in `cpu/cpu_driver.h` (`eval_phase = vib::Phase::force`, `end_phase = vib::Phase::end_of_step`) — without it `CPUDriver<vib::Phase>` does not compile, and the inherited categorical setters (`set_pair`, `set_thermo`, …) need it.

A new capability bit `Capability::Hessian = 1u << 7` is added to `facade/component_role.h` (7 bits currently used). `HessianAccumulator` advertises `capabilities() == Hessian`; `ResponseDriver::validate_components()` requires exactly one such component **and** at least one component carrying `Capability::Force`. No new `ComponentRole` is introduced — all new structural components return `role() == None`, matching `Zero*` / `KineticEnergy`.

The concrete user-facing driver is `cpu::vib::Hessian` in the new domain subdir `include/atomiverse/cpu/vib/`, mirroring the `cpu::md::MD` / `cpu::minimizer::Minimize` precedent (domain-scoped, not root-level). It auto-registers the four structural components in its constructor and exposes the result.

Frame purity is preserved: no `std::optional<FrameHessian>` is added to `cpu::Frame`. Instead a plain `cpu::vib::HessianResult` struct (`include/atomiverse/cpu/vib/hessian_result.h`) holds the dense Hessian matrix (`xt::xtensor<double,2>`, shape `[3N,3N]`), an `n_atoms` field, and a `units` tag string ("eV/Angstrom^2"); link 2/3 extend it with frequencies/modes/thermo. `cpu::vib::Hessian::run(frame)` returns a `HessianResult` by value (mirroring `OptimizeDriver::run()` returning `int`) and also caches it for `last_result()`. The four components write into a shared `[3N,3N]` workspace owned by the driver; the workspace is reachable from `LaunchContext` (driver stores it and the components access it through a driver-held handle, exactly as the driver owns its component set today).

New structural components (all `role() == None`), each reading the displacement index from `LaunchContext::step`:
- `ReferenceGeometry` — bound at `vib::Phase::setup`; on `step == 0` caches the equilibrium coordinates `R₀`; idempotent thereafter.
- `DisplaceStep` — bound at `vib::Phase::displace`; for displacement index `s = ctx.step` (`0 ≤ s < 6N`) decodes `(atom i, axis α, sign ±)` and sets `frame` coordinates to `R₀ ± h·e_{iα}`; restores `R₀` on the final pass. Uses `h = 1.0e-3 Å`.
- `HessianAccumulator` — bound at `vib::Phase::accumulate`; `capabilities() == Hessian`; after the evaluator has written forces for displacement `s`, stores the force vector and, once both `±h` evaluations for coordinate `iα` are available, writes Hessian row `H[iα,:] = −(f₊ − f₋)/(2h)` into the workspace.
- `Symmetrize` — bound at `vib::Phase::diagonalize` with order `-1` (runs on the final `6N` pass before any future diagonalizer); replaces the workspace with `(H + Hᵀ)/2`.

Layering: `cpu/vib/hessian.h` includes no `cpu/scf/` header. The force evaluator (`SCFEvaluator` or an FF evaluator) is registered by user code via `on(vib::Phase::force).add(...)`, exactly as AIMD registers `SCFEvaluator` into `cpu::md::MD`. No domain-to-domain header edge is created.

The new `src/cpu/vib/` directory is registered in the top-level `CMakeLists.txt` next to the existing explicit `src/cpu/md/*.cpp` listing.

## Files to create or modify

- `include/atomiverse/facade/stage.h` — add `ATV::vib::Phase` enum + `phase_names()`.
- `include/atomiverse/facade/component_role.h` — add `Capability::Hessian = 1u << 7`.
- `include/atomiverse/cpu/cpu_driver.h` — add `CategoryTraits<ATV::vib::Phase>` specialization.
- `include/atomiverse/cpu/response_driver.h` (new) — `cpu::ResponseDriver` category base.
- `include/atomiverse/cpu/vib/hessian_result.h` (new) — `cpu::vib::HessianResult` struct.
- `include/atomiverse/cpu/vib/hessian.h` (new) — `cpu::vib::Hessian` concrete driver + shared workspace handle.
- `include/atomiverse/cpu/vib/reference_geometry.h` (new) — `ReferenceGeometry` component.
- `include/atomiverse/cpu/vib/displace_step.h` (new) — `DisplaceStep` component.
- `include/atomiverse/cpu/vib/hessian_accumulator.h` (new) — `HessianAccumulator` component.
- `include/atomiverse/cpu/vib/symmetrize.h` (new) — `Symmetrize` component.
- `src/cpu/vib/hessian.cpp` (new) — `cpu::vib::Hessian::run` and component definitions if non-trivial.
- `CMakeLists.txt` — register `src/cpu/vib/hessian.cpp` in `ATOMIVERSE_SRCS`.
- `tests/cpu/vib/test_numerical_hessian.cpp` (new) — failing tests for the components and driver.
- `tests/CMakeLists.txt` — register the new test target `test_cpu_vib_numerical_hessian`.

## Tasks

- [ ] Write failing tests for the Response category + numerical Hessian (`tests/cpu/vib/test_numerical_hessian.cpp`), register `test_cpu_vib_numerical_hessian` in `tests/CMakeLists.txt`
- [ ] Add `ATV::vib::Phase` enum + `phase_names()` in `include/atomiverse/facade/stage.h` and `Capability::Hessian` in `include/atomiverse/facade/component_role.h`
- [ ] Add `CategoryTraits<ATV::vib::Phase>` specialization in `include/atomiverse/cpu/cpu_driver.h`
- [ ] Implement `cpu::ResponseDriver` (6N-displacement loop, no outer_status polling, `validate_components()`) in `include/atomiverse/cpu/response_driver.h`
- [ ] Implement `cpu::vib::HessianResult` in `include/atomiverse/cpu/vib/hessian_result.h` and the structural components `ReferenceGeometry` / `DisplaceStep` / `HessianAccumulator` / `Symmetrize` in their respective `include/atomiverse/cpu/vib/*.h` headers
- [ ] Implement `cpu::vib::Hessian` driver (`include/atomiverse/cpu/vib/hessian.h` + `src/cpu/vib/hessian.cpp`) wiring the four components and returning/caching `HessianResult`; register `src/cpu/vib/hessian.cpp` in `CMakeLists.txt`
- [ ] Add doxygen `@brief`/`@param`/`@return` with units to every new public symbol per the project doc style
- [ ] Verify the driver Hessian against an in-test central-difference-of-gradient oracle for H2O (same method/h/convergence) to ≤ 1e-6 in Hessian units
- [ ] Run full check + test suite

## Testing strategy

- Happy path: build a `cpu::vib::Hessian`, register a `Capability::Force` evaluator on `vib::Phase::force`, run on a small molecule, and confirm `HessianResult.hessian` has shape `[3N,3N]` and finite entries.
- Self-consistency (domain validation): for H2O at RHF or PBE with a small basis (STO-3G or 6-31G), the driver-produced Hessian matches an oracle that central-differences `scf_gradient` directly in the test (same `h = 1.0e-3 Å`, same method, same convergence) to ≤ 1e-6 eV/Å².
- Symmetry (domain validation): the returned Hessian satisfies `‖H − Hᵀ‖∞ ≤ 1e-12` (Symmetrize correctness).
- Structural / rank (domain validation): at a near-stationary geometry, the symmetrized Hessian's translational-block row sums are near zero (`‖Σ_j H[iα,jα]‖ ≤ 1e-4` per axis), i.e. it satisfies the translational-invariance sum rule — a cheap structural check; the full 6-near-zero-eigenvalue analysis is deferred to link 2.
- Validation: `ResponseDriver::validate_components()` throws when no `Capability::Force` evaluator is registered, and throws when no `Capability::Hessian` component is present.
- Edge case: a single-atom frame (`N = 1`) runs `6` displacement passes without crashing and yields a `[3,3]` Hessian.
- Regression: the full `ctest` suite stays green — the new category and `Capability::Hessian` bit must not regress any existing target.

## Out of scope

- Frequencies, mass-weighting, mass-weighted Hessian and eigen-frequency conversion — link 2 (`phase4-vibrations-02-frequencies`).
- Thermochemistry (ZPE, partition functions, entropy) — link 3 (`phase4-vibrations-03-thermochemistry`).
- Analytic (second-derivative) Hessian — not planned in this chain; the finite-difference path is the deliberate first implementation.
- CUDA `ResponseDriver` — CUDA driver taxonomy migration is tracked separately in `cuda-driver-taxonomy.md`.
- Normal-mode trajectory output / `NormalModeWriter` Recorder — link 2/3.
- Wall-time / performance budgets — the draft's wall-time targets are recorded as non-binding notes; this link's contract is correctness-only.

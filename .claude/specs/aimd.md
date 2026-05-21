---
title: First-class AIMD — set_qm composition + faithful ASPC density warm-start
status: approved
created: 2026-05-21
---

# First-class AIMD — set_qm composition + faithful ASPC density warm-start

## Summary

Ab initio molecular dynamics (Born-Oppenheimer MD) already runs in Atomiverse by adding a `cpu::scf::SCFEvaluator` to a `cpu::md::MD` driver, but the composition is hand-wired (`sim.on(md::Phase::force).add(SCFEvaluator{...})`) and the cross-step density warm-start it depends on is built on a misnamed, mathematically-incorrect predictor. This spec makes AIMD a first-class, ergonomic capability without introducing any new driver or outer category: it adds a `set_qm()` categorical setter to `cpu::md::MD` (consistent with `set_pair` / `set_integrator`, and mutually exclusive with `set_pair` because a step's force comes from a classical pair potential OR the QM evaluator, never both), and fixes `ASPCPredictor` to implement Kolafa's actual Always-Stable Predictor-Corrector — replacing the current Lagrange extrapolation in `predict()` with the true ASPC `B_j` binomial coefficients, and correcting a docstring that misattributes ASPC to second-generation Car-Parrinello. The result is a documented, validated AIMD workflow whose energy-conservation properties rest on a time-reversible predictor.

## Domain basis

**Born-Oppenheimer MD.** Classical nuclei move on the electronic ground-state surface. Each MD step solves the electronic SCF problem to self-consistency; the potential energy is the converged SCF energy `E_pot = E_SCF`, and the nuclear force is the analytic SCF gradient `F = -dE_pot/dR`. AIMD is therefore classical MD with a quantum force provider — no new dynamics. Velocity-Verlet (the existing integrator) is symplectic and time-reversible; in NVE the conserved quantity is `E_total = E_kinetic + E_pot(SCF)`. Ref: Lin, Lu, Shao, *Entropy* **16**(1), 110 (2014), DOI 10.3390/e16010110.

**Energy conservation is the defining BOMD validation.** A well-behaved NVE BOMD code shows total-energy drift of order `1e-5 eV/ps/atom` or tighter (literature PBE NVE runs reported 1.4e-5 and 4.3e-5 eV/ps/atom, called "acceptable"). A poorly-converged or non-time-reversible scheme drifts orders of magnitude more (`>= 1e-3 eV/ps/atom`, visible as monotonic heating within a few ps). Heating-rate identity: `dT/dt = (dE_total/dt) / (0.5 * N_dof * k_B)`, with `k_B = 8.617333e-5 eV/K`.

**SCF tolerance couples to conservation.** An incompletely converged SCF gives forces that are not the exact gradient of the computed energy (Hellmann-Feynman invalid; a nonzero residual-force term remains), so the dynamics becomes non-conservative and `E_total` drifts; the drift grows monotonically as the SCF tolerance is loosened. The standard requirement is SCF converged to `~1e-6 a.u.` or tighter, tied to the DIIS/commutator (gradient) error, not just the energy change `dE`. AIMD configurations must therefore use an SCF tolerance no looser than `1e-6`.

**Timestep.** For H-containing systems `dt ~ 0.25-0.5 fs` is required to resolve the ~3000-4000 cm⁻¹ X-H stretch (a 4000 cm⁻¹ period is ~8.3 fs, i.e. 10-20 steps/period).

**ASPC density warm-start.** Kolafa, *J. Comput. Chem.* **25**, 335 (2004), DOI 10.1002/jcc.10385. The Always-Stable Predictor-Corrector supplies a cheap, accurate SCF initial guess from trajectory history; it is purely an initial-guess accelerator (the SCF is still converged every step) and is **not** a dynamics scheme — it is distinct from Car-Parrinello / extended-Lagrangian BOMD, which the current docstring wrongly conflates with it via the Kühne PRL 98 066401 citation.

Predictor (history newest-first `D_n, D_{n-1}, …`), with `nvec = K + 1` history terms:

```
D_pred(n+1) = sum_{j=1..nvec} B_j * D_{n+1-j}
B_1 = (4*nvec - 2) / (nvec + 1)
B_i = (-1)^(i+1) * i * C(2*nvec, nvec-i) / C(2*nvec-2, nvec-1)   for i = 2..nvec
```

Coefficient sets (each sums to 1):

| nvec | B_j |
|------|-----|
| 2 | {2, -1} |
| 3 | {2.5, -2, 0.5} |
| 4 | {2.8, -2.8, 1.2, -0.2} |
| 5 | {3, -3.428571, 1.928571, -0.571429, 0.071429} |
| 6 | {3.142857, -3.928571, 2.619048, -1.047619, 0.238095, -0.023810} |

Corrector damping (already correct in the code): `omega = (K+2)/(2K+3)`; `D_store = omega*D_SCF + (1-omega)*D_pred`. ASPC is time-reversible to predictor order (long-time energy-conservation error of order `2K-1`); a non-time-reversible extrapolation (such as plain Lagrange for `K >= 2` — Lagrange `K=2` gives `{3,-3,1}`, while true ASPC `nvec=3` gives `{2.5,-2,0.5}`) injects a hysteresis/memory force bias that drifts energy even at fixed SCF tolerance. Density-matrix idempotency: a linear combination of idempotent `P` is not itself idempotent, so `D_pred` is used **only** as the SCF initial guess and the SCF cycle restores idempotency — Atomiverse keeps this approach (it already does so; `SCFEvaluator` feeds `D_pred` through `SCFInitialGuess` into `SCFSolver`).

**Unit handoff.** `1 Hartree = 27.211386245988 eV`; `1 Hartree/Bohr = 51.42208619 eV/Å`; constants in `core/units.h`. Force is the negative energy gradient. `SCFEvaluator` already performs this handoff correctly (it works inside `cpu::Evaluate` and `cpu::minimizer::Minimize`); this spec confirms it but does not change it.

## Design

**No new types in the driver/component taxonomy.** AIMD remains exactly `cpu::md::MD` (a `DynamicsDriver`) composed with `cpu::scf::SCFEvaluator` (an `Evaluator` component with `Capability::Force`). The work is one categorical setter, one corrected algorithm, and tests.

**`MD::set_qm`.** A new MD-specific categorical setter, declared alongside `set_integrator` in `cpu/md/md.h`. It accepts an `SCFEvaluator` (rvalue or lvalue), pins it at `ATV::md::Phase::force @ 0` (the same slot `set_pair` uses for the classical evaluator), and registers it through the existing singular-slot machinery so a second `set_qm` with a different component throws. `set_qm` and `set_pair` are **mutually exclusive**: MD owns a single "force-source" identity. The slot is implemented so that whichever of `set_pair` / `set_qm` is called first claims the force source; a subsequent call to the other throws `std::invalid_argument` naming both setters. `set_bond` / `set_angle` / `set_dihedral` are *not* blocked — bonded terms may legitimately coexist, and the mutual exclusion is specifically the pair-vs-QM force source. Because `MD::set_pair` is inherited from `CPUDriver` and uses the private `pair_set_` slot, the cross-check is implemented by adding a shared force-source guard: a protected `claim_force_source_(const char* setter)` helper on `CPUDriver` that both `set_pair` and `MD::set_qm` call before registering, throwing if a different setter already claimed it.

**`SCFEvaluator` integration.** `set_qm` does not change `SCFEvaluator`. The evaluator already: rebuilds integrals each step, runs the SCF through a pre-built `SCFSolver`, computes the analytic gradient, writes `frame.atoms.f*` and `frame.energy`, and (when `predictor()` was called) warm-starts the next step's density through `ASPCPredictor`. `set_qm` is a thin ergonomic wrapper over the existing `on(Phase::force).add(...)` path.

**`ASPCPredictor` correction.** `compute_coefficients()` currently fills the coefficient array with Lagrange extrapolation coefficients. Replace this with the Kolafa `B_j` formula above (the docstring must not say "Lagrange"). `predict()` keeps its current structure (`D_pred = sum_k coeff_k * history_k`) but consumes the corrected coefficients; the history-ordering convention (`history_` newest term, `B_1` multiplies the newest) must be made explicit and asserted. `correct()` is already the correct Kolafa damping and is unchanged. The class header docstring is reframed: ASPC is a Born-Oppenheimer SCF *initial-guess accelerator*, not a second-generation Car-Parrinello dynamics scheme; the Kühne PRL citation is removed; the `order` parameter doc no longer says "Lagrange". `ASPCPredictor::order()`'s doc comment ("Lagrange extrapolation order") is corrected to "ASPC extrapolation order K".

**Existing ASPC wiring reconciliation.** `SCFEvaluator::launch()` already calls `predict()` / `correct()` / `record()` per step, feeding the predicted density into `SCFInitialGuess`, so ASPC *is* a live cross-MD-step density guess today — not dead, not intra-SCF-only. The bug is purely the coefficient set inside `predict()`; the cross-step plumbing is correct and stays.

## Files to create or modify

- `include/atomiverse/cpu/md/md.h` — add `MD::set_qm` categorical setter + class docstring update.
- `include/atomiverse/cpu/cpu_driver.h` — add protected `claim_force_source_` guard; `set_pair` calls it.
- `include/atomiverse/cpu/scf/aspc_predictor.h` — correct docstrings (remove Car-Parrinello / Kühne citation, drop "Lagrange"), update `predict()` / `order()` / ctor docs to ASPC.
- `src/cpu/scf/aspc_predictor.cpp` — replace `compute_coefficients()` Lagrange body with Kolafa `B_j` formula; update file-header comment.
- `tests/cpu/scf/test_aspc_predictor.cpp` (new) — unit tests for the `B_j` coefficient table and predictor correctness.
- `tests/cpu/md/test_aimd_set_qm.cpp` (new) — `set_qm` ergonomics + `set_pair`/`set_qm` mutual-exclusion tests, and an NVE AIMD run via `set_qm`.
- `tests/cpu/scf/test_aimd_conservation.cpp` — extend with an ASPC warm-start SCF-iteration-count test (warm-start vs cold-start mean iterations/step).
- `tests/CMakeLists.txt` — register the two new test files.

## Tasks

- [ ] Write failing tests for ASPC Kolafa B_j coefficients (tests/cpu/scf/test_aspc_predictor.cpp) — assert computed coefficients match the verified nvec=2..6 table to 1e-12, each set sums to 1, and a linear-ramp density history is reproduced exactly by predict()
- [ ] Replace Lagrange extrapolation with Kolafa B_j coefficients in ASPCPredictor::compute_coefficients (src/cpu/scf/aspc_predictor.cpp); make history-ordering convention explicit
- [ ] Correct ASPCPredictor docstrings in include/atomiverse/cpu/scf/aspc_predictor.h per doc.style — reframe as a Born-Oppenheimer SCF initial-guess accelerator with units, remove the second-gen Car-Parrinello / Kühne PRL citation, drop "Lagrange"
- [ ] Write failing tests for MD::set_qm and set_pair/set_qm mutual exclusion (tests/cpu/md/test_aimd_set_qm.cpp) — set_qm registers an SCFEvaluator at the force phase, double set_qm with a different component throws, set_pair-then-set_qm and set_qm-then-set_pair both throw naming both setters
- [ ] Implement claim_force_source_ guard in CPUDriver and MD::set_qm in include/atomiverse/cpu/md/md.h + include/atomiverse/cpu/cpu_driver.h; update the MD class docstring per doc.style to document the QM force source
- [ ] Add an NVE AIMD energy-conservation assertion driven through set_qm in tests/cpu/md/test_aimd_set_qm.cpp — H2 STO-3G, dt=0.5 fs, assert E_total drift below the stated bound
- [ ] Extend tests/cpu/scf/test_aimd_conservation.cpp with an ASPC warm-start test — assert mean SCF iterations/step over an AIMD trajectory with predictor() enabled is strictly lower than cold-start (no predictor)
- [ ] Register tests/cpu/scf/test_aspc_predictor.cpp and tests/cpu/md/test_aimd_set_qm.cpp in tests/CMakeLists.txt
- [ ] Run full check + test suite

## Testing strategy

- **Happy path — `set_qm`.** `MD::set_qm(SCFEvaluator{config})` followed by `set_integrator` + `build()` + `run()` produces a finite, bonded-well energy on an H2 STO-3G fixture, equivalent to the existing hand-wired `on(Phase::force).add(...)` composition.
- **Edge case — mutual exclusion.** `set_pair` then `set_qm` throws `std::invalid_argument`; `set_qm` then `set_pair` throws; the message names both setters. A second `set_qm` with a *different* component throws; re-calling `set_qm` with the same component name is idempotent (consistent with the existing singular-slot semantics).
- **Edge case — ASPC history.** `predict()` returns false until `K+1` densities are recorded; after a `reset()` it again returns false until refilled.
- **Domain validation — ASPC coefficients.** Computed `B_j` for `nvec = 2..6` match the verified table to 1e-12 and each set sums to 1. A density history that is an exact linear function of step index is reproduced exactly by `predict()` (ASPC is exact for linear input — this distinguishes a correct implementation from a broken one without depending on the Lagrange-vs-ASPC difference, which a separate quadratic-input check exercises: ASPC `nvec=3` `{2.5,-2,0.5}` applied to a quadratic history gives the ASPC value, not the Lagrange `{3,-3,1}` value).
- **Domain validation — NVE energy conservation.** An H2 STO-3G NVE AIMD run driven through `set_qm` (dt = 0.5 fs, 200 steps = 100 fs) keeps the `E_total` peak-to-peak envelope within the bound used by the existing `test_aimd_conservation.cpp` (≤ 0.10 eV), consistent with literature drift ~1e-5..1e-4 eV/ps/atom scaled to a 0.1 ps trajectory. Energy must stay in the bonded H2 well.
- **Domain validation — warm-start efficiency.** Over an AIMD trajectory (≥ 30 steps past predictor fill), the mean SCF iterations/step with `predictor()` enabled is strictly lower than with cold-start. Forces from the warm-started path still agree with a fresh single-point at a mid-trajectory geometry within 1e-5 eV/Å (already covered by an existing test; the corrected `B_j` must not regress it).
- **Regression.** Full `ctest` suite stays green, including the existing AIMD-conservation and ASPC-vs-fresh-single-point tests.

## Out of scope

- Car-Parrinello / extended-Lagrangian MD (a different method that adds fictitious electronic dynamics and dissipation) — AIMD here is strictly Born-Oppenheimer.
- Path-integral MD.
- NPT-AIMD — only NVE and NVT are in scope; NVT composes unchanged via the existing `NoseHoover` / `FixLangevin` thermostats with no code change, so no NPT-specific work is added here.
- CUDA AIMD — `SCFEvaluator` is CPU-only; the CUDA SCF path is tracked separately under `cuda-driver-taxonomy`.
- Multiple-time-step / RESPA integration for AIMD.
- A categorical `set_qm` on `cpu::Evaluate` or `cpu::minimizer::Minimize` — those already accept `SCFEvaluator` via `add()`; only `MD` gets the setter because only `MD` has the pair-vs-QM force-source conflict to manage.

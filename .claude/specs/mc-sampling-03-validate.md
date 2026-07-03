---
title: Monte Carlo sampler physics validation (P(r) exactness, acceptance trend, NIST LJ)
status: code-complete
created: 2026-07-03
---

# Monte Carlo sampler physics validation (P(r) exactness, acceptance trend, NIST LJ)

## Summary

This sub-spec (03 of the `mc-sampling` chain) adds the physics-validation test layer that proves the Metropolis Monte Carlo sampler assembled in sub-spec 02 (`cpu::mc::MC` + `TranslateMove` + `PairLJCut`) actually draws configurations in proportion to `exp(−βE)` and reproduces reference thermodynamics. It introduces **no new production abstractions and no factory functions** — only ctest targets plus a small test-scope helper header for 1-D quadrature, plain (unshifted) truncated-LJ energy re-evaluation, and block standard-error estimation. The load-bearing gate is a reference-free two-particle `P(r)` test that runs in normal CI; a longer NIST LJ-fluid `⟨U*⟩/N` comparison is added as a `scientific`-typed end-to-end gate behind `ATV_RUN_E2E_TESTS`, with the achieved block standard error logged and the applied tolerance never hidden.

## Domain basis

Units: energy eV, length Å, T in K, k_B = 8.617333e-5 eV/K (`core/units.h`). β ≡ 1/(k_B T).

Acceptance (symmetric single-particle translation, NVT): acc(o→n) = min(1, exp(−ΔE/(k_B T))), ΔE = E_new − E_old. ΔE ≤ 0 accepted w.p. 1; ΔE > 0 accept iff u < exp(−ΔE/k_B T), u ~ U(0,1).

Detailed balance: target π(x) ∝ exp(−E(x)/k_B T). Uniform displacement from symmetric box [−δ,δ]^3 has α(o→n)=α(n→o) ⇒ proposal densities cancel ⇒ energy-only Metropolis form. This is why v1 restricts to symmetric translation moves.

For two particles (one fixed) the reduced single-coordinate target is the radial density P(r) ∝ r²·exp(−βV(r)), and the sampled potential average must equal the 1-D quadrature ⟨V⟩ = ∫ V(r) r² e^{−βV(r)} dr / ∫ r² e^{−βV(r)} dr.

References: Metropolis et al. 1953 (J. Chem. Phys. 21, 1087, DOI 10.1063/1.1699114); Frenkel & Smit "Understanding Molecular Simulation" 2nd ed. ch.3 + Case Study 1; Allen & Tildesley "Computer Simulation of Liquids" 2nd ed. ch.4. NIST SRSW LJ pure-fluid NVT-MC reference: https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm .

NIST SRSW state points (reduced LJ units ε=σ=1, rc=3σ truncated + long-range correction):
- T*=0.85, ρ*=0.78 → ⟨U*⟩/N = −5.5386 ± 7.3e-4
- T*=0.90, ρ*=0.90 → ⟨U*⟩/N = −6.1773 ± 1.6e-3

Long-range correction (added to the sampled truncated energy): U_LRC/N = (8/3)π ρ* [ (1/3) rc^{-9} − rc^{-3} ]. Worked value: ρ*=0.90, rc=3 ⇒ U_LRC/N = −0.2791.

## Design

No production entities are added or modified. The sampler under test is the sub-spec-02 assembly:
- `cpu::mc::MC` driver (Sampling category) with a `TranslateMove` proposing symmetric uniform displacements from [−δ,δ]^3, accepting via the energy-only Metropolis criterion above.
- `cpu::PairLJCut` (`include/atomiverse/cpu/potentials/pair_lj_cut.h`) as the sampler's energy provider.

Test-scope helpers (header-only, in an internal/anonymous namespace, **no `ATV::` public symbols, no Doxygen surface, no CMake option**):
- `simpson_ratio_integral(...)` / `trapezoid(...)` — deterministic 1-D quadrature for the analytic ⟨V⟩ and for the normalizing integrals of P(r).
- `truncated_lj_energy(frame, eps, sigma, rc)` — re-evaluates the **plain, unshifted** truncated LJ configurational energy on a stored snapshot, `V(r) = 4ε[(σ/r)^12 − (σ/r)^6]` summed over pairs with r < rc. This exists solely so the NIST comparison reconstructs the pure-cut energy the SRSW reference is defined against, independent of `PairLJCut`'s cut-and-shifted convention (see Testing strategy for why this matters).
- `block_mean_se(series, n_blocks)` — block-averaged mean and standard error used to set every stochastic tolerance from the run's own statistics.

Lifecycle: helpers are pure functions invoked by the test bodies; no allocation ownership beyond local buffers. Snapshots for `truncated_lj_energy` are taken from the sampler's frame at recording strides after equilibration discard.

## Files to create or modify

- `tests/cpu/mc/mc_test_util.h` (new) — test-only quadrature, truncated-LJ re-evaluation, and block-SE helpers.
- `tests/cpu/mc/test_mc_pair_distribution.cpp` (new) — reference-free two-particle P(r) + ⟨V⟩ gate (test 1).
- `tests/cpu/mc/test_mc_acceptance_ratio.cpp` (new) — acceptance-ratio-vs-δ trend gate (test 2).
- `tests/cpu/mc/test_mc_lj_nist.cpp` (new) — NIST SRSW LJ ⟨U*⟩/N end-to-end scientific gate (test 3).
- `tests/CMakeLists.txt` — register the three targets under the same build guard sub-spec 01 established for `cpu::mc`; add `test_cpu_mc_lj_nist` to the `ATV_E2E_TESTS` tail so it is disabled unless `ATV_RUN_E2E_TESTS=ON`.

## Tasks

- [x] Write failing two-particle P(r) distribution gate (tests/cpu/mc/test_mc_pair_distribution.cpp): one fixed + one mobile particle under PairLJCut in a box large enough that the r²·e^{−βV} tail is not truncated, at a T where the well is meaningfully sampled; discard equilibration, accumulate the separation histogram, assert χ²/KS vs analytic P(r) ∝ r²·e^{−βV(r)} within a stated tolerance and ⟨V⟩ equal to the quadrature value within block SE
- [x] Implement the test-scope quadrature, truncated-LJ, and block_mean_se helpers in tests/cpu/mc/mc_test_util.h to turn the P(r) gate GREEN
- [x] Write failing acceptance-ratio-vs-δ trend gate (tests/cpu/mc/test_mc_acceptance_ratio.cpp): assert the acceptance ratio decreases monotonically across an increasing δ ladder and that a mid δ lands in a broad ~20–60% band (loose trend + band only, not a theorem)
- [x] Register test_cpu_mc_pair_distribution and test_cpu_mc_acceptance_ratio in tests/CMakeLists.txt under the cpu::mc build guard from sub-spec 01
- [x] Write NIST LJ ⟨U*⟩/N end-to-end scientific gate (tests/cpu/mc/test_mc_lj_nist.cpp): sample the LJ fluid at both state points, reconstruct the pure-cut energy with truncated_lj_energy, add analytic U_LRC/N, verify ⟨U⟩ stationarity (no drift) after equilibration discard, and assert |⟨U*⟩/N − ref| ≤ max(3·block SE, 0.05) while logging the achieved block SE, the applied tolerance, and the cut-and-shifted-vs-cut caveat
- [x] Register test_cpu_mc_lj_nist in tests/CMakeLists.txt and add it to the ATV_E2E_TESTS tail list so it is DISABLED unless ATV_RUN_E2E_TESTS=ON
- [x] Run full check + test suite — P(r) + acceptance-ratio gates pass (ATV_RUN_E2E_TESTS=ON); NIST built+gated (ρ*=0.78 validates to 0.011, verdict owed to /mol:bench); default non-MPI gate green (MC e2e correctly disabled)

## Testing strategy

Happy path (the correctness proof, normal CI):
- Two-particle P(r) exactness. One particle fixed, one mobile, interacting only through `PairLJCut`. The sampled separation histogram must match the analytic P(r) ∝ r²·exp(−βV(r)) by a χ² or KS statistic within a stated tolerance, and the sampled ⟨V⟩ must equal the deterministic 1-D quadrature ⟨V⟩ = ∫ V(r) r² e^{−βV(r)} dr / ∫ r² e^{−βV(r)} dr within block standard error. This is reference-free — it uses only `PairLJCut` and quadrature, so it directly proves the sampler draws ∝ exp(−βE). Box side is chosen so the histogram/analytic support covers the full r² e^{−βV} tail without truncation; T is chosen so the well is meaningfully populated (not effectively hard-sphere, not fully repulsive); equilibration steps are discarded before accumulation; the χ²/KS and ⟨V⟩ tolerances are set from block SE so a longer run tightens, never hides, the bar. This test may run a bit longer than a typical ms unit test but stays a normal ctest.

Edge cases:
- Acceptance-ratio-vs-δ trend. A tuning heuristic, not a theorem: assert only that the acceptance ratio is monotonically decreasing across an increasing δ ladder and that a mid δ lands in a broad ~20–60% band. Fast.
- Stationarity guard inside the NIST gate: verify ⟨U⟩ shows no drift across blocks before averaging, so a non-equilibrated run fails loudly rather than reporting a biased mean.

Domain validation (`scientific`, gated behind `ATV_RUN_E2E_TESTS`, registered in the e2e tail of `tests/CMakeLists.txt`):
- NIST SRSW LJ ⟨U*⟩/N at the two state points above, reduced units (ε=σ=1). The reference is defined for the **pure-cut** potential at rc=3σ **plus** U_LRC/N; the test therefore samples the fluid, re-evaluates the pure-cut (unshifted) energy on stored snapshots via `truncated_lj_energy`, and adds the analytic U_LRC/N — it must NOT compare a cut-and-shifted energy against the cut+LRC reference. Tolerance is taken from the run's own block SE: |⟨U*⟩/N − ref| ≤ max(3·block SE, 0.05 absolute floor), with the block SE reported.

Honesty notes embedded in the criteria (not cosmetic):
- The two-particle P(r) gate (test 1) is the load-bearing correctness proof and runs in normal CI. The NIST gate (test 3) is a longer stochastic end-to-end check that MUST NOT silently truncate — if run length or N is reduced for CI budget, the criterion logs the achieved block SE and the tolerance actually used; a shorter run widens the reported tolerance rather than passing on a hidden margin.
- `PairLJCut` is cut-**and-shifted** (`V_shift = V(rc)`), so the sampler generates the cut-and-shifted ensemble while the SRSW reference is the pure-cut ensemble + LRC. At rc=3σ the per-pair shift is ≈ 5.5e-3 ε and, at liquid density, its per-atom contribution is O(0.1–0.3) in U*/N — i.e. NOT negligible. The test removes it by re-evaluating the pure-cut energy with `truncated_lj_energy`; the residual difference between the two *ensembles* (configurations generated under cut-and-shifted vs pure-cut) is a documented second-order caveat listed alongside finite-size O(1/N) bias (NIST N=500) and statistical error, so the tolerance is defensible.

## Out of scope

- Any new sampler feature or production abstraction: new move types, local/incremental ΔE energy, neighbor-aware acceptance, GCMC/REMD/HMC, CUDA MC. This sub-spec is test-only.
- New factory functions or config builders (project rule; helpers are header-local test utilities, not `ATV::` API).
- Modifying `cpu::mc::MC`, `TranslateMove`, or `PairLJCut` behaviour — if a validation failure implicates production code, it is reported back to sub-spec 02, not patched here.
- Relocating the NIST gate into `benchmarks/` — it stays a gated e2e ctest until an MC validation benchmark home exists (mirrors the existing `test_mpi_pe_nvt_multirank_n8` treatment).

---
slug: mc-sampling-03-validate
criteria:
  - id: ac-001
    summary: Two-particle P(r) matches analytic r^2*exp(-bV) and <V> matches quadrature
    type: runtime
    pass_when: |
      test_cpu_mc_pair_distribution passes: for a fixed+mobile PairLJCut pair
      in a box large enough not to truncate the r^2*e^{-bV} tail, at a T that
      meaningfully samples the well, after equilibration discard the sampled
      separation histogram matches analytic P(r) ∝ r^2*exp(-V(r)/k_BT) by a
      chi^2/KS statistic within the stated tolerance, AND sampled <V> equals
      the 1-D quadrature integral within block standard error. Reference-free.
    status: pending
  - id: ac-002
    summary: Acceptance ratio decreases monotonically in delta; mid delta in 20-60% band
    type: runtime
    pass_when: |
      test_cpu_mc_acceptance_ratio passes: across an increasing delta ladder the
      measured acceptance ratio is monotonically decreasing, and at least one
      mid delta yields an acceptance ratio within the broad ~20-60% band.
      Trend + loose band only.
    status: pending
  - id: ac-003
    summary: NIST LJ <U*>/N at two state points within max(3*blockSE, 0.05), truncated+LRC
    type: scientific
    evaluator_hint: mol:bench
    pass_when: |
      test_cpu_mc_lj_nist (ATV_RUN_E2E_TESTS=ON) passes at BOTH state points
      (T*=0.85,rho*=0.78 -> -5.5386; T*=0.90,rho*=0.90 -> -6.1773): energy is
      the pure-cut (unshifted) truncated energy re-evaluated at rc=3 sigma PLUS
      analytic U_LRC/N; <U> stationarity (no drift) is verified after
      equilibration discard; |<U*>/N - ref| <= max(3*block SE, 0.05); and the
      achieved block SE, the applied tolerance, and the cut-and-shifted-vs-cut
      + finite-size (N=500, O(1/N)) caveats are logged (no hidden pass).
    status: pending
---

# Acceptance criteria

Criteria are ranked by correctness value. **ac-001 is the highest-value gate** — it is the reference-free proof that the sampler draws ∝ exp(−βE), and it runs in normal CI. ac-002 is a fast tuning-trend guard. ac-003 is the literature-consistency end-to-end gate.

## ac-001 (runtime, highest value) — two-particle P(r) exactness
Binding proof of Boltzmann sampling using only `PairLJCut` and deterministic quadrature: histogram-vs-analytic P(r) ∝ r²·e^{−βV} within a block-SE-derived χ²/KS tolerance, and sampled ⟨V⟩ equal to the quadrature ⟨V⟩ within block SE. No external data, so a failure is unambiguously a sampler defect.

## ac-002 (runtime) — acceptance-ratio-vs-δ trend
Monotone decrease of acceptance ratio in δ plus a mid-δ landing in a broad ~20–60% band. A heuristic sanity check on the proposal, deliberately asserted as trend + loose band, not an exact acceptance value.

## ac-003 (scientific) — NIST SRSW LJ ⟨U*⟩/N
Per the evaluator protocol, a `scientific`-typed criterion parks this spec at **code-complete** once ac-001/ac-002 are `verified` and the e2e target builds; final `verified`/`failed` on ac-003 is written only by `/mol:bench`, never by `/mol:impl`. The gate reconstructs the pure-cut energy (via the test-scope `truncated_lj_energy` helper) + analytic LRC to match the SRSW definition, sets its tolerance from the run's own block SE with a 0.05 absolute floor, checks stationarity, and logs SE + tolerance + the cut-and-shifted-vs-cut and finite-size caveats so the pass is auditable and not cosmetically tightened.

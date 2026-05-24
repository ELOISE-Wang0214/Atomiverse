---
slug: cuda-md-review-cleanup-09-bench-and-regression
criteria:
  - id: ac-001
    summary: bench_md_step_pe_npt benchmark exists and runs
    type: runtime
    pass_when: |
      `cmake --build build --target bench_md_step` builds clean. Running
      `./build/benchmarks/bench_md_step --benchmark_filter=pe_npt` produces
      output for at least 3 N values (e.g., 10000, 30000, 50000) with a
      mean step time reported. The benchmark uses
      `cuda::NoseHooverNPT` + full PE topology (bonds+angles+dihedrals).
    status: pending
    last_checked: ""
  - id: ac-002
    summary: test_launch_purity runtime gate passes
    type: runtime
    pass_when: |
      `ctest --test-dir build -R '^test_launch_purity$' --output-on-failure`
      exits 0. The test drives a 100-step NPT MD loop with the full force
      suite and asserts `cuda::detail::sync_count() == 0` after every step.
    status: pending
    last_checked: ""
  - id: ac-003
    summary: test_pair_lj_no_d2h static grep gate passes
    type: runtime
    pass_when: |
      `ctest --test-dir build -R '^test_pair_lj_no_d2h$' --output-on-failure`
      exits 0. The test reads `src/cuda/potentials/pair_lj_cut.cu` and
      asserts zero occurrences of `cudaStreamSynchronize` and zero
      occurrences of `cudaMemcpy.*HostFromDevice|cudaMemcpy.*DeviceToHost`.
    status: pending
    last_checked: ""
  - id: ac-004
    summary: test_driver_no_per_phase_sync runtime gate passes
    type: runtime
    pass_when: |
      `ctest --test-dir build -R '^test_driver_no_per_phase_sync$' --output-on-failure`
      exits 0. A 50-step `cuda::md::MD` NVE loop reports zero
      `cuda::stream_synchronize` calls during the hot path.
    status: pending
    last_checked: ""
  - id: ac-005
    summary: test_neighborlist_skin_gating runtime gate passes (registered)
    type: runtime
    pass_when: |
      `ctest --test-dir build -R '^test_neighborlist_skin_gating$' --output-on-failure`
      exits 0. Created in link 02; this acceptance criterion confirms the
      test is registered in `tests/CMakeLists.txt` and continues to pass.
    status: pending
    last_checked: ""
  - id: ac-006
    summary: Umbrella /mol:review re-run shows zero 🚨, ≤3 🔴, ≤8 🟡 residual
    type: review
    pass_when: |
      Running `/mol:review --axis=perf,arch` over the cuda-md scope (the
      directories enumerated in `cuda-md-review-cleanup.md` Scope section)
      yields a verdict line with severity distribution where:
        - critical (🚨) count == 0
        - high (🔴) count ≤ 3
        - medium (🟡) count ≤ 8
        - low (🟢) any count
      Each remaining 🔴 / 🟡 finding is explicitly enumerated below as
      `backlog-finding-NN` with file:line and rationale-for-deferral.
    status: pending
    last_checked: ""
  - id: ac-007
    summary: GH200 PE-NPT benchmark — cluster 2 (launch-purity) speedup ≥ 1.20×
    type: runtime
    pass_when: |
      `scripts/bench_md_step_perf.sh` (or equivalent operator-run measurement)
      reports a mean-step-time speedup of ≥ 1.20× for the
      cluster-2-attributable subset of changes (measured by checkpointing
      git refs at the cluster boundary, or by feature-flag toggle if available).
      Recorded in a build-log artifact `bench_results_pe_npt.json`.
    status: pending
    last_checked: ""
  - id: ac-008
    summary: GH200 PE-NPT benchmark — cluster 3 (skin-gating) speedup ≥ 1.30×
    type: runtime
    pass_when: |
      As ac-007 but for cluster 3 (`cuda-md-review-cleanup-02-nlist-skin-gating`).
      Skin-gating alone should free ~30-50% of step time per the review's
      structural estimate; ≥ 1.30× total step speedup is the gate.
    status: pending
    last_checked: ""
  - id: ac-009
    summary: GH200 PE-NPT benchmark — clusters 5+6+7 combined speedup ≥ 1.10×
    type: runtime
    pass_when: |
      Combined attribution to clusters 5 (reductions+atomics), 6 (launch
      config+streams), and 7 (microfixes) shows ≥ 1.10× step-time speedup
      on the PE-NPT benchmark, recorded in `bench_results_pe_npt.json`.
    status: pending
    last_checked: ""
  - id: ac-010
    summary: GH200 PE-NPT benchmark — other clusters no regression (≥ 0.95×)
    type: runtime
    pass_when: |
      Clusters 1 (component-capability), 4 (Box by value), 8 (placement)
      show step-time ratio ≥ 0.95× (no regression). Placement-only changes
      are not expected to move perf; gate is no-regression only.
    status: pending
    last_checked: ""
  - id: ac-011
    summary: Inventory spec status updated to done with completion summary
    type: code
    pass_when: |
      `.claude/specs/cuda-md-review-cleanup.md` frontmatter `status:` field
      changed from `reference` to `done`. A new "Completion" section at the
      bottom summarizes:
        - Each of the 9 chain links with verified ✓ marker
        - Per-cluster perf delta from `bench_results_pe_npt.json`
        - Enumerated residual backlog (≤3 🔴 + ≤8 🟡) with file:line
    status: pending
    last_checked: ""
  # ---- Backlog enumeration (to be filled when ac-006 runs) ----
  # The following are PLACEHOLDERS for the residual backlog enumeration
  # required by ac-006. They will be populated when the umbrella review
  # actually runs. Until link 09 implementation, these IDs are reserved.
  - id: ac-backlog-residual-list
    summary: Residual finding enumeration (filled at link 09 implementation)
    type: code
    pass_when: |
      This acceptance file (when finalized at link 09 implementation) lists
      each residual 🔴 / 🟡 finding from the umbrella `/mol:review` re-run
      with file:line + rationale-for-deferral + suggested follow-up spec slug.
      Placeholder until the implementation phase fills it in.
    status: pending
    last_checked: ""

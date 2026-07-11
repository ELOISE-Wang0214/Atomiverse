# Specs index

| slug | status | summary |
|---|---|---|
| [mc-sampling-03-validate](mc-sampling-03-validate.md) | code-complete | Chain 3/3 — MC physics validation (benchmarks/mc/lj/validation). ac-001 (reference-free two-particle P(r) ∝ r²e^(−βV) + ⟨V⟩-vs-quadrature) and ac-002 (acceptance-ratio-vs-δ trend) are verified. **ac-003 (NIST SRSW LJ ⟨U*⟩/N) is the open item**: ρ*=0.78 matches (dev 0.0096) but ρ*=0.90 needs a high-T melt of the FCC seed before the target-T average — melt pass added, rerun pending. |
| [gpw-edl-aimd](gpw-edl-aimd.md) | approved (blocked) | GPW production bridge for the 1214-atom ZnO electric-double-layer AIMD. **Cannot be completed as written**: 3 of its 5 open criteria (ac-006 analytic forces, ac-010 capstone NVT, ac-013 DFT-D3) are CP2K-parity gates, and calling CP2K is forbidden in this tree. They also sit on genuine open physics — GPW analytic forces are incomplete for GGA / d-projectors / periodic images (test_cpu_gpw_forces: ZnO 0.534 vs 5e-5 tol), and the ZnO diagonalization SCF does not converge (test_cpu_gpw_ot). Retained as the tracking doc for that physics campaign, not as an in-flight spec. |

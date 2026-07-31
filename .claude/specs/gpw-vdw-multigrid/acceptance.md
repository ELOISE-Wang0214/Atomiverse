---
spec: gpw-vdw-multigrid
created: 2026-07-31
criteria:
  - id: ac-001
    summary: SC nonlocal vdW potential applied on the multigrid path
    type: code
    pass_when: |
      build_H_and_energy adds nonlocal_vdw_on_grid V_nl to the fine-grid veff
      regardless of use_prod; the "multigrid not supported" throw for
      vdw_df/vdw_df2/rvv10 is gone from the Γ RKS runner.
    status: passed
  - id: ac-002
    summary: Multigrid vs single-grid internal consistency
    type: runtime
    pass_when: |
      reg_gpw-vdw-multigrid runs the same dimer with ngrids=1 and
      ngrids=4/rel_cutoff=40 and |E_tot(mg) − E_tot(sg)| passes the measured
      tight gate documented in REFERENCE.md (1e-5 Ha class or better).
    status: passed
  - id: ac-003
    summary: E_disp on multigrid matches offline CP2K NGRIDS 4 pin
    type: runtime
    pass_when: |
      |E_disp(mg) − CP2K 2026.1 DRSLL NGRIDS 4 Dispersion energy| ≤ tolerance
      stated in REFERENCE.md (1e-4 Ha class); pin recorded with version,
      input, and output line; ctest never launches CP2K.
    status: passed
  - id: ac-004
    summary: UKS + nonlocal vdW rejects fail-fast (found rot closed)
    type: runtime
    pass_when: |
      A unit test shows Γ UKS (n_alpha/n_beta) with vdw_df/vdw_df2/rvv10
      throws with a how-to instead of silently dropping E_c^nl.
    status: passed
  - id: ac-005
    summary: Queue and roadmap updated
    type: docs
    pass_when: |
      docs/developer/cp2k-method-queue.md moves gpw-vdw-multigrid to done;
      docs/roadmap/dft.md drops the "multigrid still throws" caveat.
    status: passed
out_of_scope:
  - k-point vdW (gpw-vdw-kpoints)
  - UKS vdW capability (reject only)
  - vdW × multigrid forces/stress gates
  - CUDA
---

# Acceptance — `gpw-vdw-multigrid`

Done = Γ RKS GPW with `vdw_df`/`vdw_df2`/`rvv10` runs self-consistently on
the CP2K-style REL_CUTOFF multigrid (`ngrids>1`), internally consistent with
the single grid and pinned to offline CP2K 2026.1 NGRIDS 4 DRSLL.

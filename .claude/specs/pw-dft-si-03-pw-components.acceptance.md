---
slug: pw-dft-si-03-pw-components
criteria:
  - id: ac-001
    summary: PW_Poisson matches analytic Hartree for Gaussian charge
    type: scientific
    pass_when: |
      Gaussian ρ(r)=exp(-αr²), V_H from FFT Poisson matches analytic erf solution, rel error<1e-4.
    status: pending
  - id: ac-002
    summary: PW_LDA uniform electron gas matches VWN5 reference
    type: scientific
    pass_when: |
      ρ=0.01 e/Bohr³ → E_xc matches VWN5 reference within 1e-10 Hartree/electron.
    status: pending
  - id: ac-003
    summary: PW_Hamiltonian kinetic term correct at k=0
    type: scientific
    pass_when: |
      Single plane-wave test: eigenvalue = ½|G|² within 1e-12.
    status: pending
  - id: ac-004
    summary: PW_Davidson converges N_occ eigenpairs
    type: scientific
    pass_when: |
      Free-electron H, N_occ=8, residual<1e-8 within 50 iterations per k-point.
    status: pending
  - id: ac-005
    summary: Pulay mixing ≥2x better than simple mixing
    type: scientific
    pass_when: |
      Random initial + fixed target density, Pulay residual ≤ 0.5× simple mixing residual after 10 iters.
    status: pending
  - id: ac-006
    summary: PW_Convergence correct at tol boundary
    type: code
    pass_when: |
      tol=1e-8: diff=1e-7→not converged, diff=9e-9→converged. Warmup=2 suppresses first 2 checks.
    status: pending
  - id: ac-007
    summary: All six components share single PWWorkspace after build()
    type: code
    pass_when: |
      After driver.build(), all components' workspace pointer equal; array dims match n_G + n_grid.
    status: pending
  - id: ac-008
    summary: Full SCF converges for 2-atom Si in <20 iterations
    type: runtime
    pass_when: |
      2-atom Si, E_cut=10Ry, Gamma-only → ‖ρ_out−ρ_in‖<1e-8 within 20 SCF steps.
    status: pending
  - id: ac-009
    summary: No allocation in launch() across all components
    type: code
    pass_when: |
      Memory monitor shows zero allocations during launch() calls post-build().
    status: pending
---

---
spec: gaussian-electrode-charge
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Offline literature/MetalWalls-kernel pin recorded
    type: docs
    pass_when: |
      regressions/md/gaussian-electrode-charge/REFERENCE.md states
      oracle: literature (or metalwalls with published-kernel version),
      oracle_version, geometry, η, units (eV/Å), fixed-Q energy pin,
      toy Q Gaussian vs point-charge, and documents point-charge A
      approximation. MD-only build flags noted.
    status: passed
  - id: ac-002
    summary: Public Gaussian kernel + pair API without Config bags
    type: code
    pass_when: |
      cpu::GaussianChargeKernel and cpu::PairGaussianCoul exist under
      public include path; parameters are ctor args / members; no new
      public *Config; no make_*/build_* factories; not a new Driver category.
    status: passed
  - id: ac-003
    summary: Fixed-Q energy matches literature kernel pin
    type: runtime
    pass_when: |
      Regression asserts |E - E_pin| < 1e-9 eV for the 5-atom parallel-plate
      toy with published Gaussian kernel (include self-energy).
    status: passed
  - id: ac-004
    summary: Toy Q Gaussian vs point-charge A
    type: runtime
    pass_when: |
      Regression solves 2-site neutrality-constrained A Q = Ψ for Gaussian
      and point-charge A; |Q_gauss - Q_pin| < 1e-9 e; |Q_gauss - Q_point| > 0
      documented (point-charge A is approximate for electrodes).
    status: passed
  - id: ac-005
    summary: Unit tests for kernel + forces
    type: runtime
    pass_when: |
      tests/cpu/potentials/test_gaussian_charge*.cpp (or test_pair_gaussian_coul)
      pass: G/self limits, dimer E/F analytic, FD force consistency.
    status: passed
  - id: ac-006
    summary: Dual perf gate
    type: runtime
    pass_when: |
      Fair wall vs MW/LAMMPS ELECTRODE when binary available (ratio ≤ 10),
      OR perf_exempt with reason in REFERENCE.md (expected: no twin binary).
    status: passed
  - id: ac-007
    summary: MD-only build
    type: runtime
    pass_when: |
      Configure with ATV_BUILD_MD=ON KSPACE=ON TESTS=ON REGRESSIONS=ON and
      HF/PW/GPW/CUDA/Torch OFF; atomiverse + unit + regression green.
    status: passed
out_of_scope:
  - Full constant-potential BO Component (E1)
  - 2D Ewald / slab kspace
  - CUDA / polarizable coupling
---

# Acceptance — gaussian-electrode-charge

Done means electrode Gaussian charge kernels are public on the MD surface,
match the published Reed/MetalWalls analytic pin for energy and a toy \(Q\)
solve, and document why point-charge \(A\) is only an approximation.

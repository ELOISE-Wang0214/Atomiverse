---
spec: gpw-gaussian-plane-wave
title: Acceptance criteria — GPW (Gaussian and Plane Waves) DFT
---

# Acceptance — GPW CP2K-parity path

Reference numbers come from the CP2K repo at
`/nobackup/proj/disk/teoroo/personal/jicli594/work/cp2k` (binary at
`~/.local/bin/cp2k`, v2026.1). Regtest refs live in the matching
`TEST_FILES.toml`; tutorial refs are produced by running the deck once.

---

ac-001: GTH potential + MOLOPT basis parsers
  type: structural
  status: pending
  test_file: tests/cpu/gpw/test_cpu_gpw_cp2k_io.cpp
  description: >
    Parse GTH potentials (data/POTENTIAL: GTH-PADE-q8/q1; data/GTH_POTENTIALS:
    GTH-PBE-q6/q1) and GTH basis sets (data/BASIS_SET: DZVP-GTH-PADE;
    data/BASIS_MOLOPT: DZVP-MOLOPT-GTH) for Ar, H, O. Parsed rloc, C_i, r_l, h^l_ij
    and contracted GTO exponents/coefficients match the source files exactly
    (hand-checked fixtures).

ac-002: Uniform-grid collocation conserves charge and round-trips
  type: scientific
  status: pending
  test_file: tests/cpu/gpw/test_cpu_gpw_collocate.cpp
  description: >
    Density collocated from a known density matrix onto the uniform FFT grid
    integrates to N_electrons within grid tolerance; integrate(collocate(φ_μφ_ν))
    reproduces the analytic overlap S_μν within the CUTOFF/REL_CUTOFF discretization
    floor.

ac-003: Analytic GTH integrals match a reference single atom
  type: scientific
  status: pending
  test_file: tests/cpu/gpw/test_cpu_gpw_gth_integrals.cpp
  description: >
    Local short-range and nonlocal separable projector matrix elements for an
    isolated atom agree with CP2K's per-atom PP matrix (or an independently derived
    analytic value) to ≤ 1e-8.

ac-004: M1 — periodic LDA GPW matches CP2K Ar.inp and H2.inp
  type: scientific
  status: pending
  test_file: tests/cpu/gpw/test_cpu_gpw_m1_lda_energy.cpp
  description: >
    With CP2K's exact basis + GTH-PADE potentials + matching CUTOFF, GPW total
    energy matches CP2K tests/QS/regtest-gpw-1/Ar.inp (-21.04944231395054 Ha) and
    H2.inp (-1.06345574964568 Ha) to |ΔE| ≤ 1e-6 Ha.

ac-005: M2 — periodic PBE (GGA) GPW matches a CP2K PBE regtest
  type: scientific
  status: pending
  test_file: tests/cpu/gpw/test_cpu_gpw_m2_pbe_energy.cpp
  description: >
    GGA density-gradient path on the uniform grid; total energy of a periodic PBE
    GPW regtest input matches its CP2K reference to |ΔE| ≤ 1e-6 Ha (native XC) /
    ≤ 1e-7 Ha (runtime libxc).

ac-006: M3 — non-periodic MT Poisson reproduces the H2O getting-started tutorial
  type: scientific
  status: pending
  test_file: tests/cpu/gpw/test_cpu_gpw_m3_h2o_tutorial.cpp
  description: >
    The getting-started docs/getting-started/h2o.inp (isolated H2O, PBE,
    DZVP-MOLOPT-GTH, GTH-PBE-q6/q1, CUTOFF 400/REL_CUTOFF 50, PSOLVER MT) run via
    GPWConfig reproduces the CP2K total energy (from a local cp2k run of the deck)
    to |ΔE| ≤ 1e-5 Ha; SCF converges at EPS_SCF 1e-6 within MAX_SCF 50.

ac-007: MT Poisson correctness
  type: scientific
  status: pending
  test_file: tests/cpu/gpw/test_cpu_gpw_mt_poisson.cpp
  description: >
    Martyna-Tuckerman solver on a single Gaussian charge reproduces the analytic
    isolated-system potential/energy and is size-independent of box padding beyond
    the decoupling threshold.

ac-008: GPWEvaluator plugs into Evaluate / MD / Minimize unchanged
  type: structural
  status: pending
  test_file: tests/cpu/gpw/test_cpu_gpw_evaluator_driver.cpp
  description: >
    GPWEvaluator (role Evaluator, capabilities Energy|Force|DensityMatrix|Orbitals)
    runs under cpu::single_point, md::MD (force phase), and minimizer::Minimize with
    no driver changes, writing frame.scalars.energy (eV) and atom forces.

ac-009: M4 — analytic GPW forces match finite difference and CP2K
  type: scientific
  status: pending
  test_file: tests/cpu/gpw/test_cpu_gpw_forces.cpp
  description: >
    Analytic forces agree with the finite-difference oracle to ≤ 1e-5 Ha/bohr and
    with CP2K forces on the tutorial geometry to ≤ 1e-4 Ha/bohr; a short geom-opt of
    H2O reduces the max force below the tutorial threshold.

ac-010: M5 — tutorial + benchmark migrated under benchmarks/gpw with CP2K parity
  type: performance
  status: pending
  test_file: benchmarks/gpw/h2o/validation/bench_gpw_h2o_cp2k.cpp
  description: >
    benchmarks/gpw/h2o/validation mirrors pw/si/validation: a CP2K deck template +
    run script generate a reference CSV; the Atomiverse benchmark asserts total
    energy within tolerance and emits its own CSV. SOURCES.md documents provenance
    (CP2K version, basis/PP, libxc). examples/09_gpw.cpp runs the H2O single point.

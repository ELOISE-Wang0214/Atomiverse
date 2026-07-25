---
spec: rsh-eri-cam-b3lyp
created: 2026-07-24
criteria:
  - id: ac-001
    summary: rys_eri_erf_block ω→0 matches bare Coulomb ERI
    type: runtime
    pass_when: |
      tests/cpu/scf/test_rys_eri_erf.cpp asserts max |I_erf(ω=0) − I_bare|
      ≤ 1e-12 on covered shell classes (at least s and p).
    status: pending
  - id: ac-002
    summary: rys_eri_erf_block known (ss|ss) pin hard-coded
    type: scientific
    pass_when: |
      Unit test embeds a hard-coded (ss|ss) erf-Coulomb value (with offline
      source comment: tool, version, date) and matches within 1e-10 Ha.
    status: pending
  - id: ac-003
    summary: CAM-B3LYP RSH params are (0.19, 0.46, 0.33)
    type: code
    pass_when: |
      rsh_exchange_params(XCType::CAMB3LYP) returns alpha=0.19, beta=0.46,
      mu=0.33; non-RSH types return beta=0 and mu=0.
    status: pending
  - id: ac-004
    summary: Fock uses bare+erf mix for CAM (not single alpha only)
    type: code
    pass_when: |
      hamiltonian_builder (or equivalent) multiplies K[1/r] by alpha and
      K[erf(μr)/r] by beta when rsh_exchange_params(type).beta != 0; pure
      B3LYP path still uses hybrid_exchange_fraction only.
    status: pending
  - id: ac-005
    summary: cam-b3lyp-h2-sto3g regression matches PySCF pin
    type: runtime
    pass_when: |
      regressions/scf/cam-b3lyp-h2-sto3g (or documented twin slug) builds under
      ATV_BUILD_HF, runs without live PySCF, and |E_atv − E_pin| ≤ 1e-6 Ha
      for the REFERENCE.md offline pin.
    status: pending
  - id: ac-006
    summary: No public Config / factory for CAM path
    type: code
    pass_when: |
      Product uses XCType::CAMB3LYP + RKS members/ctors only; grep on touched
      public headers finds no *Config bag and no make_* factory for CAM.
    status: pending
  - id: ac-007
    summary: Queue and roadmap updated when done
    type: docs
    pass_when: |
      docs/developer/pyscf-method-queue.md moves cam-b3lyp out of blocked into
      done with regression path; docs/roadmap/dft.md marks CAM-B3LYP CPU ✓.
    status: pending
out_of_scope:
  - Analytic CAM gradients / Hessian
  - UKS/ROKS CAM and RSH TDDFT kernels
  - ωB97X / LC-ωPBE / HSE product functionals (engine reuse only)
  - CUDA attenuated ERI
  - Live third-party oracles in ctest
---

# Acceptance — `rsh-eri-cam-b3lyp`

**Done** means: (1) an erf-attenuated four-center ERI exists and is
identity-consistent with the bare Coulomb engine at \(\omega=0\); (2)
CAM-B3LYP is a first-class `XCType` with published \((\alpha,\beta,\mu)\);
(3) closed-shell RKS energy on a textbook STO-3G case matches an offline
PySCF pin to \(10^{-6}\) Ha; (4) queue/roadmap reflect the foothold.

## AC-001 — ω→0 identity

Regression risk is silent drift of B3LYP if the bare path is shared
incorrectly. Gate: ω=0 attenuated block == existing `rys_eri_block`.

## AC-002 — hard-coded integral pin

At least one non-trivial erf integral golden, not only the identity gate.

## AC-003 — parameter table

CAM constants must match PySCF/libxc CAM-B3LYP defaults used for the energy pin.

## AC-004 — Fock mix

Method correctness lives in the exchange operator mix, not only the enum.

## AC-005 — e2e energy

Single public-API regression under `regressions/scf/`; goldens in-tree.

## AC-006 — iron laws

No public Config bags / factories on the new path.

## AC-007 — harness close-out

Queue row and roadmap stay consistent with product truth.

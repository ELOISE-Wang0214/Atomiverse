---
spec: periodic-rhf-gamma
created: 2026-07-24
criteria:
  - id: ac-001
    summary: Periodic exchange lattice sum is cutoff-stable
    type: runtime
    pass_when: |
      tests/cpu/scf/test_periodic_exchange_lattice.cpp shows max |K(cutoff) −
      K(cutoff+1 shell)| below the documented tol for the toy system.
    status: pending
  - id: ac-002
    summary: Molecular RHF path unchanged when periodic=false
    type: runtime
    pass_when: |
      Existing molecular RHF energy regression(s) still pass; periodic flag
      default is false and does not alter bare RHF Fock assembly.
    status: pending
  - id: ac-003
    summary: Public API is members/ctor only (no Config bag)
    type: code
    pass_when: |
      Product surface is RHF (or PeriodicRHF) members such as periodic=true;
      no public *Config and no make_* factory for periodic HF.
    status: pending
  - id: ac-004
    summary: periodic-rhf-gamma regression matches offline PySCF pin
    type: runtime
    pass_when: |
      regressions/scf/periodic-rhf-gamma builds under ATV_BUILD_HF, never
      launches PySCF, and |E_atv − E_pin| ≤ tol stated in REFERENCE.md
      (≤ 1e-5 Ha foothold unless REFERENCE tightens).
    status: pending
  - id: ac-005
    summary: REFERENCE documents electrostatic twin settings
    type: docs
    pass_when: |
      REFERENCE.md records PySCF version, cell, basis, and the J/K/Enn
      convention (Ewald vs truncated) used for the pin so a re-pin is
      reproducible.
    status: pending
  - id: ac-006
    summary: Queue and roadmap mark Γ-point periodic RHF done
    type: docs
    pass_when: |
      docs/developer/pyscf-method-queue.md moves periodic-rhf-gamma to done;
      docs/roadmap/hartree-fock.md marks Γ-point periodic RHF CPU ✓.
    status: pending
out_of_scope:
  - k-points / complex Bloch orbitals
  - Periodic UHF and hybrid PBC DFT
  - Analytic forces / stress
  - CUDA
  - Live third-party oracles in ctest
---

# Acceptance — `periodic-rhf-gamma`

**Done** means a Γ-point closed-shell periodic RHF energy is available on the
public Frame path, numerically twinned to offline PySCF PBC HF for one
minimal orthorhombic case, without regressing molecular RHF, and with
harness/roadmap updated.

## AC-001 — cutoff stability

Guards silent incompleteness of the lattice sum for K.

## AC-002 — molecular non-regression

PBC work must not fork or poison the default molecular engine.

## AC-003 — API iron laws

Constructor/members only; no Config bags.

## AC-004 — e2e pin

One `regressions/scf/periodic-rhf-gamma` textbook foothold.

## AC-005 — twin documentation

Electrostatics convention is part of the scientific contract.

## AC-006 — close-out

Queue blocked row and roadmap ○→✓ stay honest.

---
spec: esmm-contract-05-accumulate-switch
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Co-bound sources sum instead of clobbering
    type: runtime
    pass_when: |
      U7 unit test passes: Evaluate with a += stub potential and a nested
      RHF-shaped source yields kEnergy == E_stub + E_es and
      kEnergyElectronic == E_es; order-independence checked by swapping
      binding order.
    status: pending
  - id: ac-002
    summary: All five ES evaluators are += on accumulator channels
    type: code
    pass_when: |
      grep over src/cpu/{gpw,planewave,semiempirical,scf} shows no in-plan
      `=` writes to kEnergy/kAtom_f*/kStress remaining except documented
      standalone entries (MFM run/forces/mp2_forces/ccsd_forces) and state
      channels (kAtom_charge, kStressKinetic); GPW's 7 kEnergy + 7 force +
      1 stress sites and PW/AM1/RM1 sites all flipped.
    status: pending
  - id: ac-003
    summary: kEnergyElectronic registered, pre-inserted, zeroed, written
    type: runtime
    pass_when: |
      U8 passes (ZeroEnergy clears it); Frame(n) pre-inserts it; MFM
      (nested), GPW, PW, AM1, RM1, DFTB2 all += their electronic energy into
      it; the zero.cpp invariant comment names the Bias term-less debt.
    status: pending
  - id: ac-004
    summary: Sole-source behavior unchanged
    type: runtime
    pass_when: |
      Full ctest green and the AIMD/BOMD/CPMD/XL regression subset passes
      with unchanged numbers (Zero*-then-+= ≡ = for sole sources).
    status: pending
  - id: ac-005
    summary: R1 pins the contract with two literature values
    type: runtime
    pass_when: |
      regressions/architecture/lj_rhf_coupled builds under ATV_BUILD_HF,
      runs ≤ 60 s, asserts kEnergy == E_LJ(analytic) + (−1.117506 Ha→eV)
      to 1e-6 relative and kEnergyElectronic equals the RHF term alone;
      REFERENCE.md cites Szabo & Ostlund Table 3.11 and the LJ source.
    status: pending
out_of_scope:
  - Bias per-term key
  - QM/MM coupling internals (Wave 2), CPMD branch restructuring (Wave 3)
  - Test-local overwrite stubs
---

# Acceptance — `esmm-contract-05-accumulate-switch`

**Done** means any combination of energy/force/stress sources — classical,
electronic, or both — composes additively in every driver, the electronic
share is observable on `kEnergyElectronic`, and no sole-source number in the
suite moved.

## AC-001 — the connection, demonstrated

The first assertion in the tree that two heterogeneous sources sum.

## AC-002 — no silent survivors

A missed overwrite site passes every sole-source test; only the grep gate
catches it.

## AC-003 — decomposition invariant heals

kEnergy == Σ per-term (+ named Bias gap), now including the electronic share.

## AC-004 — equivalence proof

The switch is invisible until a second source appears.

## AC-005 — literature-pinned witness

Both addends pinned independently; the sum pins `+=`.

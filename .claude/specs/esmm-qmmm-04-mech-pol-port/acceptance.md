---
spec: esmm-qmmm-04-mech-pol-port
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Mechanical requires an engine; DFTB2 hard default is gone
    type: runtime
    pass_when: |
      QmmmMechanical::build() throws without an engine; grep shows no
      hard-wired DFTB2 member construction inside the coupling; AM1 and
      DFTB2 both run through set_engine (embedding capability never probed);
      QM region resolves from kAtomQmZ.
    status: pending
  - id: ac-002
    summary: Mechanical classical cross-Coulomb from kAtom_charge
    type: runtime
    pass_when: |
      Unit test pins the QM–MM classical Coulomb energy + forces against the
      closed-form two-charge sum with QM charges taken from kAtom_charge;
      no qm_charges vector exists (grep); no zero-charge validation on the
      mechanical path.
    status: pending
  - id: ac-003
    summary: Polarizable permanent field is density-derived
    type: runtime
    pass_when: |
      The permanent QM field is built from the engine-written kAtom_charge
      on the private QM frame; an engine that leaves the column all-zero
      triggers a typed throw; μ = αE single-site Applequist limit passes
      against the closed-form dipole.
    status: pending
  - id: ac-004
    summary: Both couplings ride the shared machinery and the contract
    type: code
    pass_when: |
      QmmmMechanical and QmmmPolarizable each hold internal::QmSubsystem and
      internal::QmMmVdw by value (no duplicated region/link/scatter logic —
      grep for a second tag→row resolve); both += into kEnergy/kAtom_f* and
      write the QM share into kEnergyElectronic.
    status: pending
  - id: ac-005
    summary: Three regressions upgraded to quantitative pins
    type: runtime
    pass_when: |
      qmmm_h2_ar pins E_QM to the h2_dftb2 literature value + closed-form
      12-6 cross-LJ; qmmm_pol_h2o pins the Applequist single-site closed
      form; qmmm_link_ch3oh passes rigid-translation invariance ≤ 1e-9 eV
      and |Σ F| ≤ 1e-9 eV/Å including link projection; all ≤ 60 s, public
      API only.
    status: pending
out_of_scope:
  - MD conservation closure and docs (esmm-qmmm-05)
  - New polarizable physics; GEEP / PME / ONIOM; MPI; CUDA
---

# Acceptance — `esmm-qmmm-04-mech-pol-port`

**Done** means all three couplings share one region/link/scatter/vdW machinery
and one engine slot; mechanical embedding is engine-agnostic with honest
classical cross-Coulomb; polarizable embedding polarizes against the QM
density's own charges; and the three shipped regressions assert numbers, not
signs.

## AC-001 — engine-agnostic mechanical
The last hard-coded engine in the tree dies here.

## AC-002 — mechanical electrostatics
kAtom_charge is the single charge channel — no parallel vectors.

## AC-003 — polarizable upgrade
Density-derived permanent field replaces bare-Z; loud when absent.

## AC-004 — shared machinery
Three couplings, one implementation of the common 80%.

## AC-005 — quantitative pins
The qualitative-gate debt of the original multiscale wave is repaid.

---
spec: esmm-cpmd-04-ao-cpmd
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Fixed-state gradient parity at a converged determinant
    type: runtime
    pass_when: |
      tests/cpu/scf/test_scf_cpmd.cpp shows ScfGradient::compute_at_density
      (W = occ·C(CᵀFC)Cᵀ) matches ScfGradient::compute element-wise to 1e-10
      at a converged determinant, restricted AND unrestricted; the extracted
      mean_field_energy_at_density leaves DensityStage numerics bit-identical
      (existing SCF suites green).
    status: pending
  - id: ac-002
    summary: Production CPMD steps run zero SCF iterations
    type: runtime
    pass_when: |
      In all three new regressions, every post-seed MD step reports
      iterations() == 0; the seed step is one converged full SCF (throws if
      unconverged).
    status: pending
  - id: ac-003
    summary: Conserved-energy drift pinned for RHF, RKS, UHF
    type: runtime
    pass_when: |
      regressions/dynamics/{cpmd-rhf-h2,cpmd-rks-water,cpmd-uhf-o2} pass:
      relative kConserved drift within the envelope frozen in each
      REFERENCE.md; kConserved == kEnergy+kKinetic+T_e(+E_therm) consistency
      holds; each REFERENCE.md cites Car & Parrinello PRL 55, 2471 (1985)
      and Marx & Hutter (CUP 2009) §2.4-2.5 and records the tuned mu/dt.
    status: pending
  - id: ac-004
    summary: Refusal lattice is complete and guard-tested
    type: runtime
    pass_when: |
      Typed throws with one EXPECT_THROW each for: ROHF, ROKS, GHF,
      predictor(), xl_bomd, ri != None, x2c, solvent_eps > 0, mom, soscf,
      non-empty point_charges, nonzero external_field, d3bj, d4,
      n_dropped_ > 0 at seed; standalone run()/forces() with cpmd=true throw
      std::logic_error.
    status: pending
  - id: ac-005
    summary: dt arrives only via TimestepConsumer
    type: code
    pass_when: |
      MeanFieldMethod has no public dt member; timestep_required() returns
      cpmd; build_on throws when cpmd is set and no dt was pushed; grep shows
      no local dt duplication on the SCF CPMD path.
    status: pending
  - id: ac-006
    summary: d3bj/d4 force-gap verification recorded
    type: docs
    pass_when: |
      The close-out summary states whether ScfGradient lacks dispersion force
      terms while e_total_ includes the energy; if confirmed, a /mol:fix
      ticket note exists naming the non-conservative BOMD+d3bj debt (fix
      itself explicitly out of scope here).
    status: pending
  - id: ac-007
    summary: Determinant staleness documented
    type: docs
    pass_when: |
      MeanFieldMethod cpmd Doxygen and the tutorial page state that
      determinant() returns the seed-SCF wavefunction on production steps.
    status: pending
out_of_scope:
  - XL-BOMD on SCF, PW CPMD, OT-CP, multi-k, CPMD stress
  - QM/MM x CPMD (point_charges gate stays)
  - Fixing the d3bj/d4 gradient gap
---

# Acceptance — `esmm-cpmd-04-ao-cpmd`

**Done** means Car–Parrinello dynamics runs on RHF/RKS/UHF/UKS through the
nested MD force slot with zero per-step SCF iterations, conserved-energy
drift pinned for all three spin shapes, a complete typed refusal lattice,
and dt single-sourced from the integrator.

## AC-001 — the two new primitives are correct

Fixed-density energy and multiplier-form forces are the physics core; parity
at convergence is their ground truth.

## AC-002 — it is actually CPMD

Zero production SCF iterations is the defining signature (and the cpmd-water
precedent pin).

## AC-003 — E_cons is the physics gate

Drift envelopes + the CP citation make the regression scientific, not a
golden snapshot.

## AC-004 — loud unsupported paths

Every combination not productized must throw, not silently degrade.

## AC-005 — no duplicated clock

The esmm-cpmd-02 protocol is the only dt channel.

## AC-006 — iron law on discovered debt

The dispersion-force gap is named and ticketed, never silently ignored.

## AC-007 — honest observability

Observers must know the determinant is the seed state.

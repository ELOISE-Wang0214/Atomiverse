---
slug: esmm-qmmm-04-mech-pol-port
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: multiscale / qmmm (mechanical + polarizable)
queue: esmm-campaign
deferred: false
depends_on: [esmm-qmmm-03-electrostatic-rewrite]
unblocks: [qmmm-mech-pol-engine-slot]
grilled: true
---

# Spec: `esmm-qmmm-04-mech-pol-port` — Mechanical + Polarizable on the shared machinery

**Wave:** esmm Wave 2 (QM/MM redesign), stage 4 of 5.
**Status:** approved — grilled via plan-mode design review + operator Q&A, 2026-07-26.

## Summary

Port the other two couplings onto `internal::QmSubsystem` + `internal::QmMmVdw`
(from `esmm-qmmm-03`), keeping three concrete classes (distinct embedding semantics =
distinct classes; no enum-morphing bag). `QmmmMechanical` **loses its hard-wired DFTB2
member** — the engine slot is required (any `QmEngine`, embedding capability never
probed) — and reads the QM region from the `kAtomQmZ` mask column. Mechanical
embedding's classical QM–MM Coulomb reads the QM atoms' charges **straight from
`kAtom_charge`** (in mechanical embedding QM atoms legitimately keep their MM charges;
the zero-charge invariant applies ONLY to Electrostatic/Polarizable). `QmmmPolarizable`
gets the engine slot plus a physics upgrade: its **permanent QM field is built from
engine-provided per-atom charges** (the engine writes `kAtom_charge` on the private QM
frame — Mulliken for the Gaussian family, AM1 charges for NDDO) instead of bare
nuclear Z; it throws when the engine provides none (all-zero column). The three
existing multiscale regressions upgrade from qualitative gates to quantitative pins.

## Domain basis

- Mechanical embedding: E = E_QM(vacuum) + E_MM + E_QM–MM^{LJ} + E_QM–MM^{Coul,classical}
  with fixed classical charges on both sides (Senn & Thiel 2009 §2.2). The classical
  cross-Coulomb term and its forces are closed-form pair sums (minimum image,
  lj_cutoff-consistent cutoff policy documented).
- Polarizable embedding: Applequist/Thole induced dipoles μ_i = α_i E_i(permanent +
  induced), solved self-consistently by the existing public `solve_induced_dipoles`
  free function; E_pol = −½ Σ μ·E^perm (Thole, Chem. Phys. 59, 341 (1981); Thompson &
  Schenter, J. Phys. Chem. 99, 6374 (1995) — both already cited in the header).
  Permanent field now sourced from density-derived charges.

Units: eV, Å, e, eV/Å, α [Å³].

## Design

| Class | Ctor (binding shape) | Notes |
| ----- | -------------------- | ----- |
| `QmmmMechanical : Potential` | `(std::vector<LJTypeCoeff> cross_lj, double lj_cutoff)` | engine required at build (any `QmEngine`; `as_embedding()` never probed); mask column region; classical QM–MM Coulomb from `kAtom_charge` (QM side may be nonzero — **no** zero-charge validation here); `links()`, `link_scale()`, `reset_history_each_launch()`, `set_engine/set_engine_ref` as in stage 3 |
| `QmmmPolarizable : Potential` | `(std::vector<LJTypeCoeff> cross_lj, double lj_cutoff, double embed_cutoff)` + members `alpha` (per-MM-type [Å³] table), `thole_a`, `pol_tol` [e·Å], `pol_max_iter` | electrostatic-embedding base flow (zero-charge invariant, `as_embedding()` required) + induced-dipole solve; permanent QM field from the private frame's engine-written `kAtom_charge`; throws when that column is all-zero after evaluate |

Both hold `internal::QmSubsystem` + `internal::QmMmVdw` **by value** (shared machinery
= member classes, not a base — recorded in stage 3). Both accumulate (`+=`) into
`kEnergy`/`kAtom_f*` and write the QM share into `kEnergyElectronic`. DFTB2 remains a
first-class engine choice — `mech.set_engine(cpu::DFTB2{...})` — it is only the
*hard default* that dies.

## Files to create or modify

- rewrite: `cpu/multiscale/qmmm_mechanical.{h,cpp}`, `cpu/multiscale/qmmm_polarizable.{h,cpp}`
- `src/cpu/multiscale/` CMake rows if Mechanical's gate moves (DFTB2 default gone ⇒
  Mechanical itself needs only the engine the user supplies; keep the DFTB2-using
  regression under `ATV_BUILD_SEMIEMPIRICAL`)
- tests: update `tests/cpu/multiscale/test_qmmm_mechanical.cpp`,
  `test_qmmm_polarizable.cpp` (new ctors, engine slot, new throws)
- regressions upgraded in place: `regressions/multiscale/qmmm_h2_ar/`,
  `qmmm_pol_h2o/`, `qmmm_link_ch3oh/` (+ their `REFERENCE.md`)

## Tasks

- [ ] RED: mechanical unit tests — engine required (build throws without);
  mask-column region; classical QM–MM Coulomb equals the closed-form two-charge sum;
  QM charges read from `kAtom_charge` (no qm_charges vector — grep); AM1 and DFTB2
  both accepted (embedding never probed)
- [ ] RED: polarizable unit tests — permanent field from engine charges (hand-built
  frame where Mulliken charges are known); all-zero engine charge column throws;
  μ = αE single-site Applequist limit against the closed-form dipole in the analytic
  field of two point charges
- [ ] Port both couplings onto `QmSubsystem`/`QmMmVdw`
- [ ] Upgrade `qmmm_h2_ar`: E_QM pinned to the existing `semiempirical/h2_dftb2`
  literature value + cross-LJ evaluated in closed form from the 12-6 formula
  (arithmetic, not a snapshot)
- [ ] Upgrade `qmmm_pol_h2o`: single-site Applequist μ = αE closed-form pin
- [ ] Upgrade `qmmm_link_ch3oh`: rigid-body translation invariance |ΔE| ≤ 1e-9 eV and
  |Σ F| ≤ 1e-9 eV/Å **including the link-force projection** (keeps Senn & Thiel cite)
- [ ] GREEN: full suite; Doxygen with units on changed public surfaces

## Testing strategy

Unit: the RED batches (mirror paths). Regressions: the three upgrades — each replaces
a qualitative finite/sign gate with a quantitative closed-form or literature pin; no
new regression case is added (no new capability — the port is the feature, its
witnesses are the upgraded pins).

## Out of scope

- MD conservation closure + docs (`esmm-qmmm-05`)
- New polarizable physics (damped fields beyond existing Thole, anisotropic α)
- GEEP / PME / ONIOM / adaptive; CUDA; MPI mode

## Notes for implementers

- **HARD GATE:** before touching `src/cpu/scf/*` or any CMakeLists, verify
  `git status` shows them clean (operator's ~760-file libcint/C++20 work is
  uncommitted as of 2026-07-26). If dirty: STOP and report. `cpu/multiscale/` is clean.
- Breaking change by design: every mechanical-embedding user must now install an
  engine (house rule: no hard-wired default engine, no back-compat shim). Update all
  in-tree consumers in the same change.
- Polarizable's induced-dipole solver (`solve_induced_dipoles`) keeps its public
  free-function status (documented, literature-cited exception).

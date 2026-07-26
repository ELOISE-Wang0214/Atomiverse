---
slug: esmm-cpmd-04-ao-cpmd
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: electronic-propagation / cpmd
queue: esmm-campaign
deferred: false
depends_on: [esmm-cpmd-01-electronic-home, esmm-cpmd-02-timestep-protocol, esmm-cpmd-03-observability]
unblocks: [esmm-cpmd-05]
grilled: true
---

# Spec: `esmm-cpmd-04-ao-cpmd` — Car–Parrinello MD on the Gaussian SCF methods

**Slug:** `esmm-cpmd-04-ao-cpmd`
**Domain:** Gaussian-basis AO CPMD (CPU)
**Status:** approved (grilled via plan-mode design review + operator Q&A, 2026-07-26)

## Summary

CPMD ships only for GPW. `electronic::Cpmd::advance(H, S, dt, occ)` needs
exactly the matrices `MeanFieldMethod` already holds as members (`F_alpha_`,
`S_`), so AO-basis Car–Parrinello on RHF/RKS/UHF/UKS is a wiring problem plus
two missing primitives: a **fixed-density evaluation entry** (build F at a
GIVEN P without SCF iterations) and **forces at a non-eigenstate C** (the
full Lagrange-multiplier Pulay weight). Deliver both, mirror GPW's member
surface (minus `cpmd_dt` — dt arrives via `TimestepConsumer`), and pin three
regressions. Car & Parrinello, PRL **55**, 2471 (1985); Marx & Hutter,
*Ab Initio Molecular Dynamics* (CUP 2009) §2.4–2.5.

## Domain basis

Extended Lagrangian: μ C̈ = −(I − S C Cᵀ)(∂E/∂C) with holonomic CᵀSC = I
(electronic RATTLE — already inside `electronic::Cpmd`). Forces on nuclei at
the propagated (non-stationary) state use the multiplier form of the Pulay
weight: **W_σ = occ·C_σ (C_σᵀ F_σ C_σ) C_σᵀ** (restricted: W = 2·C(CᵀFC)Cᵀ),
which reduces to Σ_i ε_i c_i c_iᵀ when C are eigenvectors of F — exactly what
GPW ships in `gpw_forces_at_density`. Energy at fixed density uses the same
four-branch formula `DensityStage` evaluates each SCF iteration.

## Design

### Public members on `MeanFieldMethod<Derived>` (mirror GPW minus dt)

`bool cpmd = false; double cpmd_mu = 400.0; bool cpmd_thermostat = false;
double cpmd_ekin_target_ha = 1e-5; double cpmd_tau_e = 50.0;` + readbacks
`last_electronic_ke_ev()`, `last_electronic_thermostat_ev()`. Doxygen states:
nested-MD-force-slot only; RHF/RKS (occ=2) and UHF/UKS (dual channel, occ=1);
`run()`/`forces()` with cpmd → `std::logic_error`.
`MeanFieldMethod` implements `TimestepConsumer` (`timestep_required() ⇔ cpmd`;
`set_timestep` validates finite > 0 [fs]); reached in MD's push walk via
`dynamic_cast` on the nested `DriverNode` (esmm-cpmd-02 machinery).
Privates: `std::optional<electronic::Cpmd> cpmd_alpha_, cpmd_beta_;
double cpmd_dt_fs_ = 0, last_te_ev_ = 0, last_etherm_ev_ = 0;` +
`validate_cpmd_config_() / cpmd_seed_() / cpmd_step_(Frame&)`.

### Entry — `run_nested_erased` branch only

- `cpmd && propagator ready` → `cpmd_step_(frame)`, return 0 (the
  production-SCF-iters==0 pin).
- `cpmd && !ready` → **cold seed**: normal `drive_scf_` full SCF with forces;
  throw if `!converged()`; `cpmd_seed_()` copies the occupied columns of
  `C_alpha_` (S-orthonormal by construction, canonical-orth path included)
  into `electronic::Cpmd` (+β for UHF/UKS).
- `build_on` gains the dt gate: `cpmd && !(cpmd_dt_fs_ > 0)` → throw
  ("cpmd requires an MD integrator publishing the timestep").

### `cpmd_step_` — fixed-density evaluation, NOT through the stage Plan

`DriveOpts.max_iter = 0` cannot express this (the loop body never runs, so no
Fock is built), and the stage pipeline (DIIS, diagonalize, density-rebuild)
is exactly what CPMD must skip. Use `internal::HamiltonianBuilder` directly —
the precedent is the de-shift epilogue (`mean_field_method.cpp:901-918`):

1. `PositionGuard` unwrap; rebuild `*integrals_sp_ = IntegralPack{frame,
   basis, charge, compute_eri, ecp, ghost_atoms, false}`; `resolve_occupation`.
2. `cpmd_alpha_->align_to_overlap(ints.S)` (+β) — Löwdin transport under the
   new geometry's metric.
3. `P_a = occ·C_a C_aᵀ` (restricted: the factor-2 total on the α channel,
   matching the loop convention); `P_b` for β.
4. `HamiltonianBuilder hb; hb.prepare(...); FockResult fock = hb.step(ints,
   restricted, is_hf, xc_type_, hybrid_alpha_override, P_a, P_b)`.
5. `e_total_ = internal::mean_field_energy_at_density(...)` — extracted from
   `DensityStage::launch` (second user ⇒ extraction allowed; DensityStage
   calls it with unchanged conventions, **bit-identical**).
6. Forces: `ScfGradient::compute_at_density(frame, ints, C_occ_a, F_a,
   C_occ_b, F_b, e_ha)`; write `kAtom_f*` + `kEnergy` via `+=` — correct
   independent of esmm-contract-05, because MD binds ZeroForces/ZeroEnergy at
   pre_force. (When esmm-contract-05 lands, its ES-site sweep adds the
   `kEnergyElectronic` per-term write at this site; do not reference that key
   here while it is absent from `frame_schema.h`.)
7. `cpmd_alpha_->advance(F_alpha, S, cpmd_dt_fs_, occ)` (+β, occ=1); publish
   `kEnergyElectronKinetic` / `kEnergyElectronThermostat` (insert-guarded);
   `iterations_ = 0`.

`determinant()` is NOT updated on production steps (propagated C has no
canonical ε) — stays the seed-SCF wavefunction; documented limitation.

### `ScfGradient::compute_at_density` (new, `scf/internal/gradient.h`)

Signature per the campaign plan; internally the assembly is refactored to
accept `(P_total, P_a, P_b, W)` directly, with `compute()` forwarding from
the determinant as today. Parity unit test: at a converged determinant,
`compute_at_density == compute` to 1e-10, restricted **and** unrestricted.

### Refusal lattice (`validate_cpmd_config_`, all typed throws + guard tests)

ROHF/ROKS (shared spatial orbitals ↛ independent channels), GHF (spinor CP is
new physics), `use_predictor_` (no ASPC-CP hybrid; mirrors GPW exclusivity),
`xl_bomd` (esmm-cpmd-05 sibling, mutually exclusive), `ri != RI::None`,
`x2c`, `solvent_eps > 0`, `mom != Off`, `soscf`, non-empty `point_charges`
(**the QM/MM×CPMD seam — lifting this gate later is the entire integration**),
nonzero `external_field`, `d3bj`/`d4` (dispersion forces are absent from
`ScfGradient`, so E_cons would silently drift), `n_dropped_ > 0` at seed
(near-singular S ⇒ ill-conditioned RATTLE metric).

### Named-debt verification (first task; do NOT fix here)

Verify the suspicion that `ScfGradient` has **no d3bj/d4 force terms while
`e_total_` includes the dispersion energy** — if confirmed, existing
BOMD+d3bj forces are non-conservative (pre-existing). Record the finding in
the close-out summary and open a separate `/mol:fix` ticket note; the CPMD
refusal of d3bj/d4 keeps E_cons honest either way.

## Files to create or modify

- `include/atomiverse/cpu/scf/mean_field_method.h`, `src/cpu/scf/mean_field_method.cpp`
- `include/atomiverse/cpu/scf/internal/gradient.h`, `src/cpu/scf/internal/gradient.cpp`
- `src/cpu/scf/stages/density_stage.cpp` (+ the extracted energy helper's home
  under `scf/internal/`)
- NEW `tests/cpu/scf/test_scf_cpmd.cpp`
- NEW `regressions/dynamics/cpmd-rhf-h2/`, `cpmd-rks-water/`, `cpmd-uhf-o2/`
  (each `reg_*.cpp` + `REFERENCE.md`)
- `regressions/CMakeLists.txt`, `tests/CMakeLists.txt`
- NEW `docs/tutorial/methods/cpmd-rhf-h2.md` (+ index row)

## Tasks

- [ ] Verify the d3bj/d4 force gap in `ScfGradient`; record + ticket note
- [ ] RED: `test_scf_cpmd.cpp` — full refusal lattice (one EXPECT_THROW per
      row above) + standalone `run()`/`forces()` throw + parity test shape
- [ ] Extract `mean_field_energy_at_density` from `DensityStage` (bit-identical;
      DensityStage keeps its old-P/new-P call conventions)
- [ ] Implement `compute_at_density` with W_σ = occ·C_σ(C_σᵀF_σC_σ)C_σᵀ;
      parity green at 1e-10 (restricted + unrestricted)
- [ ] Members + TimestepConsumer + `validate_cpmd_config_` + `cpmd_seed_` +
      `cpmd_step_` + `run_nested_erased` branch + `build_on` dt gate
- [ ] Tune μ/dt empirically on H2/water/O2 (AO residual is stiffer than PW —
      GPW water needed μ=2000 at dt=0.1 fs); then pin the three regressions
- [ ] Tutorial page + roadmap row; scoped build + ctest + regressions green

## Testing strategy

**Unit** (`test_scf_cpmd.cpp`): refusal lattice; standalone-entry throws;
`compute_at_density` parity; seed-then-step returns `iterations() == 0`;
meta keys present after a step. **Regressions** (each ≤60 s, public API,
`ATV_BUILD_HF AND ATV_BUILD_MD`):

- `dynamics/cpmd-rhf-h2` — RHF/STO-3G H2, VV dt (protocol-sourced), NVE;
  pins: production SCF iters == 0; relative `kConserved` drift ≤ envelope
  (set from the first stable run, cpmd-water style); kConserved ==
  PE+KE+T_e consistency; bond-length sanity. Cite CP 1985 + Marx & Hutter.
- `dynamics/cpmd-rks-water` — RKS/LSDA STO-3G water; exercises the KS energy
  branch + XC gradient at fixed P; drift + iters==0.
- `dynamics/cpmd-uhf-o2` — UHF triplet O2/STO-3G; dual-channel occ=1; drift +
  iters==0 + ⟨S²⟩ sanity.

## Out of scope

- XL-BOMD on SCF (esmm-cpmd-05); PW CPMD; OT-CP; multi-k anything; CPMD
  stress; QM/MM×CPMD (the point_charges gate stays shut here)
- Fixing the d3bj/d4 force gap (verify + ticket only)

## Notes for implementers

- **HARD GATE:** before touching `src/cpu/gpw/gpw_evaluator.{h,cpp}`,
  `src/cpu/scf/mean_field_method.{h,cpp}`, `scf/internal/gradient.*`, or any
  CMakeLists, verify `git status` shows them clean. As of 2026-07-26 the
  operator's libcint-removal/C++20 work (~760 files) is uncommitted and
  includes them. If dirty: STOP and report.
- Drift tolerances are **envelope pins**, not literature numbers — tune
  before freezing, then cite the conservation criterion (Marx & Hutter) for
  the *form* of the gate.
- Operator ordering (2026-07-26): the esmm-cpmd wave runs BEFORE Wave 1's
  esmm-contract-05. The `+=` in cpmd_step_ does not depend on it (MD zeroes
  at pre_force); `DriveOpts::accumulate` / `kEnergyElectronic` integration
  arrives with contract-05, which will rebase over this wave's
  mean_field_method / gpw_evaluator changes — do not block on it here.

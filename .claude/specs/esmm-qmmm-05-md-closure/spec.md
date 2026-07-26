---
slug: esmm-qmmm-05-md-closure
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: multiscale / qmmm (MD closure + docs)
queue: esmm-campaign
deferred: false
depends_on: [esmm-qmmm-04-mech-pol-port]
unblocks: [qmmm-molecular-dynamics]
grilled: true
---

# Spec: `esmm-qmmm-05-md-closure` — QM/MM molecular dynamics closes; docs tell the new story

**Wave:** esmm Wave 2 (QM/MM redesign), stage 5 of 5.
**Status:** approved — grilled via plan-mode design review + operator Q&A, 2026-07-26.

## Summary

With the third-law site forces (stage 1), the engine slot (stages 2–3), and the
contract citizenship in place, QM/MM MD is finally *correct* — momentum-exact and
energy-conserving. This spec ships the conservation regression that closes the roadmap
row **"QM/MM molecular dynamics" (○ → ✓)** in `docs/roadmap/multiscale.md`, rewrites
the multiscale narrative docs onto the new API (mask column, engine slot, MM-side
force-field composition), and documents the Monte-Carlo posture
(`reset_history_each_launch`). No new library code beyond what the regression needs.

## Domain basis

Conservation as the correctness invariant of MD (Frenkel & Smit, *Understanding
Molecular Simulation*, §4.3): with exact Newton's-third-law pairs, total momentum is
conserved to machine precision and NVE energy drift scales as O(dt²) for velocity
Verlet. For QM/MM MD specifically (Senn & Thiel 2009 §3), momentum conservation is
the direct witness that the MM sites receive the reaction of the embedding operator —
the exact defect the old coupling had. SCF-convergence noise enters the energy gate,
so the drift tolerance must sit above `tol`-induced force noise (engine `tol = 1e-9`
Ha recommended for the case).

## Design — the regression

`regressions/multiscale/qmmm_md_conservation/` (`ATV_BUILD_HF AND ATV_BUILD_MD`):

- System: one QM water (RHF/STO-3G via the mask column: qm_z = {8,1,1}) + a handful
  of charged LJ point-MM atoms in a cubic box; `embed_cutoff` ≥ half the box diagonal
  (no cutoff-crossing discontinuity — the sharp-cutoff limitation is documented, not
  exercised).
- Composition is the 自由组合 witness: `cpu::md::MD` + `VelocityVerlet` +
  `PairLJCut` (MM–MM) + `CoulLong`-free tiny cell (MM–MM Coulomb via the pair styles
  the case defines) + `QmmmElectrostatic` with `set_engine(RHF{...}.predictor(3))`.
- Gates (all from meta/spans through the public API):
  1. |Σ **p**| < 1e-8 amu·Å/fs at every sampled step over 200 NVE steps (exact third
     law — no drift trend);
  2. |ΔE_total|/N < 1e-4 eV/atom over the run;
  3. halving dt halves-squared the drift (O(dt²) Verlet scaling, factor 4 ± 1);
  4. warm start engaged: production SCF iteration counts drop after step 1.
- `REFERENCE.md` cites Frenkel & Smit §4.3 and Senn & Thiel 2009 §3 and records the
  tolerances' provenance (conservation identities, not snapshots); `perf_exempt`
  with reason (no external QM/MM-MD twin wall on the worker host).

## Files to create or modify

- new: `regressions/multiscale/qmmm_md_conservation/{reg_qmmm_md_conservation.cpp,REFERENCE.md}`
- `regressions/CMakeLists.txt` (gate `ATV_BUILD_HF AND ATV_BUILD_MD`)
- `docs/roadmap/multiscale.md` — "QM/MM molecular dynamics" ○ → ✓; rewrite the
  "Textbook QM/MM" paragraph (engine slot, mask column, MM-side ownership, third law)
- `docs/tutorial/ch20-hf-properties-environment.md` — `EmbeddingSite`,
  `embed_self_energy`, `embedding_site_forces` on the SCF surface
- `docs/tutorial/ch19-hf-large-systems.md` — QM/MM section onto the new coupling API
  (mask column + engine slot + ordinary MM styles end-to-end example)
- MC posture note where the tutorial covers sampling (`reset_history_each_launch`,
  and the esmm-contract-04 ASPC-under-MC guard rationale)

## Tasks

- [ ] RED-equivalent: write the regression with the four gates; verify it FAILS
  against a build with the old coupling behavior stubs disabled is not possible —
  instead verify gate 1 fails when site-force scatter is artificially zeroed in a
  local scratch build (sanity that the gate has teeth), then restore
- [ ] Tune dt / n_steps / engine tol so the case runs ≤ 60 s single-threaded
- [ ] `REFERENCE.md` with citations + tolerance provenance + perf_exempt reason
- [ ] Register in `regressions/CMakeLists.txt`
- [ ] Roadmap flip + Textbook-QM/MM paragraph rewrite
- [ ] Tutorial ch19/ch20 updates (end-to-end example = the usage block from
  esmm-qmmm-03, extended with MD)
- [ ] GREEN: `regressions` target green; docs build clean

## Testing strategy

This spec's deliverable IS the e2e case (one feature — "QM/MM MD" — one regression).
Unit coverage came with stages 1–4. The teeth-check for gate 1 is a one-off local
experiment, not a committed test.

## Out of scope

- QM/MM in NPT / with thermostats beyond the NVE witness (compose freely; not pinned)
- QM/MM-MC regression (MC posture is documented; a sampling pin belongs to a future
  wave if wanted)
- QM/MM × CPMD (seam: lift the point_charges gate in the CPMD specs)
- PME-QM/MM, switching functions, GEEP, ONIOM, adaptive QM/MM; CUDA

## Notes for implementers

- **HARD GATE:** `regressions/CMakeLists.txt` and docs files are in the operator's
  uncommitted ~760-file sweep as of 2026-07-26 — verify clean before editing; if
  dirty: STOP and report.
- Momentum gate: sample Σp every 10 steps; assert max over samples. Use double
  accumulation (AccReal) when summing.
- Keep the box small enough that `embed_cutoff ≥ half diagonal` holds with the
  10–12 Å cutoffs used in the docs example, or shrink cutoffs consistently.
- The dt-halving leg re-runs the same trajectory at dt/2 — budget both runs inside
  the 60 s envelope (short n_steps for the scaling leg is fine; document).

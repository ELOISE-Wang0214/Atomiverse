---
title: Fix NoseHooverNPT runaway box expansion (Trotter-order investigation)
slug: nose-hoover-npt-trotter-fix
status: code-complete
created: 2026-05-20
last-updated: 2026-05-20
parent: msu-pe-validation.md
---

# Fix NoseHooverNPT runaway box expansion (Trotter-order investigation)

## Summary

`tests/cpu/md/test_lammps_npt_parity.cpp` (currently `DISABLED_`) shows that
`cpu::md::NoseHooverNPT` fails the LAMMPS `fix npt iso` parity gate on a
256-atom LJ-argon system at T=100 K, P=0, τ_P=1000 fs, dt=1 fs:

| Quantity | LAMMPS | Atomiverse | Δ |
|---|---|---|---|
| ⟨T⟩      | 99.965 K | 93.5 K | 6.4% (acceptable) |
| ⟨V⟩      | 14250 Å³ | 2.5 × 10⁶ Å³ | 17 000% (runaway) |

A 2026-05-20 session pinned the failure mode precisely: for the first ~3000
steps the system settles around V ≈ 13500 Å³ (close to LAMMPS), then a small
sustained positive ⟨ε̇⟩ ≈ 3 × 10⁻⁴ /fs slowly inflates the box. Once V exceeds
the LJ cutoff range (~50000 Å³ for 256 atoms at σ=3.405 Å, cutoff=8.5 Å) the
attractive pair virial vanishes, nothing balances the kinetic-pressure drive,
and V exponentiates without bound.

Line-by-line comparison with `src/fix_nh.cpp` from LAMMPS develop branch
(github.com/lammps/lammps) showed every closed-form expression in
`NoseHooverNPT` is algorithmically equivalent to LAMMPS:

* pressure drive   `f_ω = V·(P_int − P_ext)/W`         ✓
* MTK kinetic coupling `(tdof·k·T)/(pdim·N·W)`         ✓
* chain restoring  `(kecurrent − lkt_press)/etap_mass`  ✓
* velocity scaling `exp(−dt/2 · ω̇·(1 + 1/N))`           ✓
* box affine       `exp(dt · ω̇)`                       ✓

Yet LAMMPS converges on the same input. Spec scope: find the missing
ingredient and bring Atomiverse to parity.

## Domain basis

Standard MTK NPT (Martyna–Tobias–Klein, *J. Chem. Phys.* 1994 101, 4177;
Tuckerman 2006 *J. Phys. A* 39, 5629; Shinoda 2004 *Phys. Rev. B* 69, 134103)
defines an extended-Lagrangian Hamiltonian that is exactly conserved under
the symmetric Trotter splitting

    e^{i L_NHC_atoms · dt/2} · e^{i L_NHC_baro · dt/2}
       · e^{i L_eps_dot · dt/2} · e^{i L_v_baro · dt/2} · e^{i L_v_force · dt/2}
       · e^{i L_x_box · dt}
       · e^{i L_v_force · dt/2} · e^{i L_v_baro · dt/2} · e^{i L_eps_dot · dt/2}
       · e^{i L_NHC_baro · dt/2} · e^{i L_NHC_atoms · dt/2}

Both Atomiverse and LAMMPS implement this scheme. They differ in how the
central `e^{i L_x_box · dt}` operator is split:

| Engine | Central block |
|---|---|
| Atomiverse | one fused `propagate_x_and_box(dt)` — analytical solution `x = x·exp(ω̇·dt) + v·dt·sinhc(ω̇·dt/2)`; box scaled in one step. |
| LAMMPS | `remap(dt/2) → nve_x(dt) → remap(dt/2)` — atoms affine-rescaled with box in two half-steps that bracket the velocity-only position kick. |

For small ω̇·dt both forms agree to O((ω̇·dt)²). For our case ω̇·dt ≈ 3×10⁻⁴,
the per-step difference is ~10⁻⁸. Over 10⁴ steps the cumulative discrepancy
is ~10⁻⁴ — far too small to explain the observed 17 000% volume divergence.
That suggests the bug is *not* in the central block itself, but somewhere
in how Atomiverse stitches the chain integrator and the omega_dot update to
that block. The next-tier suspects, in priority order:

1. **Trotter ordering of `nhc_atoms` vs `nhc_baro` in the outer pair.**
   LAMMPS calls `nhc_press_integrate` *before* `nhc_temp_integrate` in
   `initial_integrate` (line 920–926) and again *after* `nhc_temp_integrate`
   in `final_integrate` (line 1012–1013). Atomiverse calls `nhc_atoms_half`
   *before* `nhc_baro_half` in `initial_integrate` (P1, P2) and reverses to
   `nhc_baro_half` → `nhc_atoms_half` in `end_of_step` (P10, P11). Both
   are symmetric Trotter splits, but the order matters for energy
   conservation: LAMMPS = baro–atoms / atoms–baro, Atomiverse = atoms–baro
   / baro–atoms. Swap to match.

2. **Stale kinetic stress in `update_eps_dot` at `end_of_step` P9.**
   LAMMPS's `final_integrate` runs `compute_pressure` *after* `nh_v_press`
   (line 998) so the pressure used in `nh_omega_dot` reflects the freshly-
   rescaled velocities. Atomiverse's P9 uses `frame.stress[0..2]` as written
   by `VirialStress` at `post_force` time — using velocities that pre-date
   both the final-half VV kick and the `scale_v_by_baro` operation.
   Investigation in 2026-05-20 noted that `KE_global_` is tracked through
   `scale_v_by_baro` correctly, but the kinetic stress trace stored in
   `frame.stress` is not. The delta is small (≈ (dt/2)·F/m magnitude on
   each component), but over thousands of steps could create the residual
   sustained ⟨ε̇⟩.

3. **Box affine vs sinhc subtlety.** LAMMPS's two half-step remaps with
   intermediate `nve_x` give an effective `x_new = x_0·exp(ω̇·dt) +
   v·dt·exp(ω̇·dt/2)`. Atomiverse's fused form gives `x_0·exp(ω̇·dt) +
   v·dt·sinhc(ω̇·dt/2)`. For positive ω̇, `exp(ω̇·dt/2) > sinhc(ω̇·dt/2)`,
   so LAMMPS gives a slightly *larger* per-step position update. Over many
   steps the integrated trajectories differ; whether this could pump
   energy into ε̇ is unclear without a Hamiltonian-conservation check.

4. use LAMMPS execution order exactly for whole MD loop: https://docs.lammps.org/Developer_flow.html

The session also ruled out (with controlled experiments) the following
non-causes: sign-convention flip (mathematically equivalent), zeroing the
`(6/N_dof)·KE` coupling (no change), `drag = 2`/`drag = 20` (no
improvement), the `3·W·ε̇²` chain restoring force (Atomiverse's `W` already
absorbs the factor of 3 relative to LAMMPS's per-axis `omega_mass`), and
larger `τ_P = 5000 fs` (reduces but does not eliminate drift).

## Design

This spec is an **investigation + targeted fix**, not a redesign of NPT
from scratch. The work splits into three phases.

### Phase A — Trotter ordering swap (most likely fix)

Swap the two outer NHC calls so the order matches LAMMPS:

| File | Change |
|---|---|
| `src/cpu/md/nose_hoover_npt.cpp` `launch()` | In `initial_integrate`: swap P1↔P2 (call `nhc_baro_half` *before* `nhc_atoms_half`). In `end_of_step`: swap P10↔P11 (call `nhc_atoms_half` *before* `nhc_baro_half`). |
| `include/atomiverse/cpu/md/nose_hoover_npt.h` docstring | Update the P1–P11 trace to reflect the new order. |

Expected impact: if Trotter symmetry breaking is the cause, parity test
volume divergence should drop from ~17 000 % into the 5 % gate.

### Phase B — Fresh kinetic stress in `update_eps_dot` (if Phase A insufficient)

If Phase A doesn't close the gap, address the stale-kinetic-stress
hypothesis by recomputing the kinetic-stress trace from current velocities
inside `update_eps_dot` rather than reading `frame.stress[0..2]`.

Cleanest implementation: split frame.stress accumulation so VirialStress
writes only the kinetic component into a new `frame.stress_kinetic[6]`
field (or `Frame::Stress`-style helper struct), leaving `frame.stress[0..6]`
as the pair-virial-only sum. Then:

```cpp
// in update_eps_dot at end_of_step phase:
const double trace_virial = f.stress[0] + f.stress[1] + f.stress[2];
const double trace_kinetic_fresh = 2.0 * KE_global_;  // KE_global_ is freshly
                                                       // tracked through P8.
const double trace_total = trace_virial + trace_kinetic_fresh;
const double P_int = trace_total / (3.0 * V);
```

Affected components:

* `src/cpu/md/virial_stress.{h,cpp}` — write to `stress_kinetic`, not `stress`.
* `src/cpu/md/barostat.cpp` (Berendsen) — read `stress + stress_kinetic`.
* `src/cpu/md/nose_hoover_npt.cpp` `update_eps_dot` — compute kinetic fresh.
* Any observer that reads stress (RDF, thermo) — combine the two fields.

Scope estimate: 8–12 files. Mechanically simple but touches the public Frame
data model. Defer to Phase B only if Phase A isn't enough.

### Phase C — Hamiltonian conservation diagnostic (fallback / validation)

Instrument the integrator with a per-step `H_NPT_extended` computation
(atom KE + atom PE + ½ W·ε̇² + Σ ½ Q_t·η_t_dot² + Σ ½ Q_b·η_b_dot² +
N_dof·kT·η_t[0] + kT·Σ η_t[k≥1] + kT·Σ η_b + p_target·V) and verify it's
conserved to single-precision floating point over 100 000 steps. If the
parity test passes but conservation drifts, expose it as a regression gate.
This is the canonical MTK correctness diagnostic from Tuckerman 2006 §III.

## Files to create or modify

Phase A (likely the whole spec):

```
src/cpu/md/nose_hoover_npt.cpp            (swap P1↔P2, P10↔P11)
include/atomiverse/cpu/md/nose_hoover_npt.h (docstring update)
tests/cpu/md/test_lammps_npt_parity.cpp   (rename DISABLED_ → Slow_; tighten
                                           ⟨V⟩ tolerance to 5%; update
                                           docstring findings)
```

If Phase A fails and Phase B is needed:

```
include/atomiverse/cpu/frame.h             (add stress_kinetic field)
src/cpu/md/virial_stress.cpp               (target new field)
src/cpu/md/barostat.cpp                    (read combined)
src/cpu/md/nose_hoover_npt.cpp             (recompute kinetic fresh)
+ ~5 observer/recorder files               (combine stress fields)
```

## Resolution (2026-05-20)

The spec's premise — *"LAMMPS converges on the same input, find the
missing Trotter ingredient"* — was **false**. Two findings reframed it:

1. **The benchmark was physically ill-posed.** At T=100 K the LJ-argon
   liquid is unstable at P=0 (below liquid–vapour coexistence) — it
   evaporates, and a LAMMPS run of the original input shows LAMMPS
   itself running away. A P(V) probe found the V≈1.1–1.2×10⁴ Å³ range
   is also anomalous (dP/dV>0, mechanically unstable). The benchmark
   was moved to FCC a=5.3 Å / P=2731 bar — a stiff, dP/dV<0 state where
   NPT is well conditioned; `in.npt` and `lammps_npt_avgs.txt`
   regenerated (`lammps.log` is gitignored).

2. **The load-bearing bug was driver-level, not Trotter ordering.**
   `DynamicsDriver::run()` called `setup_pass()` (component `build()`)
   on every invocation, so a chunked `run(burn_in); run(sample)…` wiped
   `NoseHooverNPT`'s NHC-chain state and re-randomised `InitVelocity`
   between chunks. `run()` is now continuable (build once; destroy at
   driver teardown; async-analysis drained per call via a new
   `AnalysisExecutor::quiesce()` barrier).

Final parity: Δ(T)=0.08 %, Δ(KE)=0.08 %, Δ(V)=0.04 % — all far inside
the 5 % gate. Phases A and B were implemented and are correct, but were
not the load-bearing fix; the `propagate_x_and_box` drift factor (spec
wrongly "ruled out" — its Taylor expansion was miscomputed) and a stale
`CellList` grid on box change were also fixed.

## Tasks

- [x] **Phase A**: outer-NHC Trotter ordering swapped to match LAMMPS
      (baro–atoms in `initial_integrate`; atoms–baro in `end_of_step`).
- [x] **Phase B**: `frame.stress_kinetic` field added; `VirialStress`
      routed to it; `update_eps_dot` uses fresh `2·KE_global_`;
      `Barostat` / MPI-reduce / `ZeroStress` / multirank test updated.
- [x] **Drift factor**: `propagate_x_and_box` now uses the exact MTK
      factor `(e^{ε̇dt}−1)/(ε̇dt)` (`expm1c`) instead of bare `sinhc`.
- [x] **CellList grid**: `CellList::build()` re-runs `compute_grid()` /
      `build_stencil()` when the box edge changed (NPT box scaling).
- [x] **DynamicsDriver continuable**: `run()` builds components once,
      keeps a cumulative step counter, defers `destroy()` to the
      destructor, and drains async analysis per call.
- [x] **Benchmark regenerated**: `in.npt` / `lammps_npt_avgs.txt` /
      `lammps.log` for the stable a=5.3 Å, P=2731 bar state point.
- [x] Parity test enabled (`DISABLED_` → `Slow_`), 5 % gates, docstring
      rewritten with the resolution history.
- [x] Full ctest suite green — 93/93, including the enabled parity test
      and the `NoseHooverNptRelaxesPressure` smoke test.
- [x] `msu-pe-validation.md` — no explicit NPT-parity blocker row
      exists; its NPT-dependent items (C.2/C.3, Phase D) are implicitly
      unblocked now that `NoseHooverNPT` is validated.
- [ ] **ac-006** (scientific): `H_NPT_extended` conservation diagnostic
      — manual attestation, deferred to the `npt-conservation-gate`
      follow-up spec.

## Testing strategy

**Happy path:**

* `test_cpu_md_lammps_npt_parity` (after rename to `Slow_`) passes with
  ⟨T⟩, ⟨V⟩, ⟨KE⟩ within 5 % of the LAMMPS reference for the canonical
  256-atom LJ-argon, T=100 K, P=0, τ_P=1000 fs, dt=1 fs, 50 000-step run.
* Full ctest (91 tests) stays green.

**Edge cases:**

* `NoseHooverNptRelaxesPressure` (compressed-init smoke test) continues
  to pass — Phase A operator-swap is structural and should not regress
  the smoke test, but verify.
* `test_md_set_api`, `test_md_plan_equivalence` (which add
  `NoseHooverNPT` without running it to convergence) continue to pass.
* The MPI integration test `test_npt_multirank` (if it exists) continues
  to give bit-identical 1-rank vs N-rank ⟨V⟩ and ⟨T⟩ trajectories.

**Domain validation** (`$META.science.required: true`):

* Phase C diagnostic — `H_NPT_extended` conserved to ≤ 1 × 10⁻³ eV per
  10 000 steps on the same 256-atom LJ-argon system. This is the
  canonical MTK correctness check; surface as a regression gate after
  the parity fix lands.

**Performance:**

* No kernel-timing change expected — operator swap is pure ordering.
  Phase B adds one extra `Frame` field read in observers; should be
  unmeasurable at the per-step granularity.

## Out of scope

* **Other ensembles.** `Barostat` (Berendsen) is internally consistent
  and passes its existing smoke tests; we don't touch it unless Phase B
  is needed (where it's a mechanical field-name change).
* **Anisotropic / triclinic NPT.** `NoseHooverNPT` is isotropic-only; we
  keep that scope.
* **CUDA `cuda::md::NoseHooverNPT`.** Not implemented yet; the CPU fix
  here will be ported when the CUDA migration spec
  (`cuda-driver-taxonomy`) reaches NPT.
* **MSU paper-scale run.** That's `msu-pe-validation.md` Phase D's job;
  this spec only unblocks it.
* **Compute-Hamiltonian conservation gate as a regular ctest.** Phase C
  produces the diagnostic but the gate itself is a follow-up spec
  (`npt-conservation-gate`) — not bundled here.

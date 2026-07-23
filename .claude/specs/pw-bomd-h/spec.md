# Spec: `pw-bomd-h` — plane-wave BOMD foothold (H₂ NVE)

**Slug:** `pw-bomd-h`  
**Domain:** `pw` / dynamics (plane-wave Born–Oppenheimer MD)  
**Queue:** CP2K capability loop — incomplete `[~]`  
**Branch:** `wf/pw-bomd-h`  
**Status:** open  
**Resume:** queue `[~]`; no prior open spec / unmerged `wf/pw-bomd-h` / regression stub

## Problem

Roadmap marks **plane-wave BOMD** as product-ready
(`cpu::pw::PWEvaluator` + `with_gradient(true)` bound at `md::Phase::force`),
and analytic ionic forces ship (`pw_forces_analytic`). There is **no**
end-to-end regression under `regressions/pw/` that proves the **MD product
path** (NVE + energy/force kind `PlaneWave`) for a short trajectory, nor an
offline-oracle pin for **energy drift** under this CP2K capability loop.

Queue row: `H / Si | UPF | energy drift | PW analytic forces + MD ✓; short NVE vs CP2K`.

## Goal

Ship a **dep-safe PW BOMD foothold**:

1. Public API only: `cpu::md::MD` + `VelocityVerlet` +
   `cpu::pw::PWEvaluator{PW}.with_gradient(true)` on `md::Phase::force`.
2. Minimal geometry: **two hydrogen atoms** (toy empty `upf_path`, matching
   existing PW force unit/regression footholds) in a cubic cell — fast, no
   external UPF file at ctest.
3. Short **NVE** trajectory; assert finite PE/KE and **bounded total-energy
   drift** (literature symplectic / BOMD conservation envelope).
4. Offline pin recorded in `regressions/pw/pw-bomd-h/REFERENCE.md` (CP2K twin
   when binary available; else literature + product twin with
   `perf_exempt`).

## Out of scope

- GPW BOMD (`gpw-bomd-water` — separate, dep on GPW forces pin)
- XL-BOMD / CPMD (`xl-bomd-water`, `cpmd-water`)
- Si bulk UPF long traj / production AIMD thermostats
- True GGA / spin / metallic smearing PW (other queue slugs)
- Live CP2K from ctest
- MPI multi-rank scaling
- Public `*Config` bags / `make_*` factories
- Changing `PWEvaluator` force math (already analytic)

## Public API (minimal — already shipped; pin only)

```cpp
cpu::Frame frame(/* 2 H */);
// masses, types Z=1, positions [Å], cubic simbox

cpu::pw::PW pw;
pw.ecut_ry = 5.0;
pw.k_mesh = 1;
pw.tol = 1e-4;
pw.max_iter = 12;
// upf_path empty → toy V_loc (same foothold as h_force_*)

cpu::md::MD sim;
sim.set_integrator(cpu::md::VelocityVerlet{/*dt [fs]=*/0.5});
sim.add(cpu::md::InitVelocity{/*T [K]=*/50.0, /*seed=*/7u});
sim.on(md::Phase::force)
    .add(cpu::pw::PWEvaluator{std::move(pw)}.with_gradient(true));
sim.build();
sim.run(frame, /*n_steps=*/10);
// frame.meta kEnergy [eV], kKinetic [eV]; forces written each step [eV/Å]
```

Rules:

- Parameters are **members** of `cpu::pw::PW` (no public `*Config`).
- Force binding is explicit `on(md::Phase::force).add(...)` — no categorical
  SCF setter bag.
- Frame is sole I/O currency.

## Tasks

1. **Spec** (this file) + INDEX Live row — **before** production code.
2. **Unit test** `tests/cpu/planewave/test_pw_bomd.cpp` (or equivalent under
   `tests/cpu/planewave/`): MD + `PWEvaluator.with_gradient(true)` runs ≥2
   steps on H₂ toy; PE/KE finite. **No e2e under tests/**.
3. **Offline pin + regression** `regressions/pw/pw-bomd-h/`:
   - `REFERENCE.md` — method, geometry, settings, observable, tol, oracle
     citation / unavailability, perf gate.
   - `reg_pw-bomd-h.cpp` — public headers only; short NVE; assert acceptance.
   - Register in `regressions/CMakeLists.txt` under
     `ATV_BUILD_PW` (and MD product on).
4. Build + run unit test and `reg_pw-bomd-h` until PASS.
5. **Perf gate:** fair wall vs CP2K if binary present; else
   `perf_exempt=true` with reason in REFERENCE.
6. Commit on `wf/pw-bomd-h` only (never `dev`; never push).

## Acceptance criteria

| # | Criterion | How verified |
| - | --------- | ------------ |
| A1 | Public path is `MD` + `VelocityVerlet` + `PWEvaluator.with_gradient(true)` at `md::Phase::force`; no public `*Config` / `make_*` | Code review + compile |
| A2 | Unit test: ≥2 MD steps on H₂ toy PW; PE and KE finite | `test_pw_bomd` (or named equivalent) exit 0 |
| A3 | Regression `reg_pw-bomd-h` builds with `ATV_BUILD_PW` + `ATV_BUILD_REGRESSIONS` (+ MD on) | CMake + link |
| A4 | After short NVE: PE finite, KE finite and ≥ 0 | Regression exit 0 |
| A5 | Total energy drift bounded: \(\max_t |E(t)-E_0| / \max(1, |E_0|) \le\) stated tol in REFERENCE (default **0.25** for toy SCF BOMD foothold) | Regression |
| A6 | No trajectory explosion: positions remain finite; optional soft bond-length ceiling if H₂ | Regression |
| A7 | `REFERENCE.md` cites offline CP2K (version + input summary) **or** explicit binary unavailability + literature basis for NVE/BOMD energy conservation pin | File present |
| A8 | Dual gate: correctness (A4–A6) **and** (`within_10x` **or** `perf_exempt=true` with reason) | REFERENCE + report |
| A9 | No live third-party process from ctest / regression binary | Grep / design |
| A10 | Runtime of regression ≤ 60 s on login-node class CPU | Manual wall check |

## Offline pin strategy

1. Prefer **CP2K** plane-wave / Quickstep BOMD twin (same H₂ cell, comparable
   cutoff, short NVE) run **offline** via `tools/cp2k_oracle/`; record version,
   input digest, and energy-drift envelope in REFERENCE.
2. If no CP2K (or QE) binary on the worker host:
   - **Observable pin:** NVE total-energy conservation envelope for a
     symplectic integrator when \(\mathbf{F} \approx -\nabla E\)
     (Frenkel & Smit §4.3; Allen & Tildesley §3.2; Marx & Hutter, *Ab Initio
     Molecular Dynamics*, Born–Oppenheimer MD chapter — energy conservation
     as the BOMD diagnostic).
   - Tolerances are **envelopes** (catch broken force binding / NaN / secular
     leak), not golden snapshots of Atomiverse digits.
   - `perf_exempt=true`, reason: no offline CP2K twin binary on host.
   - Do **not** invent absolute total-energy digits from memory.

## Files (expected touch list)

- `.claude/specs/pw-bomd-h/spec.md` (this)
- `.claude/specs/INDEX.md` (Live row)
- `tests/cpu/planewave/test_pw_bomd.cpp` (unit)
- `tests/CMakeLists.txt` or planewave test registration (if needed)
- `regressions/pw/pw-bomd-h/reg_pw-bomd-h.cpp`
- `regressions/pw/pw-bomd-h/REFERENCE.md`
- `regressions/CMakeLists.txt`

## Done when

All acceptance rows A1–A10 hold; `acceptance_met=true`; branch committed;
report JSON with `status=green` (or `blocked` only for speed after one
optimize pass / missing dep that cannot be fixed in-slug).

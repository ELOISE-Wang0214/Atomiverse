# Spec: `gpw-bomd-water` — GPW BOMD short trajectory (water NVE)

**Slug:** `gpw-bomd-water`  
**Domain:** `gpw` / dynamics (Gaussian-and-plane-waves Born–Oppenheimer MD)  
**Queue:** CP2K capability loop — incomplete `[~]`  
**Branch:** `wf/gpw-bomd-water`  
**Status:** open  
**Resume:** queue `[~]`; dep `gpw-forces-water` ✓; no prior open spec / unmerged
`wf/gpw-bomd-water` / regression stub on assign wave 2

## Problem

Wave-1 shipped **GPW analytic forces** (`gpw-forces-water`) and the MD product
already accepts `cpu::GPW` as an energy/force kind at `md::Phase::force`.
There is **no** end-to-end regression under `regressions/gpw/` that proves the
**MD product path** (NVE + `GPWEvaluator` with `compute_forces = true`) for a
short water trajectory, nor an offline-oracle pin for **energy drift**.

Queue row: `water | DZVP | energy drift | GPWEvaluator forces + MD ✓; need
CP2K MD reference trajectory pin`. Enables later `xl-bomd-water`.

## Goal

Ship a **dep-safe GPW BOMD foothold**:

1. Public API only: `cpu::md::MD` + `VelocityVerlet` +
   `cpu::gpw::GPWEvaluator{…, compute_forces = true}` on `md::Phase::force`.
2. Geometry: **isolated H₂O** (same CP2K getting-started coordinates as
   `gpw-forces-water` / `gpw-water-pbe-mt`), cubic cell, MT / nonperiodic.
3. Basis/PP product twin: **DZVP-MOLOPT-GTH** + GTH-PBE; XC **PBE** (or
   cheaper product-safe settings documented in REFERENCE if needed for ≤60 s).
4. Short **NVE** trajectory; assert finite PE/KE and **bounded total-energy
   drift** (literature symplectic / BOMD conservation envelope, and/or offline
   CP2K AIMD pin when binary available).
5. Offline pin in `regressions/gpw/gpw-bomd-water/REFERENCE.md`.

## Out of scope

- XL-BOMD / CPMD / ASPC product (`xl-bomd-water`, `cpmd-water`)
- GPW cell stress (`gpw-stress-water`)
- Thermostatted production AIMD / long traj
- Live CP2K from ctest
- MPI multi-rank scaling
- Public `*Config` bags / `make_*` factories
- Changing analytic force math (already gated by `gpw-forces-water`)

## Public API (minimal — composition of shipped pieces)

```cpp
cpu::Frame frame(/* 3 atoms H2O */);
// masses, types Z, positions [Å], cubic simbox, mild velocities

const cpu::gpw::Kind O{molrs::Element::O,
    cpu::gpw::GaussianBasis(molrs::Element::O, "DZVP-MOLOPT-GTH"),
    cpu::gpw::GthPotential(molrs::Element::O, "GTH-PBE-q6")};
const cpu::gpw::Kind H{molrs::Element::H,
    cpu::gpw::GaussianBasis(molrs::Element::H, "DZVP-MOLOPT-GTH"),
    cpu::gpw::GthPotential(molrs::Element::H, "GTH-PBE-q1")};

cpu::gpw::GPWEvaluator eval(
    /*kinds*/ {O, H}, /*cutoff_ry*/ /* product-safe */,
    /*max_iter*/ …, /*tol*/ …, /*mix*/ …,
    /*xc*/ cpu::gpw::XcKind::PBE, /*nonperiodic*/ true,
    /*compute_forces*/ true);

cpu::md::MD sim;
sim.set_integrator(cpu::md::VelocityVerlet{/*dt [fs]=*/ …});
sim.on(md::Phase::force).add(std::move(eval));
sim.build();
sim.run(frame, /*n_steps*/ …);
// frame.meta kEnergy [eV], kKinetic [eV]; forces each step [eV/Å]
```

Rules:

- Parameters are **constructor args / members** of `GPWEvaluator` (no public
  `*Config`).
- Force binding is explicit `on(md::Phase::force).add(...)`.
- Frame is sole I/O currency.
- Fixed cell (NVE): `GPWEvaluator` already rejects box changes after `build()`.

## Tasks

1. **Spec** (this file) + INDEX Live row — **before** production code.
2. **Warm-start (product BOMD):** `GPWEvaluator::launch` passes previous
   converged MO coefficients into `GpwDft::run(..., guess_C)` when available
   (solver already accepts `guess_C`). Unit-gate that multi-step MD converges
   and stays finite. No public Config bag.
3. **Unit test** `tests/cpu/gpw/test_gpw_bomd.cpp`: MD + `GPWEvaluator`
   (`compute_forces=true`) runs ≥2 steps on water (cheap settings OK); PE/KE
   finite. **No e2e under tests/**.
4. **Offline pin + regression** `regressions/gpw/gpw-bomd-water/`:
   - `REFERENCE.md` — method, geometry, settings, observable, tol, oracle
     citation / unavailability, perf gate.
   - `reg_gpw-bomd-water.cpp` — public headers only; short NVE; assert
     acceptance.
   - Register in `regressions/CMakeLists.txt` under `ATV_BUILD_GPW` +
     `ATV_BUILD_MD`.
5. Build + run unit test and `reg_gpw-bomd-water` until PASS (≤ 60 s wall).
6. **Perf gate:** fair wall vs CP2K if binary present; else
   `perf_exempt=true` with reason in REFERENCE.
7. Commit on `wf/gpw-bomd-water` only (never `dev`; never push).

## Acceptance criteria

| # | Criterion | How verified |
| - | --------- | ------------ |
| A1 | Public path is `MD` + `VelocityVerlet` + `GPWEvaluator` with `compute_forces=true` at `md::Phase::force`; no public `*Config` / `make_*` | Code review + compile |
| A2 | Unit test: ≥2 MD steps on water GPW; PE and KE finite | `test_gpw_bomd` exit 0 |
| A3 | Regression `reg_gpw-bomd-water` builds with `ATV_BUILD_GPW` + `ATV_BUILD_MD` + `ATV_BUILD_REGRESSIONS` | CMake + link |
| A4 | After short NVE: PE finite, KE finite and ≥ 0 | Regression exit 0 |
| A5 | Total energy drift bounded: \(\max_t |E(t)-E_0| / \max(1, |E_0|) \le\) stated tol in REFERENCE (default **0.05** relative; may relax to **0.25** only if documented for cheap-cutoff foothold) | Regression |
| A6 | No trajectory explosion: all coordinates remain finite; optional soft O–H bond ceiling | Regression |
| A7 | `REFERENCE.md` cites offline CP2K (version + MD input summary + drift envelope) **or** explicit binary unavailability + literature basis for NVE/BOMD energy conservation pin | File present |
| A8 | Dual gate: correctness (A4–A6) **and** (`within_10x` **or** `perf_exempt=true` with reason) | REFERENCE + report |
| A9 | No live third-party process from ctest / regression binary | Grep / design |
| A10 | Runtime of regression ≤ 60 s on login-node class CPU | Manual wall check |

## Offline pin strategy

1. Prefer **CP2K Quickstep BOMD / MD** twin (same H₂O cell, DZVP-MOLOPT-GTH,
   GTH-PBE, PBE, MT, short NVE) run **offline** via `tools/cp2k_oracle/`;
   record version, input digest, and energy-drift envelope in REFERENCE.
2. If no CP2K binary on the worker host:
   - **Observable pin:** NVE total-energy conservation envelope for a
     symplectic integrator when \(\mathbf{F} \approx -\nabla E\)
     (Frenkel & Smit §4.3; Allen & Tildesley §3.2; Marx & Hutter, *Ab Initio
     Molecular Dynamics*, Born–Oppenheimer MD chapter — energy conservation
     as the BOMD diagnostic). Force quality already pinned by
     `gpw-forces-water` vs CP2K 2026.1.
   - Tolerances are **envelopes** (catch broken force binding / NaN / secular
     leak), not golden snapshots of Atomiverse digits.
   - `perf_exempt=true`, reason: no offline CP2K twin binary on host.
   - Do **not** invent absolute total-energy digits or CP2K traj numbers
     from memory.

## Files (expected touch list)

- `.claude/specs/gpw-bomd-water/spec.md` (this)
- `.claude/specs/INDEX.md` (Live row)
- `include/atomiverse/cpu/gpw/gpw_evaluator.h` / `src/cpu/gpw/gpw_evaluator.cpp`
  (warm-start only if needed for A5/A10)
- `tests/cpu/gpw/test_gpw_bomd.cpp`
- `tests/CMakeLists.txt`
- `regressions/gpw/gpw-bomd-water/reg_gpw-bomd-water.cpp`
- `regressions/gpw/gpw-bomd-water/REFERENCE.md`
- `regressions/CMakeLists.txt`
- optional: `tools/cp2k_oracle/inputs/gpw-bomd-water*.inp`

## Done when

All acceptance rows A1–A10 hold; `acceptance_met=true`; branch committed;
report JSON with `status=green` (or `blocked` only for speed after one
optimize pass / missing dep that cannot be fixed in-slug).

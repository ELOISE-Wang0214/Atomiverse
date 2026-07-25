# Spec: `gpw-variable-cell-ice` — GPW variable-cell relax (ice-like chain)

**Slug:** `gpw-variable-cell-ice`
**Domain:** `gpw` (Gaussian-and-plane-waves DFT)
**Queue:** CP2K capability loop — wave-3 `[~]` resume (write-first)
**Branch:** `wf/gpw-variable-cell-ice`
**Status:** open
**Deps:** `gpw-stress-water` ✓ (strain-FD GPW stress on `kStress` /
`last_stress_ha_bohr3()`); `AnisotropicCellRelax` + `Minimize` preset ✓

## Problem

GPW ships a cell stress (strain-FD, `gpw-stress-water`), but nothing consumes
it: `GPWEvaluator::launch` **throws** when the frame box changes after
`build()` ("GPW does not support a variable cell"), so no ATV driver can relax
a GPW cell. The queue row productizes GPW stress into variable-cell
relaxation with a CP2K `CELL_OPT`-class twin (ice Ih row; see *System* for the
honest minimal system).

## Goal

1. **GPW variable cell:** `GPWEvaluator` re-targets its solver to the frame's
   current orthorhombic box on every `launch()` (GpwDft::run already rebuilds
   all geometry-dependent structure per call), keeping the MO warm start.
2. **Cell relaxation product path:** drive `cpu::minimizer::Minimize` +
   `AnisotropicCellRelax` with GPW stress to relax an orthorhombic cell to
   the zero-stress (or target-stress) cell, per-axis maskable — the
   CP2K `CELL_OPT` twin at fixed fractional coordinates.
3. **Regression** pinned against offline CP2K 2026.1 (never live at ctest):
   relaxed lattice constant + energy + residual stress on a periodic
   H-bonded water chain (ice-like motif; see *System*).

## Out of scope

- Analytic (Nielsen–Martin/Pulay) GPW stress — still strain-FD v1
  (`gpw-analytic-stress` queue row)
- Triclinic GPW cells / shear relaxation (GPW cell is orthorhombic)
- Joint atom+cell GPW relaxation (ions relax); fixed fractional coordinates
  only — same functional E(cell) in both engines
- Full ice Ih (8 H₂O / 24-atom orthorhombic ice XI cell): one GPW SCF at
  product cutoff is already outside the 60 s regression budget; the wall-time
  law wins. Documented in REFERENCE.md as follow-on after analytic stress.
- k-points, UKS stress, dispersion-corrected stress
- Fixing the classical-vs-electronic `kStress` convention divergence
  (see *Discovered debt*)
- CUDA

## Public API (binding shape)

```cpp
// 1) GPW evaluator follows the frame box (was: throw on box change)
cpu::gpw::GPWEvaluator eval({O, H}, /*cutoff_ry*/ 300.0, /*max_iter*/ 200,
                            /*tol*/ 1e-7, /*mix*/ 0.3, cpu::gpw::XcKind::PBE,
                            /*nonperiodic*/ false, /*compute_forces*/ false);
eval.compute_stress = true;          // existing member (strain-FD → kStress)
eval.stress_axes = {true, false, false}; // NEW member: FD only the relaxed axes

// 2) Cell relax: existing Optimizer, two NEW plain members (no Config bag)
cpu::minimizer::AnisotropicCellRelax cell(/*P_target*/ 0.0, /*alpha*/ …);
cell.stress_is_virial = false;  // NEW: kStress already intensive [eV/Å³]
                                // (PW/GPW evaluators); default true = classical
                                // extensive virial [eV] (divide by V)
cell.relax_x = true; cell.relax_y = false; cell.relax_z = false; // NEW mask

// 3) Termination is a Criterion (intensive-stress twin of PressureConverged)
opt.add_criterion(cpu::StressConverged{/*sigma_target*/ 0.0,
                                       /*stol [eV/Å³]*/ …,
                                       /*x,y,z*/ true, false, false});

cpu::minimizer::Minimize opt;
opt.set_optimizer_ref(cell);
opt.on(minimizer::Phase::force).add_ref(eval);
opt.build();
opt.run(frame, budget);   // Frame in / Frame out; relaxed box on frame.simbox
```

Rules: parameters are constructor args / public members (iron law — no
`*Config` bag, no factories); regression consumes public headers only;
Frame stays sole I/O (relaxed cell on `frame.simbox`, energy `kEnergy` [eV],
stress `kStress` [eV/Å³]).

## Design notes

- `GPWEvaluator::launch`: replace the box-change throw with re-targeting
  (`solver_->cell[d] = frame box [Bohr]`). Safe because `GpwDft::run`
  rebuilds grids / lattice images / S/T/Vpp from `cell` per call
  (`make_uniform_grid(cfg.cell, …)` in `run`); AO count is
  geometry-independent so the MO warm start stays valid.
- `gpw_stress(atoms, method, h, axes)` gains an axes mask (skip FD SCFs of
  masked axes → 0 in Voigt) and warm-starts each strained SCF from the
  unstrained converged MOs (`guess_C`) — pure cost cut, values unchanged
  (verified against full-axes result in the unit test).
- `AnisotropicCellRelax` keeps its virial default; the two new members make
  the `kStress` convention explicit at the read site.
- `StressConverged` reads diagonal `kStress` directly (intensive), maskable;
  lives in `cpu/driver/criteria.h` beside `PressureConverged`.

## System (honest minimal "ice")

Loop precedent: `gpw-pbe-d3bj-ice` pinned a periodic H-bonded water pair for
the ice row. Full ice Ih needs the 24-atom orthorhombic ice XI cell — a
single GPW SCF there blows the 60 s budget before the first relax step.

Twin system: **periodic proton-ordered water chain** — 1 H₂O per cell,
H-bonded to its own x-images (the ice Ih c-axis strand motif), vacuum-padded
in y/z. Variable-cell relax of the chain axis `a` at fixed fractional
coordinates and fixed b = c = 8 Å (per-axis mask — the `CELL_OPT`
`CONSTRAINT` twin).

Fractional coordinates (fixed throughout; b = c = 8 Å; anchored so the
monomer is gas-phase water at a = 3.0 Å):

| Atom | f_x | f_y | f_z |
| ---- | --- | --- | --- |
| O  | 0.000000 | 0 | 0.000000 |
| H1 | 0.319067 | 0 | 0.000000 |
| H2 | −0.079900 | 0 | 0.115830 |

**Basis/XC pivot (measured, binding):** the queue row's DZVP-MOLOPT-GTH is
not usable for a *periodic* ≤ 60 s regression — ATV's periodic MOLOPT image
sums cost ~10 min per 3-SCF relax launch on this chain (same wall
`gpw-pbesol-si` documented when it retreated to MT), vs **1.9 s** with the
compact GTH bases. The variable-cell machinery is basis-agnostic; the twin
uses the CP2K standard GTH basis + LDA so the periodic relax fits the law:

Method twin (both engines): **PADE LDA**, O **SZV-GTH-PADE** +
H **DZV-GTH-PADE** (exponents verified identical to CP2K `BASIS_SET` file),
**GTH-PADE-q6/q1**, CUTOFF 300 Ry (CP2K REL_CUTOFF 60, NGRIDS 4),
EPS_SCF 1e-7, periodic XYZ. MOLOPT re-twin is follow-on work for the GPW
perf campaign (periodic collocation image sums).

## Oracle (offline; pins hard-coded)

CP2K 2026.1 (`cp2k.psmp`, EasyBuild `2026.1-gmpflf-2026.03`, Arrhenius):

1. **E(a) scan** at fixed fractional coordinates (SCALED coords) →
   equilibrium a\* from the analytic minimum of a cubic fit; this is the
   zero of σ_xx by construction.
2. **`STRESS_TENSOR NUMERICAL`** at bracketing cells → σ_xx sign/magnitude
   pins (same-method FD twin as `gpw-stress-water`).
3. **`CELL_OPT`** (all atoms fixed, CONSTRAINT YZ) — the **headline pin**
   (verified: fixed atoms keep scaled coordinates, so it minimizes the same
   E(a); result 3.211148 Å agrees with the scan fit 3.2124 and the numerical
   stress zero 3.219 within 0.008 Å).

Inputs under `tools/cp2k_oracle/inputs/gpw-variable-cell-ice/`; numbers +
walls in `regressions/gpw/gpw-variable-cell-ice/REFERENCE.md`.

## Discovered debt (named per iron law)

- `kStress` convention divergence (pre-existing): classical pair/bond writers
  accumulate the **extensive virial [eV]**; PW (`with_stress`) and GPW
  (`compute_stress`) write **intensive stress [eV/Å³]** on the same meta key.
  `IsotropicCellRelax` / `AnisotropicCellRelax` / `PressureConverged` divide
  by V (virial convention) — with P_target = 0 the fixed point is unchanged,
  which is why `reg_h_cell_isotropic` (PW) "works". Not locally fixable
  (touches every classical writer + NPT/MC readers; classical potentials are
  carrying unrelated in-flight metalwalls work). This slug makes the
  convention **explicit at the read site** (`stress_is_virial`,
  `StressConverged` doc) and reports the divergence for a dedicated
  follow-up slug.
- `copy_gpw` in `internal/gpw_stress.cpp` drops `d4`, `plus_u` and
  `atom_spin_guess` — strained SCFs would silently lose the D4 / DFT+U /
  spin-seed terms when those are enabled. Local, fixed in this slug.

## Deliverables

| Path | What |
| ---- | ---- |
| `include/atomiverse/cpu/gpw/gpw_evaluator.h` + `src/cpu/gpw/gpw_evaluator.cpp` | variable-cell launch; `stress_axes` member |
| `include/atomiverse/cpu/gpw/gpw_solver.h` + `src/cpu/gpw/internal/gpw_stress.cpp` | `gpw_stress` axes mask + warm start; `copy_gpw` field fix |
| `include/atomiverse/cpu/minimize/anisotropic_cell_relax.h` + `src/cpu/minimize/anisotropic_cell_relax.cpp` | `stress_is_virial`, `relax_x/y/z` |
| `include/atomiverse/cpu/driver/criteria.h` | `StressConverged` |
| `tests/cpu/minimize/test_anisotropic_cell.cpp` | mask + intensive-mode + criterion units |
| `tests/cpu/gpw/test_cpu_gpw_stress.cpp` | axes-mask units (toy periodic fixture, seconds) |
| `tests/cpu/gpw/test_gpw_method.cpp` | GpwDft variable-cell run-rebuild unit |
| `regressions/gpw/gpw-variable-cell-ice/` | `reg_gpw-variable-cell-ice.cpp` + `REFERENCE.md` |
| `tools/cp2k_oracle/inputs/gpw-variable-cell-ice/` | offline twin inputs |
| `regressions/CMakeLists.txt` | `atv_add_regression(gpw gpw-variable-cell-ice)` |

---
slug: gpw-d3-atm-forces
status: approved
created: 2026-07-31
revised: 2026-07-31
domain: gpw / dispersion forces
queue: cp2k-method-loop
deferred: false
depends_on: []
unblocks: [gpw-d3-atm-forces]
grilled: false
sibling: gpw-d4-atm-forces (serialize on force-assembly ATM fold)
---

# Spec: `gpw-d3-atm-forces` — GPW analytic ionic forces with DFT-D3 three-body ATM (C9)

**Slug:** `gpw-d3-atm-forces`
**Domain:** GPW dispersion — DFT-D3 (zero damping) three-body Axilrod–Teller–Muto
**Queue row:** wave-14 `gpw-d3-atm-forces`
**Oracle:** CP2K 2026.1 offline (`TYPE DFTD3`, `CALCULATE_C9_TERM`), never live in ctest

## Problem

`gpw-d3-atm` shipped the **energy** foothold: `GpwDft::d3_atm = true` folds the
CP2K `CALCULATE_C9_TERM` three-body C9 energy into `e_disp`, pinned against
offline CP2K (−0.00081414118193 Ha, water dimer, `PERIODIC NONE`). The force
machinery exists — `accumulate_atm_c9` in `src/cpu/gpw/internal/dftd3.cpp`
adds explicit-geometry ATM forces plus the CN-chain into `dEdCN`, and the force
assembly (`gpw_forces.cpp`) passes `cfg.d3_atm` through to `dftd3_pbe` — but
the **ATM force path has no dedicated oracle pin**: no regression asserts the
ATM three-body gradient against a number CP2K printed, and the only FD gate
(`GpwDftD3Atm.AnalyticForcesMatchFiniteDifference`) exercises the isolated
cluster branch, never the periodic image enumeration (CASE 2/3, `fac0 = 1/3`).

## Deliverable

A dedicated **CP2K forces pin + FD dual gate** for the D3 ATM gradient, on the
same water dimer as the energy foothold, in both boundary conditions:

1. **Offline CP2K pins** (`tools/cp2k_oracle/inputs/gpw-d3-atm-forces/`):
   four `RUN_TYPE ENERGY_FORCE` twins on the `gpw-d3-atm` geometry —
   - isolated (`PERIODIC NONE`, 12 Å): `CALCULATE_C9_TERM .TRUE.` /
     `.FALSE.`;
   - small periodic cell (`PERIODIC XYZ`, 6 Å, `R_CUTOFF 10.0`): same pair.
   D3 is geometry-only (it never enters the KS potential), so per boundary
   condition the difference `F(C9) − F(two-body)` cancels the electronic force
   bitwise-to-SCF-noise and isolates **exactly the ATM three-body gradient** —
   the quantity this capability adds. Record E_disp and the `FORCES|` tables in
   `REFERENCE.md` (binary + version + input + output line).
2. **Regression** `regressions/gpw/gpw-d3-atm-forces/reg_gpw-d3-atm-forces.cpp`
   (out-of-library consumer, ≤ 60 s), gates:
   - `E_disp(iso, ATM)` and `E_disp(per, ATM)` from standalone `dftd3_pbe`
     vs the CP2K pins;
   - ATM-only force `dftd3_pbe(c9=on).forces − dftd3_pbe(c9=off).forces`
     vs the CP2K `F(C9) − F(two-body)` pin, per boundary condition;
   - **product wiring**: full GPW SCF via public `GPWEvaluator`
     (`d3 = true`, `d3_atm = true`, `compute_forces = true`) minus the
     `d3_atm = false` force assembly equals the standalone ATM difference
     (only SCF residual enters);
   - total product forces vs CP2K's ATM-run `FORCES|` table
     (grid-noise-bounded tolerance, as `gpw-d3-forces`);
   - translation invariance `|Σ_A F_A|` on the isolated system.
3. **Unit FD gate (periodic)**: extend `tests/cpu/gpw/test_cpu_gpw_dftd3.cpp`
   with an analytic-vs-central-FD ATM force test on a **periodic** cell, so the
   image-class enumeration (CASE 2/3, `fac0 = 1/3`) is FD-gated, not only the
   home-cell branch. (The isolated FD test already exists.)

No solver code change is expected; if any gate exposes a defect in
`accumulate_atm_c9` or the force fold, fix it in the same change (iron law:
no silent debt).

## Tolerances (set from measurement, justified in REFERENCE.md)

| Gate | Tolerance |
| --- | --- |
| E_disp vs CP2K (iso + per, ATM on) | 1e-9 Ha |
| ATM-only force vs CP2K difference | ≥ 20× measured floor, ≤ signal/20 |
| Product wiring | 1e-7 Ha/Bohr |
| Total F vs CP2K | 1e-2 Ha/Bohr (oracle grid noise bound, as `gpw-d3-forces`) |
| `\|Σ F\|` isolated | 1e-3 Ha/Bohr |
| Unit FD (periodic ATM) | 1e-7 Ha/Bohr |
| Wall | ≤ 60 s |

## Acceptance (met 2026-07-31)

- [x] Four offline CP2K 2026.1 runs recorded (inputs in-tree, outputs quoted in
      REFERENCE.md with version + command).
- [x] `reg_gpw-d3-atm-forces` green: iso 5.0e-10 / per 9.8e-10 Ha/Bohr vs pins
      (tol 1e-8), wiring 2.2e-17, wall 34.5 s.
- [x] New periodic ATM FD unit test green (2.2e-13 Ha/Bohr).
- [x] All 12 `test_cpu_gpw_dftd3` tests green.
- [x] Perf: T_atv/T_cp2k = 34.5/10.5 ≈ 3.3 — within the 10× fair wall
      (no exemption needed).

**Defects found and fixed by the pin** (`accumulate_atm_c9`,
`src/cpu/gpw/internal/dftd3.cpp`): periodic CASE 3 image-triple double-count
(AB pairs must replicate CP2K's half neighbor list) and missing CP2K `pbc()`
anchoring of the hard-gated truncated triple sum (per-atom ATM forces of the
two molecules came out swapped). Details: regression REFERENCE.md.

## Non-goals

- D4 ATM forces (`gpw-d4-atm-forces`, sibling — serialized).
- ATM stress / virial (separate queue rows).
- k-point force assembly changes (`gpw-dftu-kpoint-forces` owns that region).

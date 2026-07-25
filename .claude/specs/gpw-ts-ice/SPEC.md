# Spec: gpw-ts-ice — Tkatchenko–Scheffler dispersion (GPW)

**Slug:** `gpw-ts-ice`  
**Branch:** `wf/gpw-ts-ice`  
**Domain:** GPW / periodic DFT  
**Oracle:** offline libMBD (CP2K stock has no TYPE TS)  
**Status:** implement

## Goal

Ship free-atom **Tkatchenko–Scheffler (TS)** pairwise dispersion on GPW:

- Public switch `GpwDft::ts` / `GPWEvaluator::ts` (mutually exclusive with
  `d3` / `d3bj` / `d4`).
- `GpwResult::e_disp` [Hartree] + fold into `e_total`.
- In-tree free-atom tables (`basis_data/ts_tables.h`, libMBD `ts_vdw_params`).
- PBE damping \(s_R=0.94\), \(d=20\); pure-RS lattice sum, default rcut 50 Å.
- Optional volume-ratio vector (empty ⇒ free-atom); Hirshfeld from density is
  **follow-on** (not this slug). MBD → `gpw-mbd-ice`.

## Acceptance

1. Unit `test_cpu_gpw_ts`: isolated E_disp matches libMBD to 1e-12 Ha; periodic
   pure-RS pin; analytic forces vs FD; `GpwDft.ts` reports e_disp.
2. Regression `reg_gpw-ts-ice`: periodic water dimer E_disp within **1e-8 Ha** of
   −6.099630773386013e-04 Ha (REFERENCE.md).
3. No public `*Config` / `make_*` factories; Frame I/O unchanged.
4. **perf_exempt**: geometry-only hard gate (same policy as D3BJ/D4).

## Non-goals

- Hirshfeld density-dependent volumes, MBD@rsSCS, vdW-DF.
- CP2K live in ctest.
- CUDA.

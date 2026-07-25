# Spec — `gpw-rvv10-ice`

## Goal

Ship **nonlocal rVV10 correlation** on Γ-point GPW (post-SCF energy on the
final total density), offline-pinned against CP2K `NON_LOCAL TYPE RVV10` and
the community kernel table.

## Scope

- Public flags: `GpwDft::rvv10` / `GPWEvaluator::rvv10` (mutually exclusive with
  `d3` / `d3bj` / `d4` / `ts` / `mbd` / `vdw_df`).
- Internal: `cpu/gpw/internal/rvv10.{h,cpp}` + embedded
  `rvv10_kernel_data` (CP2K `rVV10_kernel_table.dat`, float32 zlib).
- Defaults: b = 6.3, C = 0.0093 (Sabatini / CP2K).
- Product path: **energy-only**, post-SCF (no potential back-propagation).
- Regression: `regressions/gpw/gpw-rvv10-ice` water dimer 12 Å cell, CUTOFF
  200 Ry, DZVP-MOLOPT-SR-GTH + GTH-PBE, PBE + rVV10.
- Unit: `tests/cpu/gpw/test_cpu_gpw_rvv10.cpp` (q₀ + synthetic energy finite /
  β > 0).

## Out of scope

- Self-consistent rVV10 potential / forces / stress.
- vdW-DF2 (LMKLL) — separate slug `gpw-vdw-df2-ice`.
- Periodic MBD / XDM.
- CUDA.

## Acceptance

1. Unit tests green (`test_cpu_gpw_rvv10`).
2. Regression green vs offline pin in `REFERENCE.md` (E_disp Ha).
3. Public energy still eV on Frame when evaluated through evaluator; internal
   Ha on `GpwResult::e_disp` (same as other disp methods).
4. Perf: `perf_exempt` (E_disp pin parity; same policy as `gpw-vdw-df-ice`).

## Oracle

- CP2K 2026.1 (`install` on Arrhenius tree), `TYPE RVV10`, same geometry.
- Input notes: `tools/cp2k_oracle/inputs/gpw-rvv10-ice/`.
- Kernel provenance: CP2K `data/rVV10_kernel_table.dat`.

## Status

acceptance_met = true (unit + regression green; CP2K 2026.1 offline twin;
perf_exempt).

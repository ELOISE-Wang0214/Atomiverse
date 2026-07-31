# gpw-d4-atm-stress — analytic fixed-P cell stress with DFT-D4 ATM folded

Queue row: `docs/developer/cp2k-method-queue.md` (GapCensus wave 2, shift
2026-07-31). Hold cleared when `gpw-d4-atm-forces` merged (closed-form ATM
derivatives already in `cpu::scf::dftd4`).

## Problem

The analytic fixed-P stress path (`gpw_stress` in
`src/cpu/gpw/internal/gpw_stress.cpp`) already folds the ATM three-body term:
the strained-cell dispersion rebuild uses `s9 = cfg.d4_atm ? 1.0 : 0.0`
(Γ line 528, multi-k line 870). The energy path folds too (`gpw-d4-atm`).
**No dedicated pin gates the ATM stress contribution** — the existing
`gpw-d4-stress` regression runs with ATM off, and its analytic-vs-FD dual gate
(residual ~7×10⁻⁷ Ha/Bohr³) cannot resolve the ATM stress increment
(~10⁻¹⁰ Ha/Bohr³ class for the 12 Å water dimer).

## Oracle

**dftd4 3.7.0** strained-cell energies (same no-CP2K-DFTD4 policy as
`gpw-d4-atm-forces` / `gpw-d4-forces`: Arrhenius CP2K 2026.1 is built without
DFTD4, and the hard pin is the same libdftd4 backend CP2K would call).

Pins = printed dispersion energies at the ±strained cells, produced offline:
for each axis a ∈ {x,y,z}, scale cell row a **and** the Cartesian
a-components by (1±h), h = 2×10⁻³ (exactly the `scale_axis` fixed-fractional
convention of `gpw_stress`), run

- `dftd4 -f pbe --json ... POSCAR_strained` (s9 = 1, ATM on)
- `dftd4 -f pbe --mbdscale 0.0 --json ... POSCAR_strained` (two-body)

12 strained pins + the 2 base pins already in
`tools/cp2k_oracle/inputs/gpw-d4-atm/`. Geometry: the `gpw-d4-atm-forces`
periodic water dimer (12 Å cubic, O O H H H H POSCAR order).

The oracle's periodic ATM **gradient** non-conservativity (~2.4×10⁻⁸,
documented in `gpw-d4-atm-forces/REFERENCE.md`) is irrelevant here by
construction: every gate is on **printed energies**, never on its gradient.

## Deliverable

`regressions/gpw/gpw-d4-atm-stress/` (REFERENCE.md + reg binary), public API
only. Gates:

1. **Strained-energy parity (hard pin):** ATV
   `cpu::scf::dftd4_energy(..., s9=1)` at each of the 6 strained geometries
   matches the dftd4 3.7.0 printed energy to ≤ 1×10⁻⁹ Ha (expected ~6×10⁻¹¹
   Bohr-conversion constant, cancels in differences).
2. **Dispersion stress parity (decisive, dedicated pin):** for each axis,
   σ_aa^disp = −(1/Ω)·[E(+h) − E(−h)]/(2h) from ATV s9=1 energies vs the same
   FD over the oracle's printed energies — tight (≤ 1×10⁻¹² Ha/Bohr³ class,
   constant offset cancels). Same for the **ATM-only** contribution
   σ_aa^ATM = σ_aa^disp(s9=1) − σ_aa^disp(s9=0), plus a non-vacuity floor
   (max|σ^ATM| above noise) so the gate fails if the ATM term is dropped from
   the strained rebuild.
3. **Product dual gate:** `GPWEvaluator` with `d4=true, d4_atm=true,
   compute_stress=true` (analytic fixed-P) vs the multi-SCF strain-FD twin
   (`gpw_stress_fd`) on the same d4+ATM method — ≤ ~1×10⁻⁶ class residual
   (sibling `gpw-d4-stress` measured 7.4×10⁻⁷; tol 2×10⁻⁵ hard). Product
   `e_disp` must equal the ATM base pin (≤ 1×10⁻⁸ Ha) — fails if `d4_atm`
   stops reaching the method.
4. **Product ATM fold wiring:** analytic `gpw_stress` with `d4_atm=true`
   minus with `d4_atm=false` (same warm start) equals the standalone
   σ^ATM of gate 2 — proves `cfg.d4_atm` reaches the strained dispersion
   rebuild inside the stress path (tolerance set from measured noise;
   dispersion is additive so the electronic Lagrangian cancels).

Unit: one strained-cell s9=1 energy pin in
`tests/cpu/gpw/test_cpu_gpw_dftd4.cpp` (single function call, fast).

## Budget / perf

Product analytic ~21 s + FD twin (x-axis only, warm-started) + warm-started
2b analytic — target ≤ ~60 s wall. `perf_exempt=true`, same reason as
siblings: Arrhenius CP2K lacks DFTD4, no fair QS stress wall; the hard gate
is dftd4 printed-energy parity.

## Acceptance

- [ ] 12 strained-cell pins under `tools/cp2k_oracle/inputs/gpw-d4-atm/`
      (re-runnable: binary + version + input + output line in README).
- [ ] `regressions/gpw/gpw-d4-atm-stress/` green with gates 1–4.
- [ ] Unit pin added; scoped tests green.
- [ ] REFERENCE.md cites dftd4 3.7.0 and records first-green numbers.

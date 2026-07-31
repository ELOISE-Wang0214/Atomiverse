---
slug: gpw-mbd-kpoint-stress
status: approved
created: 2026-07-31
revised: 2026-07-31
domain: gpw / cell stress / dispersion
queue: cp2k-method-loop
deferred: false
depends_on: []
unblocks: [gpw-mbd-kpoint-stress]
---

# Spec: `gpw-mbd-kpoint-stress` — multi-k GPW analytic stress folds periodic MBD@rsSCS

**Slug:** `gpw-mbd-kpoint-stress`
**Domain:** GPW cell stress (CPU), multi-k fixed-\(P_k\) Lagrangian
**Queue row:** `gpw-mbd-kpoint-stress` (GapCensus wave 13, assigned wave 1
2026-07-31)

## Summary

The Γ analytic fixed-\(P\) stress Lagrangian (`energy_fixed_P`) rebuilds
\(E_{\mathrm{disp}}^{\mathrm{MBD}}\) at every strained cell
(`regressions/gpw/gpw-mbd-stress`, shipped wave 14). The **multi-k** k-summed
Lagrangian (`energy_fixed_P_k`) never folded MBD; since silently dropping
\(E_{\mathrm{disp}}\) from the strain derivative is a wrong tensor, the guard
landed with `reg_gpw-mbd-stress` makes `gpw_stress` **throw** on
`mbd && k_mesh != Γ` (`src/cpu/gpw/internal/gpw_stress.cpp:351`). This spec
delivers the declared follow-on: fold MBD@rsSCS into the multi-k fixed-\(P_k\)
stress Lagrangian, identically to the Γ path.

## Physics

MBD@rsSCS with free-atom (or fixed volume-ratio) parameters is
**geometry-only**: \(E_{\mathrm{disp}}(\{Z\},\{\mathbf R\},\mathbf h)\)
does not depend on the SCF density or on the electronic k mesh. Under the
fixed-\(P_k\) Hellmann–Feynman cell theorem the dispersion contribution to
\(\partial E_L/\partial\varepsilon_{aa}\) is therefore obtained exactly as at
Γ: rebuild \(E_{\mathrm{disp}}\) at each strained cell (positions at fixed
fractional coordinates, cell axis scaled, MBD Ewald + its own
`mbd_k_mesh` lattice sum re-evaluated) inside the FD of the
non-self-consistent Lagrangian. Nothing else in the fold is k-dependent.

## Change surface

`src/cpu/gpw/internal/gpw_stress.cpp` only (multi-k path; Γ pins untouched —
coordination note for the two Γ D3 stress workers in the same file):

1. `energy_fixed_P_k` mutual-exclusion count includes `cfg.mbd`
   (d3 / d3bj / d4 / ts / xdm / mbd mutually exclusive, same as Γ).
2. Delete the `cfg.mbd` refuse block (the line-351 throw).
3. `energy_fixed_P_k` dispersion block gains the MBD branch, verbatim
   semantics of the Γ branch: `mbd_k_mesh` all-zero ⇒ isolated free-atom
   MBD; non-zero ⇒ Ewald + MBD k-mesh (`mbd_rsscs_pbe`). The k path already
   refuses `nonperiodic`, so no extra cell guard is needed.

No public surface change: `GpwDft::mbd` + `k_mesh` + `gpw_stress` /
`GPWEvaluator{compute_stress}` already exist; the throw simply becomes a
working fold.

## Tests

- `tests/cpu/gpw/test_cpu_gpw_stress.cpp`: `MbdMultiKStressThrows` asserts
  the refusal this spec deletes → **replace** with
  `MbdMultiKAnalyticMatchesFd` (cheap periodic H₂, mesh 1×1×2, `mbd = true`,
  periodic MBD 2×2×2; analytic vs multi-SCF strain-FD ≤ 5e-5 Ha/Bohr³ —
  same fixture + tier as `KpointAnalyticMatchesFdDiagonal`).

## Regression (e2e) — `regressions/gpw/gpw-mbd-kpoint-stress/`

Periodic water dimer, 12 Å cube (the `gpw-mbd-periodic` / `gpw-mbd-stress`
twin geometry), PBE, DZVP-MOLOPT-SR-GTH, GTH-PBE, **electronic k mesh
1×1×2**, MBD@rsSCS free-atom with **MBD k-mesh 2×2×2**. Gates:

| Gate | Oracle | Tol |
| --- | --- | --- |
| A. \(E_{\mathrm{disp}}\) from the multi-k SCF payload equals the **Γ pin** (k-independence cross-check) | offline pymbd pin of `gpw-mbd-periodic` (−1.4965720704989138e−3 Ha) | 1e-8 Ha |
| B. \(E_{\mathrm{disp}}\) **under strain** (±h per axis, fixed fractional coords) | offline **pure-Python pymbd** `mbd_energy` (libmbd `88069c6`) at each strained geometry — new pins in REFERENCE.md | 1e-8 Ha |
| C. analytic multi-k stress vs multi-SCF strain-FD (`gpw_stress_fd`) | internal dual gate | 2e-5 Ha/Bohr³ (Γ MBD reg tier) |
| D. \(\sigma(\mathrm{mbd}) − \sigma(\mathrm{plain})\) equals the externally pinned MBD strain-FD term \(-\Delta E_{\mathrm{disp}}/(2h\Omega)\) | pins from gate B | 1e-9 Ha/Bohr³ |
| E. off-diagonals exactly 0; multi-k payload present (`n_kpts ≥ 1`) | structural | — |

Pin provenance: pymbd is run **offline at pin time** (never at ctest);
reproduce note + script under
`tools/cp2k_oracle/inputs/gpw-mbd-kpoint-stress/`. A pin is a number pymbd
printed — never seeded from ATV.

## Perf

`perf_exempt` — CP2K 2026.1 stock `PAIR_POTENTIAL TYPE` has no MBD, so there
is no fair oracle wall for an MBD stress run (same policy as
`gpw-mbd-stress` / `gpw-mbd-kpoint-forces`). Product wall gated ≤ 60 s.

## Acceptance

- Multi-k `gpw_stress` with `mbd = true` no longer throws and passes the
  dual gate (unit + regression).
- All five regression gates green; Γ regressions untouched and green
  (`reg_gpw-mbd-stress` unaffected — Γ path not edited).
- Scoped build + test only; builds wiped on exit.

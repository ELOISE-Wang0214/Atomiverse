# Spec: `pw-lda-si`

| Field | Value |
| ----- | ----- |
| Slug | `pw-lda-si` |
| Branch | `wf/pw-lda-si` |
| Domain | plane-wave DFT (spin-unpolarized LDA) |
| Observable | total energy \(E_\mathrm{tot}\) [Ry] and Frame energy [eV] |
| Geometry | Si diamond, 2-atom fcc primitive cell, \(a = 5.43\) Å |
| Basis / PP | NC-UPF ONCVPSP LDA (PseudoDojo) `tests/data/pw/Si_lda_nc.upf` |
| XC | LDA (engine Slater + PW92; PP generated with ONCVPSP iexc=3 PZ) |
| Settings (twin) | \(E_\mathrm{cut}=15\) Ry, Monkhorst–Pack \(k_\mathrm{mesh}=2\) (\(2\times2\times2\)) |
| Public API | `cpu::pw::PW` members + `cpu::pw::PWEvaluator` / `cpu::single_point` (no public `*Config`) |
| Regression | `regressions/pw/pw-lda-si/` |
| Queue row | `docs/developer/cp2k-method-queue.md` `## in_progress` `[~] pw-lda-si` |

## Goal

**Re-pin** product plane-wave LDA bulk silicon (k-mesh energy) under the CP2K
capability loop: one offline ABINIT-class / literature pin in
`regressions/pw/pw-lda-si/REFERENCE.md`, one public-API e2e regression asserting
\(E_\mathrm{tot}\), dual-gate (correctness + perf). Enables follow-on
`pw-analytic-stress-si`. Engine path already ships (`cpu::pw::PW` + UPF + MP
mesh); this slug is primarily pin + regression productization, not a redesign.

## Non-goals

- Analytic stress / forces / GGA-PBE / metallic smearing / spin / BOMD / bands
  (other queue slugs).
- Growing public `*Config` bags or `make_*` / `build_*` factories.
- Live CP2K / ABINIT / QE at ctest time.
- MPI multi-rank scaling (later `mpi-pw-scaling-4rank`).
- Production-converged ecut (~60 Ry ONCV) if it exceeds the ≤ 60 s regression
  budget — twin settings stay tutorial-class with an explicit incomplete-basis
  window vs the converged ABINIT-class limit.

## Tasks

1. **Write this spec + acceptance** on `wf/pw-lda-si` before production edits.
2. **Offline pin** in `regressions/pw/pw-lda-si/REFERENCE.md`:
   - PseudoDojo ONCVPSP LDA Si UPF provenance.
   - ABINIT-class same-pseudo deltafactor \(E_\mathrm{min}\) from PseudoDojo
     `Si.djrepo` (cite version / URL / value).
   - Twin settings (a, ecut, k-mesh, tol) and absolute tolerance with rationale
     (incomplete basis vs converged ABINIT-class \(E_\mathrm{min}\)).
   - Optional CP2K Quickstep twin input under `tools/cp2k_oracle/inputs/`
     (note: CP2K GTH-PADE absolute energies are **not** comparable; EOS only).
3. **E2E regression** `regressions/pw/pw-lda-si/reg_pw-lda-si.cpp`:
   public headers only; Si diamond Frame; `PW` + `PWEvaluator` +
   `cpu::single_point`; assert converged + \(E_\mathrm{tot}\) [Ry] vs pin
   + Frame [eV] consistent with Ry→eV conversion.
4. **Wire CMake**: `atv_add_regression(pw pw-lda-si)` under `ATV_BUILD_PW`;
   compile-def for UPF path to in-tree `tests/data/pw/Si_lda_nc.upf`.
5. **Build + run** until PASS (`ATV_BUILD_REGRESSIONS=ON`, `ATV_BUILD_PW=ON`).
6. **Perf gate**: fair wall vs CP2K if binary available; else
   `perf_exempt=true` with reason in REFERENCE.md.
7. Commit on `wf/pw-lda-si` only (never `dev`, never push).

## Acceptance criteria

See `acceptance.md` (binding).

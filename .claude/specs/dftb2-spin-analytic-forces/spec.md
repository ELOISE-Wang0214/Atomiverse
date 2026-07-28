# dftb2-spin-analytic-forces — collinear open-shell SCC-DFTB2 analytic forces

Status: done
acceptance_met: true
Loop: tb-method-loop (wave 2)
Queue row: `dftb2-spin-analytic-forces` | Closed-form collinear open-shell nuclear F
(replace re-SCC FD) | dftb2/spin | O2 off-eq | F vs FD | **oracle=dftb+|literature**
Branch: `wf/dftb2-spin-analytic-forces`

## Problem

Wave-1 landed **collinear spin-DFTB2 energy** (`multiplicity`, dual α/β Fock,
Mulliken m; `reg_spin-dftb2-o2`) and **restricted** closed-form nuclear forces
on MioSketch (`reg_dftb2-analytic-forces`). Open-shell forces still use
**re-SCC central finite difference** of the full energy (6N SCC solves) in
`launch_mio` when `occ.unrestricted`. That is a product-depth gap: MD / Minimize
on radicals (O₂, OH, NO, …) need one-shot analytic F at the SCC density.

## Design (ARCH LAW compliant)

- No new public type, no `*Config`, no `make_*` / `build_*`.
- Keep `cpu::semiempirical::DFTB2` members: `multiplicity`, `with_forces`, …
- Extend private unrestricted SCC path to retain **α/β** eigenvectors and
  eigenvalues after the last SCC cycle (same role as restricted `C`/`eps`).
- New private `analytic_forces_unrestricted` (Hellmann–Feynman + Pulay + fixed
  Δq SCC γ + repulsive). Spin channel enters only via spin-dependent onsite
  shifts \(\varepsilon_A^\sigma = \sum_B\gamma_{AB}\Delta q_B + \sigma W_A m_A\)
  on \(\partial S\); \(W_A\) is geometry-independent so fixed-\(m\) spin energy
  \(-\tfrac12\sum W_A m_A^2\) contributes **no** extra nuclear term.
- Formula boundary matches DFTB+ / Hourahine force module for collinear
  spin-polarized SCC-DFTB2 (Elstner + Köhler spin channel).
- `sk_set = ThreeOb` remains closed-shell only (unchanged throw).
- Replace the open-shell re-SCC FD branch in `launch_mio` with the analytic path.

### Energy (existing, for reference)

\[
E = \sum_{\sigma i} f_{i\sigma}\varepsilon_{i\sigma}
  - \tfrac12\sum_{AB}\gamma_{AB}\Delta q_A\Delta q_B
  - \tfrac12\sum_A W_A m_A^2
  + E_{\mathrm{rep}}
\]

### Nuclear gradient (fixed \(P^\sigma,\Delta q,m\))

\[
\frac{\partial E}{\partial R_k}
  = \sum_{\sigma}
    \mathrm{Tr}\!\big(P^\sigma(\partial_k H_0 + \varepsilon^\sigma\!\cdot\partial_k S)
    - W^\sigma \partial_k S\big)
  - \sum_{A<B}\Delta q_A\Delta q_B\,(\partial_k\gamma_{AB})
  + \partial_k E_{\mathrm{rep}}
\]

with \(P^\sigma_{\mu\nu}=\sum_i^{\mathrm{occ}\,\sigma} C_{\mu i}^\sigma C_{\nu i}^\sigma\),
\(W^\sigma_{\mu\nu}=\sum_i^{\mathrm{occ}\,\sigma}\varepsilon_i^\sigma C_{\mu i}^\sigma C_{\nu i}^\sigma\),
and atom-pair shift \(\varepsilon^\sigma_{AB}=\tfrac12(\varepsilon_A^\sigma+\varepsilon_B^\sigma)\).

## Acceptance criteria

1. Open-shell MioSketch forces use **closed-form** path (no 6N re-SCC FD in
   `launch_mio` for unrestricted).
2. O₂ triplet (multiplicity 3), off-eq bond (e.g. 1.30 Å): energy finite;
   forces finite and non-zero; translation \(\lvert\sum_i F_i\rvert\le 10^{-6}\)
   eV/Å per Cartesian; force antisymmetry on the bond axis.
3. Analytic \(F\) vs central **re-SCC** FD of energy-only (\(h=10^{-4}\) Å)
   within **1.0 eV/Å** on gated components (same class as multi-element
   restricted residual; same-element O–O should be tighter).
4. Unit test in `tests/cpu/semiempirical/test_dftb2.cpp` (single-method,
   seconds) covering open-shell analytic vs FD.
5. Regression `regressions/semiempirical/dftb2-spin-analytic-forces/` —
   public API only, `REFERENCE.md` with `oracle: dftb+|literature` and
   `oracle_version`, hard-coded gates, ≤ 60 s; CMake targets
   `reg_dftb2-spin-analytic-forces` / `run_reg_dftb2-spin-analytic-forces`.
6. Offline oracle note under `tools/tb_oracle/inputs/dftbplus/` (no live
   DFTB+ in ctest; mio-style sketch SK ⇒ formula/module-boundary pin).
7. SE-only dual gate green; `perf_exempt=true` (no DFTB+ binary on pin host;
   analytic path replaces 6N re-SCC).
8. Header / roadmap / queue prose updated: open-shell MioSketch F = analytic.

## Out of scope

- ThreeOb open shell / ThreeOb analytic F (`dftb2-analytic-forces-threeob`)
- Noncollinear spin, DFTB3 forces, analytic stress, k-mesh
- Absolute DFTB.org mio-1-1 energy identity (`dftb2-dftbplus-energy-pin`)

## Files (expected)

| Path | Role |
| ---- | ---- |
| `src/cpu/semiempirical/dftb2.cpp` | retain α/β orbs; analytic open-shell F |
| `include/.../dftb2.h` | doc: open-shell analytic on MioSketch |
| `tests/cpu/semiempirical/test_dftb2.cpp` | unit gate |
| `regressions/semiempirical/dftb2-spin-analytic-forces/` | e2e + REFERENCE |
| `regressions/CMakeLists.txt` | register reg |
| `tools/tb_oracle/inputs/dftbplus/dftb2_spin_analytic_forces_note.md` | offline pin |
| `docs/roadmap/semiempirical.md` | deepen spin F row |
| `docs/tutorial/methods/dftb2-spin-analytic-forces.md` | short method note |

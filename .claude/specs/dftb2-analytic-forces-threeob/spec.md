---
slug: dftb2-analytic-forces-threeob
status: done
created: 2026-07-27
domain: cpu/semiempirical
queue: tb-method-loop
oracle: dftb+
method: dftb2
deferred: false
depends_on: [dftb2-3ob-organic, dftb2-analytic-forces]
unblocks: [dftb2-analytic-stress]
---

# Spec: `dftb2-analytic-forces-threeob` — closed-form F on `DftbSkSet::ThreeOb`

**Slug:** `dftb2-analytic-forces-threeob`  
**Domain:** `cpu/semiempirical` (DFTB2 / 3ob-3-1)  
**Queue:** tb-method-loop (Phase A deepen after wave-1)  
**Oracle family:** **DFTB+ only** (never OpenMOPAC / xtb)

## Summary

Wave-1 landed closed-form nuclear forces on **restricted MioSketch** and a
**ThreeOb energy** path with static 3ob-3-1 SK. ThreeOb `with_forces=true`
still uses **re-SCC central FD** (6N full SCC solves). This slug replaces that
with the **same closed-form** force assembly as MioSketch, using 3ob radial
tables + repulsive splines (Ha/Bohr internally, eV/Å on Frame).

\[
F = -\nabla\big(\mathrm{Tr}(P H_0) + E_{\mathrm{SCC}} + E_{\mathrm{rep}}\big)
\]

at fixed SCC \(\Delta q\), \(P\), \(\varepsilon\) (Hellmann–Feynman + Pulay +
\(\partial\gamma/\partial R\) + \(\partial E_{\mathrm{rep}}\)). Formula boundary =
DFTB+ / Elstner (Hourahine et al. 2020). Absolute DFTB+ force numbers are
optional; dual gate is analytic vs re-SCC FD of the **same** ATV energy on
off-eq organics (H2O or CH4 stretch).

## Domain basis

- Elstner et al., Phys. Rev. B **58**, 7260 (1998).
- Hourahine et al., J. Chem. Phys. **152**, 124101 (2020) (DFTB+).
- Gaus et al., J. Chem. Theory Comput. **9**, 338 (2013) — 3ob-3-1 SK.
- Prior: `dftb2-analytic-forces` (MioSketch closed-form), `dftb2-3ob-organic`
  (CH4 energy pin vs DFTB+ 24.1).

## Design

### Public API

- Unchanged members: `sk_set`, `scc_tol`, `max_scc`, `with_forces`, …
- `sk_set = ThreeOb` + `with_forces = true` → closed-form F (not re-SCC FD).
- No public `*Config`, no `make_*` / `build_*`. Frame remains sole I/O.
- Doxygen: ThreeOb forces are analytic on closed shell (open-shell still
  rejected on ThreeOb).

### Implementation contract

After SCC on the 3ob path (Bohr geometry, Ha energy):

1. Keep final \(C\), \(\varepsilon\), \(\Delta q\).
2. SK radial channels: linear table interp; \(\partial/\partial R\) from
   segment slope \(+(b-a)/\mathrm{grid}\) (+ Slater–Koster angular as Mio).
3. \(\gamma(R)\) Ohno–Klopman in Ha/Bohr; \(d\gamma/dR = -R\,\gamma^3\).
4. Repulsive: short-range exp derivative + spline polynomial derivative.
5. Convert \(F\) Ha/Bohr → eV/Å for Frame.

### Out of scope

- Open-shell ThreeOb / spin forces (`dftb2-spin-analytic-forces`).
- DFTB3 third-order force terms.
- Cell stress (`dftb2-analytic-stress`).
- MD/MINIMIZER (SE-only cmake).

## Files

| Path | Action |
| ---- | ------ |
| `src/cpu/semiempirical/dftb2.cpp` | Analytic ThreeOb force path; drop re-SCC FD |
| `include/.../dftb2.h` | Doc: ThreeOb closed-shell forces analytic |
| `tests/cpu/semiempirical/test_dftb2_3ob.cpp` | Unit: analytic vs FD |
| `regressions/semiempirical/dftb2-analytic-forces-threeob/` | E2E + REFERENCE |
| `regressions/CMakeLists.txt` | register regression |
| `tools/tb_oracle/inputs/dftbplus/` | offline force pin note |

## Acceptance

See `acceptance.md`.

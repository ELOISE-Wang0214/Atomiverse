---
slug: dftb2-analytic-forces
status: in-progress
created: 2026-07-27
domain: cpu/semiempirical
queue: tb-method-loop
oracle: dftb+
method: dftb2
deferred: false
depends_on: []
unblocks: [dftb2-analytic-stress, dftb3-water-mio]
---

# Spec: `dftb2-analytic-forces` — analytic SCC-DFTB2 nuclear gradient

**Slug:** `dftb2-analytic-forces`  
**Domain:** `cpu/semiempirical` (DFTB2)  
**Queue:** tb-method-loop (Phase A ladder)  
**Oracle family:** **DFTB+ only** (never OpenMOPAC / xtb for this slug)

## Summary

Replace the re-SCC **central finite-difference** nuclear forces in
`cpu::semiempirical::DFTB2` with a **closed-form** gradient of the Elstner energy:

\[
E = \mathrm{Tr}(P H_0) + E_{\mathrm{SCC}} + E_{\mathrm{rep}}
\]

with SCC Hamiltonian \(H = H_0 + \tfrac12 S_{\mu\nu}(\Delta V_A + \Delta V_B)\),
\(E_{\mathrm{SCC}} = \tfrac12\sum_{AB}\Delta q_A\gamma_{AB}\Delta q_B\),
and \(E_{\mathrm{rep}} = \sum_{A<B} A_{AB}e^{-B_{AB}R}\).
(The prior `sum nε − E_SCC` form differs by a geometry-dependent
\(\sum n_A\Delta V_A\) term; Tr(PH0)+E_SCC is the Elstner/DFTB+ energy.)

Forces \(F = -\nabla E\) [eV/Å] must match central FD of the **same** energy
to a tight tolerance on off-equilibrium H2O (and unit-tested on H2). Absolute
DFTB+ numerical identity is **not** required here: in-tree SK tables are
mio-**style** sketches, not bit-identical mio-1-1 (see
`dftb2-dftbplus-energy-pin` / `dftb2-3ob-organic`). The oracle is the DFTB+
**force module boundary** (Hellmann–Feynman + Pulay + \(\partial\gamma/\partial R\)
+ repulsive) as documented in Elstner / Hourahine.

## Domain basis

- Elstner et al., Phys. Rev. B **58**, 7260 (1998) — SCC-DFTB2 energy.
- DFTB+ force construction: \(\mathrm{Tr}(P\,\partial H_0) - \mathrm{Tr}(W\,\partial S)\)
  plus SCC \(\gamma\) geometry derivative and repulsive pair derivative
  (Hourahine et al., J. Chem. Phys. **152**, 124101 (2020); DFTB+ docs).
- Existing ATV foothold: `DFTB2::launch` energy path + FD forces
  (`src/cpu/semiempirical/dftb2.cpp`); header brief already claims “analytic”
  while the body is re-SCC FD — **doc/implementation debt**.

## Design

### Public API

- Unchanged: `DFTB2` members `scc_tol`, `max_scc`, `with_forces`; Frame I/O;
  no `*Config`, no `make_*` / `build_*`.
- `with_forces = true` writes \(F=-\nabla E\) via the closed-form path (not FD).
- Comment / Doxygen: “analytic (Hellmann–Feynman + Pulay + SCC \(\gamma\) + repulsive)”.

### Force assembly (implementation contract)

After SCC convergence (fixed \(\Delta q\), \(P\), eigenvectors, \(\varepsilon\)):

1. Energy-weighted density \(W_{\mu\nu}=\sum_i n_i\varepsilon_i C_{\mu i}C_{\nu i}\).
2. Pair SK derivatives \(\partial H_0/\partial R\), \(\partial S/\partial R\) from the
   same radial + Slater–Koster angular forms used for energy.
3. Electronic: \(\sum_{\mu\nu}\big[P_{\mu\nu}(\partial H_0 + \tfrac12(\Delta V_A+\Delta V_B)\partial S)
   - W_{\mu\nu}\partial S\big]\) plus fixed-\(\Delta q\) SCC shift
   \(\sum_A q_A\partial(\Delta V_A)/\partial R - \partial E_{\mathrm{SCC}}/\partial R\).
4. Repulsive: \(\partial_R\big(A e^{-BR}\big)\).
5. \(F_A = -\partial E/\partial R_A\); net force \(\sum_A F_A \approx 0\).

### Out of scope

- DFTB0 / DFTB3 / spin / PBC stress (other slugs).
- Replacing sketch SK with true mio-1-1 / 3ob.
- Absolute energy pin vs live DFTB+ binary (no DFTB+ module on this host;
  numerical twin deferred to `dftb2-dftbplus-energy-pin`).
- MD/MINIMIZER product paths (SE-only cmake keeps them OFF).

## Files

| Path | Action |
| ---- | ------ |
| `src/cpu/semiempirical/dftb2.cpp` | Replace FD force block with analytic assembly |
| `include/.../dftb2.h` | Doc: analytic, not FD |
| `tests/cpu/semiempirical/test_dftb2.cpp` | Unit: analytic vs FD on H2 / H2O component |
| `regressions/semiempirical/dftb2-analytic-forces/` | E2E Evaluate path + FD gate + REFERENCE |
| `regressions/CMakeLists.txt` | `atv_add_regression(... dftb2-analytic-forces)` |

## Acceptance

See `acceptance.md`. Dual gate: correctness (FD match) + perf_exempt if no
DFTB+ wall available for fair ratio.

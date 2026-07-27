---
slug: spin-dftb2-o2
status: done
created: 2026-07-27
domain: cpu/semiempirical
queue: tb-method-loop
method: dftb2
oracle: dftb+|literature
deferred: false
depends_on: []
unblocks: [noncollinear-spin-dftb]
grilled: false
---

# Spec: `spin-dftb2-o2` — collinear spin-polarized SCC-DFTB2 (O₂)

## Summary

Extend `cpu::semiempirical::DFTB2` with **collinear unrestricted** spin
polarization so open-shell molecules (O₂ triplet) have finite SCC energy,
Mulliken charge Δq, and Mulliken spin population m = Nα−Nβ. Deliver via
**Evaluate + Semiempirical** only (SE-only cmake). Family oracle = **DFTB+**;
absolute energy pin is the ATV mio-style Hamiltonian twin (SK not
DFTB.org-bit-identical — same hygiene as `h2o_dftb2_mio`) with literature
method citations.

## Domain basis

Collinear spin-DFTB (Köhler, Frauenheim et al.; DFTB+ spin channel):

- Valence electrons → \(n_\alpha, n_\beta\) from charge + multiplicity \(2S+1\).
- Separate densities \(P^\alpha, P^\beta\); Mulliken \(q_A^\sigma\), \(m_A = q_A^\alpha - q_A^\beta\), \(\Delta q_A = q_A - n_A^{\mathrm{val}}\).
- Hamiltonian \(H_{\mu\nu}^\sigma = H^0_{\mu\nu} + \tfrac12 S_{\mu\nu}(\varepsilon_A^\sigma + \varepsilon_B^\sigma)\) with
  \(\varepsilon_A^\sigma = \sum_B \gamma_{AB}\Delta q_B + \sigma\, W_A m_A\) (\(\sigma = \pm 1\)).
- Energy: \(E = \sum_{\sigma,i} f_{i\sigma}\varepsilon_{i\sigma} - \tfrac12\sum_{AB}\gamma_{AB}\Delta q_A\Delta q_B - \tfrac12\sum_A W_A m_A^2 + E_{\mathrm{rep}}\).
- Closed shell (`multiplicity == 1`, even \(N_e\)) keeps the existing restricted path bit-stable.

Static atom-averaged spin constants \(W_A\) live on `DftbElementRow` (sketch
order of magnitude; no runtime `.skf`).

## Public API (no `*Config`, no factories)

```cpp
cpu::semiempirical::DFTB2 dftb;
dftb.charge = 0;           // [e]
dftb.multiplicity = 3;     // O2 triplet
dftb.scc_tol = 1e-10;
dftb.with_forces = false;  // energy / spin gate first
eval.on(evaluate::Phase::evaluate).add(std::move(dftb));
// after run: dftb.last_energy(), last_mulliken_dq(), last_mulliken_spin()
```

Members only; parameters are constructor defaults + public fields.

## Acceptance

| Gate | Criterion |
| ---- | --------- |
| Unit | O₂ triplet energy finite; \(\lvert\sum_A m_A - 2\rvert < 0.2\); closed-shell H₂ still finite |
| Unit | Odd-electron or invalid multiplicity throws; unsupported Z still throws |
| Regression | `reg_spin-dftb2-o2`: O₂ @ 1.21 Å, mult=3; finite E; hard pin ±0.05 eV; \(\sum m \approx 2\); per-atom \(\lvert m\rvert > 0.3\) |
| REFERENCE | `oracle: literature` (+ family DFTB+); method citations; geometry; pin numbers; units eV |
| Build | SE-only cmake; scoped `--target atomiverse test_cpu_dftb2 test_cpu_dftb2_sk reg_spin-dftb2-o2` |
| Perf | `perf_exempt` if no offline DFTB+ binary (no fair wall); else \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\le 10\) |

## Out of scope

- Noncollinear spin, SOC, DFTB+U, TD-DFTB
- Analytic spin-DFTB forces (FD re-SCC ok if `with_forces`)
- True mio-1-1 / DFTB+ absolute energy identity (needs SK fidelity slug)
- MD/MINIMIZER/HF cmake ON

## Files

- `.claude/specs/spin-dftb2-o2/spec.md` (this)
- `include/.../dftb2.h`, `src/.../dftb2.cpp`
- `include/.../basis_data/dftb_mio_tables.h` (`spin_w`)
- `tests/cpu/semiempirical/test_dftb2.cpp` (spin cases)
- `regressions/semiempirical/spin-dftb2-o2/{reg_*.cpp,REFERENCE.md}`
- `regressions/CMakeLists.txt` registration

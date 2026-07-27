---
slug: dftb2-analytic-stress
status: in-progress
created: 2026-07-27
domain: cpu/semiempirical
queue: tb-method-loop
oracle: dftb+
method: dftb2
deferred: false
depends_on: [dftb2-analytic-forces]
unblocks: [dftb2-periodic-gamma, dftb3-analytic-forces]
---

# Spec: `dftb2-analytic-stress` — analytic SCC-DFTB2 cell stress

**Slug:** `dftb2-analytic-stress`  
**Domain:** `cpu/semiempirical` (DFTB2)  
**Queue:** tb-method-loop (Phase A ladder — DFTB2 depth)  
**Oracle family:** **DFTB+ only** (never OpenMOPAC / xtb for this slug)

## Summary

Ship a **closed-form cell stress** for restricted MioSketch SCC-DFTB2, writing
intensive Voigt stress to Frame `kStress` [eV/Å³]:

\[
\sigma_{\alpha\beta} = -\frac{1}{\Omega}\frac{\partial E}{\partial\varepsilon_{\alpha\beta}}
\]

at the SCC density (fixed \(\Delta q\), \(P\), \(W\)), reusing the same pair
Hellmann–Feynman + Pulay + \(\partial\gamma/\partial R\) + repulsive terms as
`dftb2-analytic-forces`. Gate against multi-SCC **affine strain-FD** of the
same total energy (DFTB+ stress module boundary / mechanical convention as
PW/GPW). Absolute DFTB+ numerical identity is **not** required: in-tree SK is
mio-**style** sketch, not bit-identical mio-1-1.

## Domain basis

- Elstner et al., Phys. Rev. B **58**, 7260 (1998) — SCC-DFTB2 energy.
- DFTB+ stress construction (Hourahine et al., J. Chem. Phys. **152**, 124101
  (2020)): pair contraction of electronic + SCC \(\gamma\) + repulsive
  derivatives under affine strain of nuclear coordinates + cell.
- ATV convention: intensive \(\sigma\) [eV/Å³] Voigt order
  `(xx,yy,zz,xy,xz,yz)` — same as PW `with_stress` / GPW `compute_stress`
  (`σ = −(1/Ω) ∂E/∂ε`).
- Prerequisite: restricted MioSketch analytic forces already ship
  (`reg_dftb2-analytic-forces`).

## Design

### Public API

- New member `bool with_stress = false` on `cpu::semiempirical::DFTB2`.
- When true (restricted MioSketch): after SCC, accumulate analytic stress into
  `frame.meta.get(kStress)` [eV/Å³] (additive; Evaluate's `ZeroStress` zeros first).
- Optional `last_stress()` → `std::array<double,6>` of last intensive Voigt.
- No `*Config`, no `make_*` / `build_*`. Frame remains sole I/O currency.
- Open-shell and ThreeOb: analytic stress **out of scope** this slug; if
  `with_stress` is set on those paths, fall back to multi-SCC strain-FD of \(E\)
  (or document/throw — prefer strain-FD for product completeness on open shell,
  analytic only for restricted Mio).

### Stress assembly (implementation contract)

After SCC (restricted MioSketch), for each pair \(A<B\) with bond vector
\(\mathbf{R}=\mathbf{R}_B-\mathbf{R}_A\) and pair force contribution \(\mathbf{F}\)
on \(A\) (same \(dE_{\mathrm{bond}}\) as nuclear forces):

1. \(\partial E/\partial\varepsilon_{\alpha\alpha} \mathrel{+}= F_\alpha R_\alpha\) (diag).
2. Off-diagonal (symmetric strain parameter \(h\) with \(\varepsilon_{\alpha\beta}=\varepsilon_{\beta\alpha}=h\)):
   \(\partial E/\partial h \mathrel{+}= F_\alpha R_\beta + F_\beta R_\alpha\).
3. \(\sigma = -(\partial E/\partial\varepsilon)\,/\,\Omega\).

Electronic SK / Pulay / \(\gamma\) / repulsive content lives **inside** \(\mathbf{F}\)
(already assembled by analytic forces). Do **not** re-introduce re-SCC per strain
axis on the restricted Mio product path.

### Out of scope

- True lattice-sum periodic DFTB stress product depth (`dftb2-periodic-gamma` /
  `dftb2-kpoints`).
- DFTB3 / spin-analytic stress deepen slugs.
- Absolute mio-1-1 DFTB+ stress numbers (sketch SK).
- MD/NPT product (SE-only cmake keeps MD OFF).

## Files

| Path | Action |
| ---- | ------ |
| `include/.../dftb2.h` | `with_stress`, `last_stress()`, Doxygen |
| `src/cpu/semiempirical/dftb2.cpp` | Analytic stress contraction; write `kStress` |
| `tests/cpu/semiempirical/test_dftb2.cpp` | Unit: analytic σ vs strain-FD on H2 |
| `regressions/semiempirical/dftb2-analytic-stress/` | E2E Evaluate + FD gate + REFERENCE |
| `regressions/CMakeLists.txt` | register regression |
| `tools/tb_oracle/inputs/dftbplus/dftb2_analytic_stress_note.md` | offline pin note |

## Acceptance

See `acceptance.md`. Dual gate: correctness (strain-FD match) +
`perf_exempt=true` if no DFTB+ wall on pin host.

# Reference — collinear spin-DFTB2 analytic nuclear forces (O₂)

```text
oracle: dftb+
oracle_version: formula-level (Hourahine et al. JCP 152, 124101 (2020);
  Elstner PRB 58, 7260 (1998); Köhler collinear spin-DFTB; no DFTB+ binary
  on pin host)
family_oracle: dftb+
```

## Method

Collinear unrestricted SCC-DFTB2:

\[
E = \sum_{\sigma i} f_{i\sigma}\varepsilon_{i\sigma}
  - \tfrac12\sum_{AB}\gamma_{AB}\Delta q_A\Delta q_B
  - \tfrac12\sum_A W_A m_A^2
  + E_{\mathrm{rep}}
\]

Nuclear forces are the **closed-form** gradient at the SCC density (dual α/β
Hellmann–Feynman + Pulay energy-weighted density + spin-dependent Hamiltonian
shifts on \(\partial S\) + fixed-\(\Delta q\) \(-\partial E_{\mathrm{SCC}}\) +
repulsive) — **not** a re-SCC finite difference. \(W_A\) is geometry-independent;
fixed-\(m\) spin energy contributes no extra nuclear term.

This is the DFTB+ collinear spin force module boundary.

## Hamiltonian / tables

- ATV class: `cpu::semiempirical::DFTB2`
- SK / Hubbard / spin_W / repulsive: static `basis_data/dftb_mio_tables.h`
  (mio-**style** sketch H/C/N/O). **Not** bit-identical to DFTB.org mio-1-1.

## Geometry (Å)

O₂ triplet, **stretched** bond so \(|F|\neq 0\):

| atom | Z | x | y | z |
| ---- | - | - | - | - |
| O | 8 | 0.00 | 0.00 | 0.00 |
| O | 8 | 1.30 | 0.00 | 0.00 |

Isolated molecule, cubic box 30 Å. `charge = 0`, `multiplicity = 3`
(\(n_\alpha=7\), \(n_\beta=5\)), `scc_tol = 1e-10`.

## Pins (hard-coded gates)

1. **Energy:** finite, \(|E| > 10^{-3}\) eV.
2. **Forces:** all components finite; \(|F|^2 > 10^{-8}\) at off-eq geometry.
3. **Translation:** \(|\sum_i F_i| \le 10^{-6}\) eV/Å per Cartesian.
4. **Action–reaction:** \(F_{0,x} + F_{1,x} = 0\) (within \(10^{-6}\)).
5. **FD identity:** analytic \(F\) vs central **re-SCC** FD of the same \(E\)
   (\(h=10^{-4}\) Å) on \(F_{0,x}\) and \(F_{1,x}\) within **0.05 eV/Å**
   (same-element O–O; multi-element residual class remains ~1 eV/Å).

Offline pin note:
`tools/tb_oracle/inputs/dftbplus/dftb2_spin_analytic_forces_note.md`.

## Perf

```text
perf_exempt: true
perf_exempt_reason: no DFTB+ binary / fair wall on pin host; analytic path
  replaces 6N full re-SCC force evaluations with one post-SCC contraction
```

## Units

Public ATV: energy **eV**, positions **Å**, forces **eV/Å**.

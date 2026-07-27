# spin-dftb2-o2 — collinear spin-polarized SCC-DFTB2 on O₂

```text
oracle: literature
oracle_version: Elstner1998+Kohler2001-method; ATV-mio-style-twin-2026-07-27
family_oracle: dftb+
```

**Family numerical oracle** for DFTB is **DFTB+**. Absolute energy identity
with a live DFTB+ binary is **not** claimed here: ATV static SK tables are
mio-**style** sketches (`basis_data/dftb_mio_tables.h`), not DFTB.org
mio-1-1 bit-identical. This regression pins the **ATV Hamiltonian twin**
and physics gates (triplet Mulliken spin). Offline DFTB+ pins land when
SK fidelity ships (`dftb2-3ob-organic` / `dftb2-dftbplus-energy-pin`).

No live OpenMOPAC / DFTB+ / xtb inside ctest.

## Method

Collinear unrestricted SCC-DFTB2:

- Charge SCC: \(\gamma_{AB}\) (Ohno–Klopman-like) + Hubbard \(U\) (Elstner
  et al., *Phys. Rev. B* **58**, 7260 (1998),
  DOI [10.1103/PhysRevB.58.7260](https://doi.org/10.1103/PhysRevB.58.7260)).
- Spin channel: atom-averaged \(W_A\) shifts
  \(\varepsilon_A^\sigma = \sum_B\gamma_{AB}\Delta q_B + \sigma W_A m_A\)
  (Köhler / Frauenheim collinear spin-DFTB; DFTB+ spin module shape).
- Energy:
  \(E = \sum_{\sigma i} f_{i\sigma}\varepsilon_{i\sigma}
  - \tfrac12\sum_{AB}\gamma_{AB}\Delta q_A\Delta q_B
  - \tfrac12\sum_A W_A m_A^2 + E_{\mathrm{rep}}\) [eV].

Public API: `cpu::semiempirical::DFTB2` with members `charge`,
`multiplicity`, `scc_tol`, `max_scc`, `with_forces`. Mulliken
\(\Delta q\) / \(m\) via `last_mulliken_dq()` / `last_mulliken_spin()`.

## Geometry

O₂, bond length 1.21 Å (near experimental), cubic box 30 Å (isolated):

| Atom | Z | x [Å] | y | z |
|------|---|-------|---|---|
| O    | 8 | 0.00  | 0 | 0 |
| O    | 8 | 1.21  | 0 | 0 |

Charge 0, multiplicity **3** (triplet ground state; \(n_\alpha=7\),
\(n_\beta=5\) for 12 valence electrons).

SCC tol \(10^{-10}\), max 80 iterations. Forces off (energy + spin gate).

## Pins / gates

| Quantity | Value | Tol | Unit |
| -------- | ----- | --- | ---- |
| Total energy \(E\) | **−142.1382039484** | 0.05 | eV |
| \(\sum_A m_A\) (Mulliken spin) | 2.0 | 0.25 | e |
| Per-atom \(\lvert m\rvert\) | > 0.3 each | — | e |
| \(\lvert m_0 - m_1\rvert\) | < 0.2 | — | e |

Observed (ATV twin, 2026-07-27): \(m = [1.0, 1.0]\), \(\Delta q = [0, 0]\).

## Perf

```text
perf_exempt: true
perf_exempt_reason: no DFTB+ binary on host; absolute energy is ATV SK twin,
  not a fair wall against family oracle
```

## Units

Public ATV: energy **eV**, positions **Å**, Mulliken populations **e**.

# Reference — extended-Lagrangian (iEL/SCF) PIM dipole dynamics (map §2)

oracle: metalwalls
oracle_version: MetalWalls JOSS 5, 2373 (2020) — map §2 "XL dipoles"; no MW binary on agent host → literature-form pins (iEL/SCF: Albaugh, Demerdash & Head-Gordon, JCP 143, 174104 (2015); recursion: Niklasson et al., JCP 130, 214109 (2009), Table I)
perf_exempt: true
perf_exempt_reason: No MetalWalls binary on PATH for a fair XL/mass-zero dipole wall twin, and LAMMPS has no iEL/SCF induced-dipole twin (its polarizable paths are Drude/QEq). Same pin path as pim-induced-dipole-bulk.

## Physics

Extended-Lagrangian induced-dipole dynamics for the bulk Polarizable Ion
Model (`cpu::PimXlDipole`): an auxiliary dipole set \(\boldsymbol{\mu}^{aux}\)
is propagated by the time-reversible dissipative Niklasson recursion

\[
\mu^{aux}_{n+1} = 2\mu^{aux}_n - \mu^{aux}_{n-1}
  + \kappa\,(\mu^{SCC}_n - \mu^{aux}_n)
  + \alpha \sum_{k=0}^{K} c_k\, \mu^{aux}_{n-k}
\]

(Table-I coefficients, K = 5) and serves as the SCC initial guess each MD
step — the iEL/SCF scheme, i.e. the classical-polarization analogue of
XL-BOMD. Energy and forces on Frame are the owned `PimInducedDipole`'s
SCC-converged Hellmann–Feynman values, so the scheme changes the *guess
dynamics*, never the physics. MetalWalls' own mass-zero route is MaZe
(SHAKE-like constrained propagation) — deferred exactly as `ConstPotMaZe`
defers it for electrode charges.

## Pins

### 1. Stationary two-site fixed point (literature-form, hard-coded)

Geometry of `pim-induced-dipole-bulk/REFERENCE.md` (q = ±1 e,
α = 1.0/1.5 Å³, r = 3 Å, thole_a = 0.39, L = 50 Å). The XL fixed point at a
stationary geometry IS the SCF solution, cold and warm launches alike:

| Observable | Value | Tol |
| ---------- | ----- | --- |
| \(E_{\mathrm{pol}}\) [eV] | **−0.12040946414707571** | 1e-9 |
| \(F_{A,x}\) [eV/Å] | **+0.040341072868572643** | 1e-8 |
| \(F_{B,x}\) [eV/Å] | **−0.040341072868572643** | 1e-8 |
| warm SCC iterations | < cold, and ≤ 2 | exact |
| warm aux deviation [e·Å] | < 1e-10 | — |

Independent closed-form re-derivation (no Atomiverse code):
`tools/metalwalls_oracle/inputs/metalwalls/extended-lagrangian-dipoles/pin_two_site.py`
(agrees to ~1e-14 relative).

### 2. Dynamics — defining iEL/SCF properties (self-consistency gates)

Bulk polarizable toy: 8 ions, perturbed 2×2×2 rocksalt-parity sub-lattice
(spacing 4 Å) centered in L = 20 Å; masses 22.99/35.45 amu; velocity-Verlet
150 × 0.5 fs from rest; PIM cutoff 10 Å (no cutoff crossings → smooth,
conservative toy force field); SCC tol 1e-12 in both twins.

| Gate | Bound |
| ---- | ----- |
| \(\max_t \max_i \|\mu^{XL}_i(t) - \mu^{SCF}_i(t)\|\) | < 1e-6 e·Å |
| total warm SCC iterations | < cold-twin total (strict) |
| \(\|E_{tot}(t_{end}) - E_{tot}(0)\|\), XL twin | < 5e-6 eV |
| \(\|E_{tot}(t_{end}) - E_{tot}(0)\|\), SCF twin | < 5e-6 eV |

\(E_{tot} = \mathrm{KE} + E_{\mathrm{pol}}\) [eV].

**API:** `cpu::PimXlDipole{cutoff, alpha_by_type, thole_a, xl_order, tol,
max_iter}` — public headers only; direct `Context` + `build`/`launch`
(maze-charge-dynamics precedent).

## Provenance

- Albaugh, Demerdash & Head-Gordon, *J. Chem. Phys.* **143**, 174104 (2015) —
  iEL/SCF extended-Lagrangian guess dynamics for classical mutual induction.
- Niklasson et al., *J. Chem. Phys.* **130**, 214109 (2009), Table I —
  dissipative recursion coefficients (κ, α, c_k).
- Niklasson, *Phys. Rev. Lett.* **100**, 123004 (2008) — time-reversible XL.
- Coretti, Bonella & Ciccotti, *J. Chem. Phys.* **149**, 191102 (2018) —
  MaZe mass-zero constrained dynamics for polarizable models (MW's route).
- Thole, *Chem. Phys.* **59**, 341 (1981); Applequist, Carl & Fung, *JACS*
  **94**, 2952 (1972) — damped point-dipole SCC (stationary pin).
- MetalWalls JOSS **5**, 2373 (2020) — capability map §2 "XL dipoles".

## Fair wall (perf gate)

No MetalWalls binary available on the agent host → **perf_exempt**. When an
MW twin appears, record \(T_{\mathrm{atv}}/T_{\mathrm{MW}}\) here and require
\(\le 10\).

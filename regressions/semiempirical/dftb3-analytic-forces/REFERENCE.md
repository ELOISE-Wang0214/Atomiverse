# Reference — DFTB3 analytic nuclear forces (H2O)

```text
oracle: dftb+
oracle_version: 24.1 — formula level only (no DFTB+ binary on the pin host)
```

## Method

Third-order self-consistent-charge DFTB (**DFTB3**):

- Gaus, Cui, Elstner, *J. Chem. Theory Comput.* **7**, 931 (2011),
  DOI [10.1021/ct100684s](https://doi.org/10.1021/ct100684s) — the E⁽³⁾ /
  Γ_AB formulation and its nuclear gradient
- Second-order backbone: Elstner et al., *Phys. Rev. B* **58**, 7260 (1998),
  DOI [10.1103/PhysRevB.58.7260](https://doi.org/10.1103/PhysRevB.58.7260)
- Force-module boundary: DFTB+ 24.1, Hourahine et al., *J. Chem. Phys.*
  **152**, 124101 (2020), DOI [10.1063/1.5143190](https://doi.org/10.1063/1.5143190)

Family-locked oracle package is **DFTB+** (never OpenMOPAC, never xtb).

## What changed on this branch

`cpu::semiempirical::DFTB3` used to produce forces by **re-SCC central finite
difference** — `6N` extra full self-consistency solves per force evaluation.
That path is deleted. Forces are now the **closed-form gradient at the converged
SCC density**:

\[
\frac{\partial E}{\partial R}
= \sum_{\mu\nu} P_{\mu\nu}\frac{\partial H^0_{\mu\nu}}{\partial R}
+ \sum_{\mu\nu} P_{\mu\nu}\tfrac12(V_A+V_B)\frac{\partial S_{\mu\nu}}{\partial R}
- \sum_{\mu\nu} W_{\mu\nu}\frac{\partial S_{\mu\nu}}{\partial R}
+ \frac{\partial (E^{(2)}+E^{(3)})}{\partial R}\Big|_{\Delta q}
+ \frac{\partial E_{\mathrm{rep}}}{\partial R}
\]

with \(W_{\mu\nu}=\sum_o n_o\varepsilon_o C_{\mu o}C_{\nu o}\) and
\(V=\partial(E^{(2)}+E^{(3)})/\partial\Delta q\). The explicit-geometry SCC term
per unordered pair \(A\neq B\), bond \(\mathbf{R}=\mathbf{R}_B-\mathbf{R}_A\):

\[
\frac{\partial}{\partial R}\Big[\gamma_{AB}\Delta q_A\Delta q_B
 + \tfrac13\big(\Delta q_A^2\Delta q_B\Gamma_{AB}
              + \Delta q_B^2\Delta q_A\Gamma_{BA}\big)\Big]
\]

\(\Gamma_{AB}\neq\Gamma_{BA}\), so **both** third-order orderings carry a radial
derivative — the same asymmetry that makes \(V^{(3)}\) need two contractions.
For the in-tree DFTB+ Elstner kernel
(\(\gamma_{AB}=(R^2+\tau^2)^{-1/2}\), \(\tau=\tfrac12(1/U_A+1/U_B)\)):

\[
\frac{\partial\gamma_{AB}}{\partial R}=-R\gamma_{AB}^3,\qquad
\frac{\partial\Gamma_{AB}}{\partial R}
 = -\frac{3R\tau}{(R^2+\tau^2)^{5/2}}\cdot\frac{1}{2U_A^2}\,U^{\mathrm{d}}_A
\]

On-site \(\gamma_{AA}=U_A\) and \(\Gamma_{AA}=U^{\mathrm{d}}_A\) carry no
geometry ⇒ zero gradient.

**Why no coupled-perturbed term.** The SCC Hamiltonian
\(H=H^0+\tfrac12 S(V_A+V_B)\) is exactly \(\partial E/\partial P_{\mu\nu}\)
(Mulliken \(\partial\Delta q_A/\partial P_{\mu\nu}
=\tfrac12 S_{\mu\nu}(\delta_{\mu\in A}+\delta_{\nu\in A})\)), so the converged
density is stationary and the fixed-Δq expression above **is** the total
derivative. The FD identity below is the experiment that proves it.

## Hamiltonian / tables

- ATV class: `cpu::semiempirical::DFTB3` (public members only; no `*Config`)
- SK / Hubbard \(U\) / \(U^{\mathrm{d}}\) / repulsive: static
  `basis_data/dftb_mio_tables.h` + transcribed mio-1-1 H–H
  (`internal/mio_hh_table.h`). mio-**style** sketch for H/C/N/O — **not**
  bit-identical to DFTB.org mio-1-1.
- The SK block and its bond-vector derivative come from ONE shared kernel
  (`internal/mio_sk_block.h`, `MioSketchPair`), also used by DFTB2, so the
  gradient is the derivative of the H⁰/S the SCC actually converged.

## Geometry (Å)

| atom | Z | x | y | z |
| ---- | - | - | - | - |
| O | 8 | 0.00 | 0.00 | 0.00 |
| H | 1 | 1.05 | 0.00 | 0.00 |
| H | 1 | −0.26 | 1.00 | 0.00 |

Isolated molecule, cubic box 20 Å (minimum-image molecular path).
`scc_tol = 1e-10`, `max_scc = 50`, `mix_beta = 0.3`, `with_forces = true`.

## Pins (hard-coded gates)

### 1. FD identity — the physics oracle (formula-exact, no binary needed)

Analytic \(F\) vs central **re-SCC** finite difference of the same total energy
(each displaced geometry runs a *full new SCC*, `h = 10^{-4}` Å):

| component | analytic [eV/Å] | re-SCC FD [eV/Å] | \|d\| [eV/Å] | tol |
| --------- | --------------- | ---------------- | ------------ | --- |
| \(F_{O,x}\) | 9.648629765142 | 9.648630016201 | **2.5e-07** | 1e-5 |
| \(F_{H1,y}\) | 0.513931632212 | 0.513931627211 | **5.0e-09** | 1e-5 |

Tolerance `1e-5` eV/Å is ~40× the measured residual and is set by the
**central-difference truncation error** at `h = 1e-4` Å, not by any physical
uncertainty: a converged SCC-DFTB3 gradient is stationary, so the only
disagreement left is \(O(h^2)\).

### 2. SCC consistency — the FD residual must fall like \(h^2\)

A missing, mis-signed, or frozen-Δq-inconsistent term leaves an
**\(h\)-independent** offset against the FD; the true gradient leaves only
truncation. Measured on \(F_{O,x}\):

| \(h\) [Å] | \|analytic − FD\| [eV/Å] |
| --------- | ------------------------ |
| 8e-4 | 1.608e-05 |
| 2e-4 | 1.005e-06 |

ratio **16.0** = \(4^2\) — textbook second order. Gate:
`err(2e-4) < 0.25·err(8e-4) + 2e-9`. This is the gate that a "close enough"
loose-tolerance comparison cannot fake.

### 3. Translational invariance

\(|\sum_i F_i| \le 10^{-10}\) eV/Å per Cartesian. Every term is a pair force
with \(F_B=-F_A\), so this is exact by construction (round-off only) — it is
deliberately **four orders tighter** than the DFTB2 sibling's 1e-6, which was
sized for a lattice-image walk.

### 4. Third order must move the *force*, not just the energy

\(F_{O,x}(\mathrm{DFTB3}) - F_{O,x}(\mathrm{DFTB2}) = 3.917\times10^{-3}\)
eV/Å at this geometry (gate: > 1e-6). Without a live
\(\partial\Gamma_{AB}/\partial R\) term this difference would still be nonzero
(the densities differ) but the FD identity in §1/§2 would fail — the two gates
are not redundant.

### 5. Provisional ATV-twin determinism pins

| Quantity | Value | Unit | Tol |
| -------- | ----- | ---- | --- |
| \(E_{\mathrm{tot}}\) | −91.710362597147 | eV | 1e-6 |
| \(F_{O,x}\) | 9.648629765142 | eV/Å | 1e-6 |
| \(F_{H1,y}\) | 0.513931632212 | eV/Å | 1e-6 |

Measured 2026-07-29, Release, `OMP_NUM_THREADS=1`. Informational at the same
point: \(E^{(3)} = -1.296771\times10^{-3}\) eV.

**These are self-comparison pins, not oracle values** — same status as the
energy pin in `regressions/semiempirical/dftb3-water-mio/`. They pin
determinism of the ATV DFTB3 Hamiltonian twin; they claim **no** DFTB+
identity.

## Gaps

- **No DFTB+ binary exists on this host** (`module spider dftb` → nothing;
  `dftb+`, `mopac`, `xtb` all absent), so no absolute DFTB+ 24.1 force could be
  generated and none is claimed. The oracle here is the **formula** (Gaus 2011)
  plus the internal-consistency FD identity, which is exact and
  parameter-independent. Replace §5 with a DFTB+ 24.1 + mio-1-1 force once the
  offline pin infrastructure lands under `tools/tb_oracle/inputs/dftbplus/`.
- The in-tree SK tables are mio-**style** sketches for every pair except H–H,
  so even a DFTB+ run would not reproduce §5 without a real `.skf` reader.
  SK fidelity is a separate slug (`dftb2-3ob-organic`,
  `dftb2-dftbplus-energy-pin`).
- **Molecular only.** DFTB3 walks minimum-image pairs; there is no
  lattice-image Γ-PBC or k-mesh gradient (that is DFTB2's
  `dftb2-periodic-forces` / `dftb2-multi-k-forces`). A periodic DFTB3 force is
  not implemented and not gated.
- **No DFTB3 cell stress**, **no open shell**, **no ThreeOb tables**, **no D3**
  on the DFTB3 path. Each is an owed slug, not a silent limitation of this one.
- This regression exercises one geometry of one molecule, and water pairs
  O(s,p) only with H(s) — its gradient never touches the **pp-σ/π** SK
  derivative. That branch is covered by the unit gate
  `tests/cpu/semiempirical/test_dftb3.cpp::AnalyticForceMatchesFdCOWithPBlocksOnBothCentres`
  (off-axis C–O, all 6 components vs re-SCC FD ≤ 1e-6 eV/Å). An N-containing
  DFTB3 force case is still owed.

## Perf

`perf_exempt: true` — no DFTB+ binary on the pin host, so no fair wall ratio.
Algorithmically the analytic path replaces \(6N\) full re-SCC solves with one
post-SCC contraction (this 3-atom case: 18 solves → 0).

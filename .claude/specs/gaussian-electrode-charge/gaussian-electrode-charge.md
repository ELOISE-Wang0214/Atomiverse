---
slug: gaussian-electrode-charge
status: approved
created: 2026-07-23
track: metalwalls-method-loop
chapter: E2
grilled: false
oracle: literature
oracle_reason: >
  Gaussian electrode charges are MetalWalls specialty physics (map E2).
  No MetalWalls binary and no LAMMPS ELECTRODE package on this host;
  offline pins are analytic evaluations of the published Reed / MetalWalls
  kernel (same equations MW implements).
---

# gaussian-electrode-charge

## Summary

Ship **Gaussian charge kernels** for electrode atoms on the classical MD
surface (map slice **E2**). Electrode atoms carry normalized 3D Gaussians of
width parameter \(\eta\) [Å⁻¹] rather than point charges; the pair Green
function is \(\mathrm{erf}(\eta_{ij} r)/r\) with
\(\eta_{ij}=\eta_i\eta_j/\sqrt{\eta_i^2+\eta_j^2}\).

Public API:

| Symbol | Role |
| ------ | ---- |
| `cpu::GaussianChargeKernel` | Value-type math: \(\eta_{ij}\), \(G\), \(G_{\mathrm{self}}\), \(dG/dr\) |
| `cpu::PairGaussianCoul` | Potential: Gaussian–Gaussian / Gaussian–point / point–point Coulomb [eV, Å] |

No public `*Config` bag. No `make_*` / `build_*` factories. Not a new Driver
category. Charge solve for full constant-potential BO remains **E1**
(`const-pot-bo-parallel-plate`); this slug only supplies the **kernels** and
documents \(E\) / toy \(Q\) vs point-charge \(A\).

## Domain basis

Reed, Lanning, Madden, *J. Chem. Phys.* **126**, 084704 (2007);
MetalWalls JOSS (Marin-Laflèche et al., 2020) / Scalfi constant-potential
reviews. Units: \(E=\texttt{units::coulomb}\,q_i q_j\,G\) with
`units::coulomb = 14.399645` eV·Å/e².

\[
\rho_i(\mathbf{r})=q_i\frac{\eta_i^3}{\pi^{3/2}}
  e^{-\eta_i^2|\mathbf{r}-\mathbf{R}_i|^2},\quad
\eta_{ij}=\frac{\eta_i\eta_j}{\sqrt{\eta_i^2+\eta_j^2}},
\]

\[
G_{ij}(r)=\frac{\mathrm{erf}(\eta_{ij} r)}{r}\ (r>0),\quad
G_{ii}=\eta_i\sqrt{\frac{2}{\pi}},
\]

\[
U=\texttt{coulomb}\Big[
  \sum_{i<j}q_i q_j G_{ij}(r_{ij})
  +\tfrac12\sum_i q_i^2 G_{ii}\Big].
\]

Gaussian–point (\(\eta_j\to\infty\)): \(\eta_{ij}\to\eta_i\). Point–point:
\(G=1/r\), no self term.

**Point-charge \(A\) is an approximation** for electrodes: replacing
\(G=\mathrm{erf}(\eta r)/r\) by \(1/r\) (and dropping \(G_{ii}\)) changes both
energy and the linear system for \(Q\). Documented in REFERENCE.md.

## Architecture

- Math lives in `include/atomiverse/cpu/potentials/gaussian_charge.h`
  (header-only; pure, no Frame).
- Force/energy component: `PairGaussianCoul` under `cpu/potentials/`
  (all-pairs with min-image + cutoff, same shape as `PairYukawa` — no
  NeighborConsumer required for the small-system foothold).
- Per-type \(\eta\): constructor
  `PairGaussianCoul(cutoff, eta_by_type)` with `eta ≤ 0` ⇒ point charge.
- Optional member `include_self` (default true) for Gaussian self-energy.
- Future const-pot BO (E1) reuses `GaussianChargeKernel` for private \(A\);
  do **not** weld solve into `KSpace*`.

## Geometry (pin)

**Parallel-plate toy (4 electrode atoms + 1 ion), open / large cubic box
(no PBC image for the pin pairs — box 100 Å, cutoff 40 Å):**

| Role | type | \(\eta\) [Å⁻¹] | \(q\) [e] | position [Å] |
| ---- | ---- | ------------- | --------- | ------------ |
| L0 | 1 | 1.0 | (solved / fixed) | (0,0,0) |
| L1 | 1 | 1.0 | | (2,0,0) |
| R0 | 1 | 1.0 | | (0,0,10) |
| R1 | 1 | 1.0 | | (2,0,10) |
| ion | 0 | point | +1 | (1,0,5) |

**Fixed-charge energy pin** (literature kernel): electrode charges
\(q_L=+0.25\) on L0,L1; \(q_R=-0.25\) on R0,R1; ion \(+1\).

**Toy \(Q\) pin** (2-site plates, no ion): two Gaussian sites at \(z=0\) and
\(z=L=10\) Å, \(\eta=1.0\) Å⁻¹, target potentials \(\Psi_L=+\Delta/2\),
\(\Psi_R=-\Delta/2\) with \(\Delta=1\) V-equivalent (in ATV units the
right-hand side is \(\Psi\) such that \(A Q = \Psi\) with
\(A_{ij}=\texttt{coulomb}\,G_{ij}\)). Neutrality \(Q_L+Q_R=0\).
Compare to point-charge \(A\) (\(G=1/r\), no self) for the same \(\Psi\).

## Offline pin

- **oracle:** literature (Reed 2007 + MetalWalls kernel; no live MW/LAMMPS
  ELECTRODE on agent host)
- **oracle_version:** Reed JCP 126 084704 (2007) eqs; MetalWalls JOSS 2020
- Analytic values recorded in
  `regressions/md/gaussian-electrode-charge/REFERENCE.md`
- Dual-gate: energy tol \(1\times10^{-9}\) eV; \(Q\) tol \(1\times10^{-9}\) e;
  Gaussian vs point-charge \(Q\) difference documented (must be nonzero)

## Tasks

- [x] Spec + acceptance (this file + `.acceptance.md`)
- [x] Offline pin → REFERENCE.md + optional oracle input notes
- [x] `GaussianChargeKernel` (η_eff, greens, self_greens, dgreens_dr)
- [x] `PairGaussianCoul` Potential (type→η, self energy, Frame I/O)
- [x] Unit tests under `tests/cpu/potentials/`
- [x] E2e regression `regressions/md/gaussian-electrode-charge/`
- [x] CMake: unit test + `atv_add_regression(md gaussian-electrode-charge)` under `ATV_BUILD_MD`
- [x] MD-only build green; unit + regression PASS
- [x] Perf note: perf_exempt (no MW/LAMMPS ELECTRODE twin binary)

## Testing strategy

- **Unit:** analytic \(G\), \(G_{\mathrm{self}}\), dimer energy/force vs closed form;
  FD force check; point limit \(\eta\to\infty\) recovers \(1/r\).
- **Regression:** 5-atom fixed-\(Q\) energy vs literature pin; 2-site toy
  BO \(Q\) Gaussian vs point-charge \(A\).

## Out of scope

- Full multi-electrode BO / CG solve (E1 / later E2)
- 2D Ewald / slab correction (sibling slug)
- Thomas–Fermi, MaZe, PIM coupling
- CUDA electrode path
- Absorbing kernels into `KSpaceEwald` / `KSpacePME`

---

**Librarian note:** map E2; Component on MD/Evaluate; confidence high (closed-form kernels).

# Reference — Thomas–Fermi finite metallicity electrode (gold capacitor toy)

oracle: metalwalls  
oracle_version: Scalfi et al. JCP 153, 174704 (2020) model as in MetalWalls (TF
hardness \(A_{ii}\mathrel{+}= l_{\mathrm{TF}}^{2} d/\epsilon_0\)); offline discrete
pin via `tools/metalwalls_oracle/inputs/metalwalls/thomas-fermi-electrode/pin_gold_cap.py`
(no live MetalWalls binary on the agent host)

## Science

Thomas–Fermi electrode metallicity (map **E7**): Scalfi, Dufils, Reeves,
Rotenberg, Salanne, *J. Chem. Phys.* **153**, 174704 (2020),
DOI [10.1063/5.0028232](https://doi.org/10.1063/5.0028232) (arXiv:1910.13341).
**Not** available as a LAMMPS ELECTRODE feature (perfect conductor only) →
oracle priority binds to **MetalWalls**.

Diagonal hardness (metal units \(e\), Å, eV; \(k_e=\) `units::coulomb`):

\[
A_{ii}\;\mathrel{+}= 4\pi\,k_e\,l_{\mathrm{TF}}^{2}\,d.
\]

Perfect metal is \(l_{\mathrm{TF}}=0\).

## System (gold capacitor toy)

Empty parallel-plate capacitor (no electrolyte):

| Quantity | Value |
| --- | --- |
| Au lattice \(a_0\) | 4.07 Å |
| Site density \(d = 4/a_0^3\) | **0.059330329963** Å⁻³ |
| Plane spacing \(a = a_0/2\) | 2.035 Å |
| Planes per electrode | 3 (1 site per plane at cell centre) |
| \(L_x = L_y\) | \(a_0 = 4.07\) Å |
| First-plane gap \(L\) | 20.0 Å |
| \(\Psi_{\mathrm{bot}},\Psi_{\mathrm{top}}\) | 0 V, 1 V (\(\Delta\Psi=1\) V) |
| Gaussian \(\eta = 1/0.56\) | **1.785714285714** Å⁻¹ (MetalWalls \(w=0.56\) Å) |
| 2D \(xy\) image half-range | 40 |
| Electroneutrality | on |

Atom \(z\) [Å]: bot `{5.0, 7.035, 9.07}`, top `{29.07, 31.105, 33.14}`
(types 1 / 2). Box \(L_z\) large enough that the 2D-image electrode path does
not wrap electrodes through \(z\).

## Offline pin (discrete ConstPotBO algebra = MW TF model)

Script: `tools/metalwalls_oracle/inputs/metalwalls/thomas-fermi-electrode/pin_gold_cap.py`

| \(l_{\mathrm{TF}}\) [Å] | \(Q_{\mathrm{top}}\) [e] | \(C=\|Q\|/\Delta\Psi\) [e/V] |
| --- | --- | --- |
| **0.0** (perfect metal) | **+4.620945546104e−03** | **4.620945546104e−03** |
| **2.0** (finite metallicity) | **+3.934601758679e−03** | **3.934601758679e−03** |

Continuum cross-check (paper Eq. 8 with \(L_{\mathrm{vac}}=L-a\)):
\(C_{\mathrm{cont}}(0)\approx 5.10\times 10^{-3}\), \(C_{\mathrm{cont}}(2)\approx 4.17\times 10^{-3}\)
[e/V] on this area — same trend (discrete 1-site/plane lattice under-samples
lateral density vs continuum; dual gate is the **discrete** pin above).

## What Atomiverse asserts

Public API: `cpu::md::ConstPotBO` members `l_tf`, `site_density`,
`capacitance()`; still a Component on `md::Phase::pre_force`.

| Check | Tol |
| --- | --- |
| \(\lvert\sum Q_i\rvert\) | \(\le 10^{-10}\) e |
| \(Q_{\mathrm{top}}(l_{\mathrm{TF}}=0)\) vs pin | **≤ 1%** rel |
| \(Q_{\mathrm{top}}(l_{\mathrm{TF}}=2)\) vs pin | **≤ 1%** rel |
| \(C(l_{\mathrm{TF}}=2) < C(l_{\mathrm{TF}}=0)\) | exact |
| \(C\) dual-gate vs pin at both \(l_{\mathrm{TF}}\) | **≤ 1%** rel |
| \(\mathrm{sign}(Q_{\mathrm{top}})=+1\) at \(\Delta\Psi=+1\) V | exact |

## Perf

```
perf_exempt: true
perf_exempt_reason: no MetalWalls binary on agent host for fair TF twin wall; algebraic BO+TF solve is sub-second
ratio_atv_over_oracle: n/a
within_10x: exempt
```

## Repro of pin (offline only — never ctest)

```bash
python3 tools/metalwalls_oracle/inputs/metalwalls/thomas-fermi-electrode/pin_gold_cap.py
# optional: MetalWalls empty-capacitor TF input if mw is installed
```

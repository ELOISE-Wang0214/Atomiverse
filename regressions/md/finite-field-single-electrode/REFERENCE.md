# Reference — finite-field + single electrode (3D, Dufils-style)

oracle: metalwalls  
oracle_version: literature analytic — Dufils et al., Phys. Rev. Lett. 123, 195501 (2019) finite-field + single constant-potential electrode; Reed/Siepmann–Sprik Gaussian electrode kernels (same as MetalWalls JOSS 5 2373 (2020) / Atomiverse ConstPotBO)

## Why not LAMMPS

LAMMPS ELECTRODE `ffield on` implements the **two-electrode** finite-field path
(\(E_z=\Delta\Psi/L_z\), requires exactly two groups). This slug is the
**single-wall** Dufils product (MetalWalls-only geometry). Binding oracle rule
→ metalwalls / literature analytic when no MW binary is available.

## Physical system (vacuum single-site pin)

Cubic orthorhombic box, 3D PBC, **no electrolyte**, one Gaussian electrode site.

| Quantity | Value |
| --- | --- |
| Box | \(L = 20\) Å cubic |
| \(\eta\) | \(2.0\) Å⁻¹ |
| \(\Psi\) | \(0\) V |
| \(\mathbf{E}\) | \((0,\,0,\,0.1)\) V/Å |
| Electrode | type 1 at \(\mathbf{r}=(0,0,10)\) Å |
| Electroneutrality | **off** (free total \(Q\)) |
| \(xy\) images | 0 (primary cell) |

## Offline analytic pin

Elastance self-term (Reed Gaussian; metal units `units::coulomb = 14.399645`):

\[
A_{ii} = k_e\,\eta\sqrt{2/\pi} = 22.978508853084353\ \mathrm{V}/e
\]

BO with external field \(\Phi_{\mathrm{ext}}=-\mathbf{E}\cdot\mathbf{r}\):

\[
A_{ii} Q = \Psi + E_z z
\quad\Rightarrow\quad
Q_{\mathrm{pin}} = \frac{0.1\times 10}{A_{ii}} = \mathbf{0.04351892485250505}\ e
\]

Polarization density:

\[
P_z = \frac{Q z}{V},\quad V=L^3=8000\ \mathrm{Å}^3
\quad\Rightarrow\quad
P_{z,\mathrm{pin}} = \mathbf{5.439865606563132\times 10^{-5}}\ e/\mathrm{Å}^2
\]

Field force on a test ion \(q=+1\) at any position (via `FixEfield`):

\[
F_z = q E_z = 0.1\ \mathrm{eV/Å}
\]

## What Atomiverse asserts

`reg_finite-field-single-electrode`:

1. `ConstPotBO` with one electrode, `Ez=0.1`, `electroneutrality=false` →
   \(|Q - Q_{\mathrm{pin}}| \le 10^{-10}\) e.
2. `ClassicalPolarization` after charge write →
   \(|P_z - P_{z,\mathrm{pin}}| \le 10^{-12}\) e/Å².
3. Dual-gate force: `FixEfield{0,0,0.1}` on a \(q=+1\) ion →
   \(|F_z - 0.1| \le 10^{-12}\) eV/Å.
4. MD binding path (one VV step) reproduces \(Q\) within the same tol.

## Perf

```
perf_exempt: true
perf_exempt_reason: no MetalWalls binary on agent host; single-electrode finite-field is MW-only (LAMMPS ffield requires two electrodes) — no fair twin wall
ratio_atv_over_oracle: n/a
within_10x: exempt
```

## Repro of pin (offline only — never ctest)

```python
import math
coulomb = 14.399645
eta, Ez, z, L = 2.0, 0.1, 10.0, 20.0
Aii = coulomb * eta * math.sqrt(2 / math.pi)
Q = (Ez * z) / Aii
Pz = Q * z / (L ** 3)
print(Aii, Q, Pz)
```

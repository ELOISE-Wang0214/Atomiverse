# Reference — coupled PIM dipoles ⊕ electrode \(Q\) (E5)

oracle: metalwalls  
oracle_version: MetalWalls multipole + Gaussian electrode BO form (Coretti et al. JCP 2022 polarizable interfaces; Marin-Laflèche et al. JOSS 5 2373 (2020)); offline pin script `tools/metalwalls_oracle/inputs/metalwalls/dipole-electrode-joint/pin_joint.py` (no MetalWalls binary on agent host)  
perf_exempt: true  
perf_exempt_reason: No MetalWalls binary on PATH for a fair joint Q–μ wall twin; pin is the published MW multipole ⊕ Gaussian-electrode algebra (same path as thomas-fermi / bulk-PIM offline pins).

## Physics

Self-consistent joint of:

1. **Constant-potential BO** electrode charges \(Q\) with Gaussian elastance \(A\)
   (same kernels as `cpu::md::ConstPotBO`).
2. **PIM induced dipoles** \(\boldsymbol{\mu}\) on polarizable electrolyte sites
   with Thole damping (same kernels as `cpu::PimInducedDipole`).

Coupling:

\[
A Q = \Psi - \Phi^{\mathrm{perm}} - \Phi^{\mu},\quad
\Phi_i^{\mu} = k_e\sum_j\frac{\boldsymbol{\mu}_j\cdot\mathbf{r}_{ij}}{r_{ij}^3}
\]

\[
\boldsymbol{\mu}_i = \alpha_i\mathbf{E}_i^{\mathrm{tot}}(q_{\mathrm{all}},\boldsymbol{\mu})
\]

Outer alternating solve until \(\max|\Delta Q|\) and \(\max|\Delta\mu|\) below
`joint_tol`. **Oracle is MetalWalls** — LAMMPS has no SCF PIM ⊕ const-pot
joint twin (ELECTRODE is non-polarizable; Drude is a different model).

## Geometry — offline pin

Orthorhombic \(L_x=L_y=10\) Å, \(L_z=20\) Å; \(\eta=2\) Å⁻¹; \(\Psi_{\mathrm{bot}}=0\),
\(\Psi_{\mathrm{top}}=1\) V; electroneutrality on; `xy_image_range=0` (min-image);
`thole_a=0.39`.

| Site | type | \(q\) [e] | \(\alpha\) [Å³] | \((x,y,z)\) [Å] |
| ---- | ---- | --------- | --------------- | --------------- |
| bot graphene wall | 1 | solved | 0 | (0, 0, 2) |
| top graphene wall | 2 | solved | 0 | (0, 0, 18) |
| IL cation | 3 | +1 | 1.0 | (0, 0, 8) |
| IL anion | 4 | −1 | 1.5 | (0, 0, 12) |

### Pin values (hard-coded in regression)

| Observable | Value |
| ---------- | ----- |
| \(Q_{\mathrm{bot}}\) [e] | **−7.3285008958883610e−02** |
| \(Q_{\mathrm{top}}\) [e] | **+7.3285008958883610e−02** |
| \(\mu_{\mathrm{cation},z}\) [e·Å] | **+5.8951916148990660e−02** |
| \(\mu_{\mathrm{anion},z}\) [e·Å] | **+8.7432701391015569e−02** |
| \(\mu_{*,x}\), \(\mu_{*,y}\) | 0 |
| \(\lvert\sum Q\rvert\) | 0 (within \(10^{-12}\)) |
| \(\alpha=0\) \(Q_{\mathrm{top}}\) [e] | **+7.5339627935767609e−02** |

**Tol (dual gate):**  
relative \(|Q-Q_{\mathrm{pin}}|/|Q_{\mathrm{pin}}| \le 10^{-6}\);  
\(|\mu_z-\mu_{\mathrm{pin}}| \le 10^{-8}\) e·Å; transverse \(|\mu|<10^{-10}\);  
neutrality \(|\sum Q|\le 10^{-10}\) e.

**API:** `cpu::md::DipoleElectrodeJoint{2.0, {0,0,0,1.0,1.5}}` +
`add_electrode(1,0).add_electrode(2,1)`.

## Provenance

- Coretti et al., *J. Chem. Phys.* (2022) — polarizable electrolyte + metallic
  electrodes (MetalWalls `dipoles_and_electrodes`).  
- Marin-Laflèche et al., JOSS **5**, 2373 (2020).  
- Thole, *Chem. Phys.* **59**, 341 (1981).  
- Offline script: `tools/metalwalls_oracle/inputs/metalwalls/dipole-electrode-joint/pin_joint.py`.  
- Map: MetalWalls capability **E5** (`dipole-electrode-joint`).

## Fair wall (perf gate)

No MetalWalls binary available on the agent host → **perf_exempt**.  
When a MW twin appears, record \(T_{\mathrm{atv}}/T_{\mathrm{MW}}\) here and
require \(\le 10\).

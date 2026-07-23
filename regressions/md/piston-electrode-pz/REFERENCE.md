# Reference — piston electrode constant \(P_z\) (E6)

```text
oracle: metalwalls
oracle_version: algorithm MetalWalls rigid-piston constant-P_z
  (Marin-Laflèche et al., JOSS 5, 2373 (2020), DOI 10.21105/joss.02373;
   capability map E6). No MetalWalls binary available offline on worker;
   pin is the MW force-identity F_ext = ±P A on a rigid electrode + mean
   <P_z> → P_target equilibration (self-consistent probe → target).
```

## Map slice

**E6 — electrode mechanics:** piston electrode + constant normal pressure
\(P_z\). Rigid piston ≠ generic isotropic barostat. Orthogonal to E1 charge
solve (pair forces only in this foothold).

## Algorithm pin (offline, not a live binary)

MetalWalls treats electrodes as **rigid pistons** for constant-\(P_z\)
runs (JOSS summary; map §1 / §5). Atomiverse `cpu::md::PistonElectrode`:

\[
F_z^{\mathrm{atoms}}=\sum_{i\in t} F_{z,i},\quad
A=L_x L_y,\quad
F_{\mathrm{ext}}=\begin{cases}
-P_{\mathrm{target}} A & \text{Upper}\\
+P_{\mathrm{target}} A & \text{Lower}
\end{cases}
\]

\[
F_{\mathrm{net}}=F_z^{\mathrm{atoms}}+F_{\mathrm{ext}}(+F_{\mathrm{damp}}),\quad
F_{z,i}\leftarrow F_{\mathrm{net}}\,m_i/M,\quad F_{x,i}=F_{y,i}=0.
\]

Instantaneous diagnostic (Upper): \(P_z=F_z^{\mathrm{atoms}}/A\).
Equilibrium: \(\langle P_z\rangle\to P_{\mathrm{target}}\).

## Geometry (regression)

| Species | type | role |
| --- | ---: | --- |
| Ar-like fluid | 1 | mobile LJ fluid (16 atoms, two layers) |
| lower wall | 2 | frozen (`FixFreeze`) at \(z=2\) Å |
| upper wall | 3 | **piston** (`PistonElectrode::Upper`) at \(z=8\) Å |

- Box: orthorhombic \(L_x=L_y=10\) Å, \(L_z=20\) Å → \(A=100\) Å²
- Walls: \(2\times2\) atoms each
- Fluid layers at \(z_{\mathrm{lo}}+2.2\) and \(z_{\mathrm{hi}}-2.2\) (repulsive contact)
- Pair: LJ \(\varepsilon=0.10\) eV, \(\sigma=3.0\) Å, \(r_c=6\) Å
- NVT fluid: VelocityVerlet + NoseHoover \(T=300\) K, \(\tau_T=50\) fs, \(dt=0.5\) fs
- Piston: \(P_{\mathrm{target}}=\) probe \(P_z>0\) at compressed geometry (self-consistent offline pin); \(\tau_{\mathrm{damp}}=20\) fs
- Warm 400 + prod 600 steps (chunk-sampled \(P_z\))

## Observables and tolerances

| Check | Bound | Rationale |
| --- | --- | --- |
| Probe \(P_z\) | \(>10^{-8}\) eV/Å³ | Compressed fluid pushes upper wall |
| \(\lvert\langle P_z\rangle-P_{\mathrm{target}}\rvert/P_{\mathrm{target}}\) | **≤ 0.10** | MW constant-\(P_z\) equilibration (10% textbook gate) |
| Upper wall coplanar | \(\Delta z\le 10^{-9}\) Å | Rigid piston |
| Upper \(\lvert v_x\rvert,\lvert v_y\rvert\) | ≤ \(10^{-12}\) Å/fs | \(z\)-only motion |
| Lower wall \(\Delta z\) | ≤ \(10^{-10}\) Å | `FixFreeze` |
| Energy finite | true | Smoke stability |

Unit tests pin the static identities \(P_z=F_z^{\mathrm{atoms}}/A\) and
\(F_{\mathrm{net}}=0\) when \(P_{\mathrm{target}}=P_z^{\mathrm{atoms}}\) to
machine precision (no MD statistics).

## Perf

```text
perf_exempt: true
perf_exempt_reason: no MetalWalls binary on the offline worker PATH;
  cannot fair-bench vs MW twin for this MW-only rigid-piston path.
  (LAMMPS fix wall/piston is shock-wave kinematics, not electrode P_z.)
```

## Units

Public surface: eV, Å, fs, amu, eV/Å, eV/Å³ (pressure). 1 bar ≈ 6.2415×10⁻⁷ eV/Å³.

## API

```cpp
cpu::md::PistonElectrode piston(/*type=*/3, /*P [eV/Å³]=*/P,
                                cpu::md::PistonElectrode::Side::Upper);
piston.tau_damp = 20.0; // optional [fs]
sim.add(std::move(piston)); // or add_ref for diagnostics
sim.add(cpu::md::FixFreeze{/*lower type=*/2});
```

No public `*Config` bag; no `make_*` factory.

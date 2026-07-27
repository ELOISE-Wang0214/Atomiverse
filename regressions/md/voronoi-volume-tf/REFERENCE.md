# Reference — Voronoi / site-volume companion for Thomas–Fermi metallicity

oracle: metalwalls  
oracle_version: MetalWalls TF model per Scalfi et al. JCP 153, 174704 (2020)
(per-species electrode inputs `thomas_fermi_length` + `voronoi_volume`,
gitlab `ampere2/metalwalls` branch `release`); offline discrete pin via
`tools/metalwalls_oracle/inputs/metalwalls/voronoi-volume-tf/pin_voronoi_tf.py`
(no MetalWalls binary on the agent host)

## Science

Map **§1 "Voronoi volume (TF)"** — the companion of the wave-2 E7 foothold
`thomas-fermi-electrode`. Scalfi, Dufils, Reeves, Rotenberg, Salanne,
*J. Chem. Phys.* **153**, 174704 (2020),
DOI [10.1063/5.0028232](https://doi.org/10.1063/5.0028232): the TF
metallicity self-energy is per-site,

\[ U_{\mathrm{TF}} = \sum_i \frac{l_{\mathrm{TF}}^2}{2\,\epsilon_0\,v_i}\,q_i^2
   \;\Rightarrow\;
   A_{ii} \mathrel{+}= 4\pi k_e\,\frac{l_{\mathrm{TF}}^2}{v_i}, \]

with \(v_i\) the **Voronoi volume** of the electrode site — MetalWalls'
per-species `voronoi_volume` input. For FCC gold \(v = a_0^3/4\) exactly
(4 sites tile the conventional cube; the FCC Voronoi cell is the rhombic
dodecahedron). LAMMPS ELECTRODE has no TF/Voronoi twin (perfect conductor
only) → oracle priority binds to **MetalWalls**.

## System

1. **Bulk FCC Au**, 2×2×2 conventional cells (32 sites), \(a_0 = 4.07\) Å,
   cubic box \(2a_0\), fully periodic — input to `cpu::VoronoiVolume`.
2. **Gold capacitor toy** — identical to `regressions/md/thomas-fermi-electrode`:
   3+3 electrode planes (1 site/plane), \(L_x=L_y=a_0\), plane spacing
   \(a_0/2\), first-plane gap 20 Å, \(\Delta\Psi = 1\) V, Gaussian
   \(\eta = 1/0.56\) Å⁻¹, 2D xy image half-range 40, electroneutrality on,
   \(l_{\mathrm{TF}} = 2\) Å.

## Offline pins (hard-coded; discrete ConstPotBO algebra = MW TF model)

Script: `tools/metalwalls_oracle/inputs/metalwalls/voronoi-volume-tf/pin_voronoi_tf.py`
(numpy 2.5.1, run offline 2026-07-27). Wave-2 perfect-metal pin from
`thomas-fermi-electrode`.

| Quantity | Pin | Units |
| --- | --- | --- |
| FCC Voronoi volume \(v = a_0^3/4\) | **16.854785750000** | Å³ |
| \(Q_{\mathrm{top}}\), \(l_{\mathrm{TF}}=2\) Å, \(v_i = a_0^3/4\) | **+3.934601758679e−03** | e |
| \(Q_{\mathrm{top}}\), \(l_{\mathrm{TF}}=2\) Å, \(v_i = a_0^3/8\) | **+3.509853715789e−03** | e |
| \(Q_{\mathrm{top}}\), \(l_{\mathrm{TF}}=0\) (perfect metal) | **+4.620945546104e−03** | e |

\(Q_{\mathrm{top}}(v_i = a_0^3/4)\) equals the wave-2 fixed-density pin
(\(d = 4/a_0^3\)) to all digits — the volume and density forms are the same
algebra, asserted below as an equivalence invariant.

## What Atomiverse asserts

Public API: `cpu::VoronoiVolume::volumes(frame)` (periodic tessellation
volumes [Å³]) and `cpu::md::ConstPotBO::site_volume` (per-row \(v_i\) [Å³],
TF hardness \(4\pi k_e l_{\mathrm{TF}}^2/v_i\)); ConstPotBO stays a
`pre_force` Component.

| Check | Tol |
| --- | --- |
| every bulk-FCC \(v_i\) vs \(a_0^3/4\) | ≤ 1e−9 rel |
| \(\sum v_i\) vs box volume \(8a_0^3\) | ≤ 1e−10 rel |
| \(\lvert\sum Q_i\rvert\) (each solve) | ≤ 1e−10 e |
| \(Q_{\mathrm{top}}(v_i)\) vs MW-model pin | ≤ 1 % rel |
| \(Q_{\mathrm{top}}(v_i)\) vs fixed-density twin (in-binary) | ≤ 1e−10 rel |
| \(C(v_i)\) vs pin | ≤ 1 % rel |
| \(Q_{\mathrm{top}}(v_i/2)\) vs pin | ≤ 1 % rel |
| \(C(v/2) < C(v) < C(l_{\mathrm{TF}}{=}0)\) | exact ordering |
| \(Q_{\mathrm{top}}(l_{\mathrm{TF}}{=}0)\) vs wave-2 pin | ≤ 1 % rel |

## Perf

```
perf_exempt: true
perf_exempt_reason: no MetalWalls binary on agent host for a fair TF/Voronoi wall twin; LAMMPS has no TF electrode feature; run is a sub-second algebraic solve + 32-site tessellation
ratio_atv_over_oracle: n/a
within_10x: exempt
```

## Repro of pin (offline only — never ctest)

```bash
python3 tools/metalwalls_oracle/inputs/metalwalls/voronoi-volume-tf/pin_voronoi_tf.py
```

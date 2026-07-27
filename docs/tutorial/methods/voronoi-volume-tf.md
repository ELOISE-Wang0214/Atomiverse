# Voronoi site-volume companion for Thomas–Fermi metallicity

*(Technical-report draft; narrative house style TBD — rewrite later.)*

| Field | Value |
| ----- | ----- |
| Slug | `voronoi-volume-tf` |
| Regression | `regressions/md/voronoi-volume-tf/` |
| Public API | `cpu::VoronoiVolume::volumes(frame)` [Å³]; `cpu::md::ConstPotBO::site_volume` (per-row \(v_i\) [Å³]) |
| Roadmap / queue | `docs/developer/metalwalls-method-queue.md` → map §1 "Voronoi volume (TF)" (wave 5); E7 companion |
| Oracle | **MetalWalls** TF model / Scalfi et al. 2020 (offline algebraic pin; LAMMPS has no TF twin) |

## Summary

The wave-2 E7 foothold (`thomas-fermi-electrode`) hard-codes one global site
density \(d\) in the Thomas–Fermi hardness. MetalWalls instead takes the
**Voronoi volume of the electrode site** as a per-species input
(`voronoi_volume` next to `thomas_fermi_length`). This slug lands the
companion: `cpu::VoronoiVolume` computes periodic Voronoi cell volumes from
geometry, and `cpu::md::ConstPotBO` accepts per-site volumes via the public
member `site_volume`, applying

\[
A_{ii}\;\mathrel{+}= 4\pi\,k_e\,\frac{l_{\mathrm{TF}}^{2}}{v_i}.
\]

`ConstPotBO` stays a Component on `cpu::md::MD` at `pre_force`; no new Driver
category, no Config bag, no charge solve inside `KSpace*`.

## Method

**Physics.** Scalfi, Dufils, Reeves, Rotenberg, Salanne, *J. Chem. Phys.*
**153**, 174704 (2020), DOI
[10.1063/5.0028232](https://doi.org/10.1063/5.0028232): the semiclassical TF
self-energy is per-site, \(U_{\mathrm{TF}} = \sum_i
l_{\mathrm{TF}}^2 q_i^2 / (2\epsilon_0 v_i)\), with \(v_i\) the Voronoi
volume of atom \(i\) in the electrode crystal. The fixed-density foothold is
the special case \(v_i = 1/d\); for FCC gold \(v = a_0^3/4\) exactly.

**Tessellation.** `cpu::VoronoiVolume` computes cells by convex bisector
clipping: the initial cell is the box cuboid centred on the site (the
orthorhombic self-image bisectors), clipped by perpendicular-bisector
half-spaces of all neighbor images that can reach it, nearest first with
early termination. \(\sum_i v_i = V_{\mathrm{box}}\) by construction of a
complete tessellation. Orthorhombic boxes only; throws on triclinic — a
setup-time companion, not a hot-path Component (MW likewise treats
`voronoi_volume` as a static input).

## Usage

```cpp
// 1. Site volumes from the bulk electrode crystal (MetalWalls workflow).
cpu::VoronoiVolume vv;
std::vector<double> v = vv.volumes(bulk_gold);   // [Å³]; each = a0³/4 on FCC

// 2. Feed the TF hardness per site (indexed by the run frame's rows).
cpu::md::ConstPotBO conp(/*eta=*/1.0 / 0.56);
conp.add_electrode(/*type=*/1, /*Psi=*/0.0);
conp.add_electrode(/*type=*/2, /*Psi=*/1.0);
conp.l_tf = 2.0;                                  // [Å]
conp.site_volume.assign(frame_rows, v[0]);        // per-row v_i [Å³]
sim.on(md::Phase::pre_force).add(std::move(conp));
```

Exactly one of `site_density` / `site_volume` may be set when `l_tf > 0`;
non-electrode rows of `site_volume` are ignored.

## Regression gates (`reg_voronoi-volume-tf`)

- Bulk FCC Au 2×2×2 (32 sites): every \(v_i = a_0^3/4 = 16.85478575\) Å³
  (≤ 1e−9 rel), \(\sum v_i = 8a_0^3\) (≤ 1e−10 rel).
- Gold capacitor toy, \(l_{\mathrm{TF}}=2\) Å: \(Q_{\mathrm{top}}\) with the
  computed \(v_i\) matches the MetalWalls-model pin ≤ 1 % **and** the
  fixed-density twin ≤ 1e−10 rel; electroneutrality ≤ 1e−10 e.
- Metallicity ordering \(C(v/2) < C(v) < C(l_{\mathrm{TF}}{=}0)\); halved
  volume and perfect-metal pins ≤ 1 %.
- Perf: `perf_exempt` (no MetalWalls binary; LAMMPS has no TF feature).

Pins: `tools/metalwalls_oracle/pin_notes_voronoi-volume-tf.md`.

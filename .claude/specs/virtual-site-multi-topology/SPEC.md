# virtual-site-multi-topology

**Branch:** `wf/virtual-site-multi-topology`  
**Chapter:** E9 residual (molecular water / virtual sites beyond TIP4P M-site pair)  
**Oracle:** `lammps` (mandatory — same physics as `lj/cut/tip4p/cut`; general sites reproduce that path)  
**Status:** open  
**Map:** MetalWalls E9; queue row `virtual-site-multi-topology`  
**Dep:** `four-site-tip4p-foothold` ✓ (specialized pair); this lands the **general** multi-topology table

## Goal

Ship a **general multi-topology virtual-site** Component: massless sites are
**Frame atoms** reconstructed each step from parent positions, with forces on the
site redistributed to parents (Feenstra / GROMACS-style). At least two
constructions:

1. **`Weighted3`** — TIP4P M-site: \(\mathbf{r}_s=(1-\alpha)\mathbf{r}_0+\frac\alpha2\mathbf{r}_1+\frac\alpha2\mathbf{r}_2\)
2. **`OutOfPlane3`** — TIP5P lone-pair form: \(\mathbf{r}_s=\mathbf{r}_0+a\mathbf{d}_1+b\mathbf{d}_2+c\,(\mathbf{d}_1\times\mathbf{d}_2)\)

Standard pair styles (`PairLJCut` + `PairCoulCut`) evaluate on Frame sites; the
virtual-site Component only **places** and **partitions** — it is **not** a
specialized `tip4p` pair. Frame I/O only; no public `*Config` bags; no `make_*`
/ `build_*` factories. Units: eV, Å on the public surface.

## Non-goals

- Replacing / deleting `PairLJCutTip4pCut` (stays as the specialized twin)
- Long-range TIP4P / PPPM tip4p
- Bonded flexible water calibration
- Electrode + water coupling
- CUDA path; SCF / PW / GPW
- Full GROMACS `virtual_sitesn` catalogue (COM-of-N, 2-site, …) — only Weighted3 + OutOfPlane3

## Public API (binding)

```cpp
// include/atomiverse/cpu/dynamics/virtual_sites.h
cpu::md::VirtualSites vs;

// TIP4P M: alpha = qdist / (cos(θ/2) * r_OH)
const double alpha = cpu::md::VirtualSites::tip4p_alpha(0.15, 0.9572, theta);
vs.add_weighted3(/*site=*/tag_M, /*O=*/tag_O, /*H1=*/tag_H1, /*H2=*/tag_H2, alpha);

// TIP5P L sites (Mahoney–Jorgensen geometry coeffs)
double a, b, c;
cpu::md::VirtualSites::tip5p_out_of_plane_coeffs(
    /*r_ol=*/0.70, /*theta_lol=*/109.47 * deg, /*r_oh=*/0.9572,
    /*theta_hoh=*/104.52 * deg, /*sign=*/+1, a, b, c);
vs.add_out_of_plane3(tag_L1, tag_O, tag_H1, tag_H2, a, b, c);
// sign=-1 → L2

// Evaluate single-point: place before pairs, redistribute after
cpu::Evaluate eval;
eval.on(evaluate::Phase::prepare).add_ref(vs, /*order=*/+100);   // place
eval.add_potential(cpu::PairLJCut{...});
eval.add_potential(cpu::PairCoulCut{cut}.set_special_bonds(0,0,0));
eval.on(evaluate::Phase::finalize).add_ref(vs, /*order=*/-100);  // redistribute
```

- `name()` → `"virtual_sites"`.
- Default MD bindings: place at `md::Phase::post_integrate` (+100) and
  `md::Phase::pre_force` (+50); redistribute at `md::Phase::post_force` (+50).
- Parameters = method args / static helpers; **no** `VirtualSiteConfig` bag.
- Charge lives on the site atom; parents carry topology for specials.

## Oracle

| Field | Value |
| ----- | ----- |
| **Primary** | **lammps** `pair_style lj/cut/tip4p/cut` (TIP4P dimer energy + parent forces after partition) |
| **Offline path** | Reuse / twin under `tools/metalwalls_oracle/inputs/lammps/virtual-site-multi-topology/` (or shared four-site pin) |
| **TIP5P** | Geometry dual-gate (r_OL, θ_LOL) + force-partition conservation (ΣF parents = F_site; torque identity); no LAMMPS tip5p pair on host → not a second energy pin |
| **REFERENCE** | `oracle: lammps` + `oracle_version:` |
| **ctest** | No live LAMMPS |

## Geometry / observables

**Gate A — TIP4P dimer via general sites (LAMMPS pin):** same geometry as
`four-site-tip4p-foothold` (O–O = 2.8 Å, classic Jorgensen TIP4P). Frame has
**4 atoms/water** (O,H1,H2,M). Observables: total \(E\) [eV], forces on O/H
after redistribute [eV/Å], OM = qdist.

**Gate B — TIP5P construction:** one rigid water; dual L sites; check
\(r_{OL}=0.70\) Å and \(\theta_{LOL}=109.47°\) within tol; apply unit force on
each L and verify force sum / Feenstra-style partition identities.

## Tasks

1. [x] Write this SPEC + acceptance.md (first).
2. [x] Implement `cpu::md::VirtualSites` under `cpu/dynamics/` (header + cpp; CMake list).
3. [x] Unit tests `tests/cpu/dynamics/test_virtual_sites.cpp` (Weighted3 place + partition, tip4p_alpha, OutOfPlane3 tip5p geometry, force sum).
4. [x] Offline pin + `REFERENCE.md` + e2e `regressions/md/virtual-site-multi-topology/`.
5. [x] Wire regression + unit test CMake.
6. [x] MD-only configure/build; unit + regression PASS.
7. [x] Perf gate vs LAMMPS tip4p twin — ratio **0.0837** ≤ 10.
8. [x] Commit on `wf/virtual-site-multi-topology` only.

## Acceptance

| ID | Criterion | Tol |
| -- | --------- | --- |
| A1 | Oracle is LAMMPS tip4p path; pin documented offline | — |
| A2 | Offline twin under `tools/metalwalls_oracle/inputs/lammps/` | — |
| A3 | `REFERENCE.md` has `oracle: lammps` and `oracle_version:` | — |
| A4 | TIP4P dimer \(E\) via VirtualSites+LJ+Coul matches pin within **1e-8 eV** | 1e-8 eV |
| A5 | Dual-gate geometry: OM = qdist within **1e-10 Å** after place | 1e-10 Å |
| A6 | O/H forces match LAMMPS pin within **1e-6 eV/Å**; M force zero after redistribute | 1e-6 eV/Å |
| A7 | TIP5P: \(r_{OL}\) and \(\theta_{LOL}\) within **1e-8** relative/absolute as stated in REFERENCE | dual-gate |
| A8 | Force partition: \(\sum F_{\mathrm{parents}}=F_{\mathrm{site}}\) (pre-zero) within **1e-12** | 1e-12 |
| A9 | MD-only cmake: HF/PW/GPW/CUDA/Torch **OFF** | — |
| A10 | Unit tests under `tests/` only; one e2e regression | — |
| A11 | Perf: \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\le 10\) or `perf_exempt` | ≤10× |

## Build (binding)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release \
  -DATV_BUILD_MD=ON -DATV_BUILD_KSPACE=ON \
  -DATV_BUILD_TESTS=ON -DATV_BUILD_REGRESSIONS=ON \
  -DATV_BUILD_HF=OFF -DATV_BUILD_PW=OFF -DATV_BUILD_GPW=OFF \
  -DATV_BUILD_CUDA=OFF -DATV_BUILD_TORCH=OFF \
  -DBLA_VENDOR=OpenBLAS
```

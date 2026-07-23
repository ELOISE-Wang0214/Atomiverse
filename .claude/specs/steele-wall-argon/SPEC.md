# steele-wall-argon

**Branch:** `wf/steele-wall-argon`  
**Chapter:** E8 (salt / wall potentials)  
**Oracle:** `lammps` (mandatory — true Steele is `fix wall/lj1043`, not `wall/lj93`)  
**Status:** open  
**Map:** MetalWalls E8 wall potential; queue row `steele-wall-argon`

## Goal

Ship a minimal **Steele 10-4-3 continuum wall** foothold as a public
`cpu::md::` Bias Component on `cpu::md::MD` / `cpu::Evaluate` — planar wall
interacting with fluid atoms via the integrated LJ form used by LAMMPS
`fix wall/lj1043`. Frame I/O only; no public `*Config` bags; no `make_*` /
`build_*` factories. Units: eV, Å on the public surface.

> **Name note:** the wave-2 queue parenthetical wrote `wall/lj93` Steele 10-4-3.
> LAMMPS `wall/lj93` is the **9-3** integrated half-space; the **10-4-3 Steele**
> form is `wall/lj1043` (Steele, *Surf. Sci.* **36**, 317 (1973); LAMMPS doc).
> This slug implements **10-4-3** / `wall/lj1043` as the method name requires.

## Non-goals

- `wall/lj93` (9-3) and `wall/lj126` (12-6) — separate footholds if needed
- Moving / variable wall position (`v_` styles), `units lattice`, `fld yes`
- Multi-wall simultaneous registration in one component (caller may add N walls)
- Soft-core walls, granular walls, region walls
- CUDA path; SCF / PW / GPW
- Full graphite multi-layer / corrugation (smooth continuum only)

## Public API (binding)

```cpp
// Planar Steele 10-4-3 wall — matches LAMMPS fix wall/lj1043 (shifted at rc).
// axis: 0=x, 1=y, 2=z.  side: -1 = lo (atoms with coord > wall_pos),
//                          +1 = hi (atoms with coord < wall_pos).
cpu::md::FixWallLJ1043 wall(
    /*axis=*/2,
    /*side=*/-1,
    /*wall_pos [Å]=*/0.0,
    /*epsilon [eV]=*/0.010338,
    /*sigma [Å]=*/3.40,
    /*cutoff [Å]=*/12.0);

// Evaluate single-point (or MD post_force via default binding):
cpu::Evaluate eval;
eval.on(evaluate::Phase::evaluate).add(std::move(wall));
eval.build();
eval.run(frame);
// frame.meta.get(kEnergy) [eV]; wall forces on kAtom_f*
```

Constructor arguments establish all parameters (no `WallConfig`).  
Energy is **shifted** so \(U(r_c)=0\) (LAMMPS wall offset).

### Formula (LAMMPS `wall/lj1043`)

\[
U(r)=2\pi\varepsilon\left[
  \frac{2}{5}\left(\frac{\sigma}{r}\right)^{10}
  -\left(\frac{\sigma}{r}\right)^{4}
  -\frac{\sqrt{2}\,\sigma^{3}}{3\bigl(r+(0.61/\sqrt{2})\sigma\bigr)^{3}}
\right]-U(r_c),\qquad 0 < r < r_c
\]

with \(r=|x_{\mathrm{axis}}-x_{\mathrm{wall}}|\) on the fluid side of the plane.
Force is \(-\nabla U\) along the wall normal (LAMMPS sign convention).

## Tasks

1. [x] `FixWallLJ1043` header + impl under `cpu/dynamics/` (Bias, default `md::post_force`)
2. [x] Unit tests: analytic E/F at several \(r\), zero beyond cutoff, reject invalid params
3. [x] Offline LAMMPS pin + twin input under `tools/metalwalls_oracle/inputs/lammps/`
4. [x] `REFERENCE.md` with `oracle: lammps` + `oracle_version:`
5. [x] E2E regression `regressions/forcefield/steele-wall-argon/`
6. [x] Wire CMake regression + unit test targets
7. [x] MD-only build; regression PASS; unit tests PASS
8. [x] Perf gate: `perf_exempt` (single-atom pin; startup-dominated)

## Acceptance

| ID | Criterion | Tol |
| -- | --------- | --- |
| A1 | Oracle is LAMMPS `fix wall/lj1043` (Steele 10-4-3); pin path documented offline | — |
| A2 | Offline twin under `tools/metalwalls_oracle/inputs/lammps/` (in.steele_wall_argon or dir) | — |
| A3 | `REFERENCE.md` has `oracle: lammps` and `oracle_version:` | — |
| A4 | Dual-gate energy: single Ar at \(z=4\) Å matches pin within **1e-10 eV** | 1e-10 eV |
| A5 | Dual-gate force: \(F_z\) matches pin within **1e-10 eV/Å**; \(F_x=F_y=0\) | 1e-10 eV/Å |
| A6 | Unit analytic check at \(r\in\{3,4,5\}\) Å within **1e-12** of closed form | 1e-12 |
| A7 | MD-only cmake: HF/PW/GPW/CUDA/Torch **OFF** | — |
| A8 | Unit tests under `tests/` only; one e2e regression | — |
| A9 | Perf: \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\le 10\) fair wall, or `perf_exempt` | ≤10× |

## Observable

- Wall energy \(E\) [eV]
- Force on Ar \(F_z\) [eV/Å]

## Geometry (regression pin)

| Field | Value |
| --- | --- |
| Atoms | 1 Ar (type 0), mass 39.948 amu |
| Position | \((10,10,4)\) Å |
| Box | \(20\times20\times30\) Å, PBC \(x,y\); wall normal \(z\) |
| Wall | \(zlo=0\), \(\varepsilon=0.010338\) eV, \(\sigma=3.40\) Å, \(r_c=12\) Å |
| Pair | none (wall only) |

## Build (binding)

```bash
bash -lc '
  module load buildenv-gcc/2026.03-mpich
  export OPENBLAS_ROOT="$OPENBLASDIR"
  export RUSTUP_TOOLCHAIN=stable-x86_64-unknown-linux-gnu
  export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
  cmake -B build -DCMAKE_BUILD_TYPE=Release \
        -DATV_BUILD_MD=ON -DATV_BUILD_KSPACE=ON \
        -DATV_BUILD_TESTS=ON -DATV_BUILD_REGRESSIONS=ON \
        -DATV_BUILD_HF=OFF -DATV_BUILD_PW=OFF -DATV_BUILD_GPW=OFF \
        -DATV_BUILD_CUDA=OFF -DATV_BUILD_TORCH=OFF \
        -DBLA_VENDOR=OpenBLAS
  cmake --build build -j$(nproc) --target atomiverse reg_steele-wall-argon test_cpu_fix_wall_lj1043
'
```

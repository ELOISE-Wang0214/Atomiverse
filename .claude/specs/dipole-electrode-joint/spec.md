# Spec: dipole-electrode-joint

**Slug:** `dipole-electrode-joint`  
**Branch:** `wf/dipole-electrode-joint`  
**Map slice:** **E5** — Coupled induced-dipole ⊕ electrode \(Q\) joint solve  
**Queue:** `docs/developer/metalwalls-method-queue.md`  
**Oracle:** **metalwalls** (MW-only; LAMMPS has no SCF PIM ⊕ const-pot joint twin)  
**Depends:** E1 `ConstPotBO` ✓, E4 `PimInducedDipole` ✓  
**Status:** open (implement)

## Goal

Land a **Component** that solves **self-consistently** metallic electrode
charges \(Q\) (constant-potential BO) **and** electrolyte induced dipoles
\(\boldsymbol{\mu}\) (bulk PIM / Thole), the MetalWalls main hardness of
JCP 2022 polarizable interfaces (`dipoles_and_electrodes`). Frame is the
sole I/O currency; no new Driver category; charge solve stays **outside**
`KSpace*`; no public `*Config` bag / `make_*` factory.

## Public API (target)

```cpp
cpu::md::DipoleElectrodeJoint joint(/*eta [1/Å]=*/2.0,
                                    /*alpha_by_type=*/{0., 0., 0., 1.0, 1.5});
joint.add_electrode(/*type=*/1, /*Psi [V]=*/0.0);
joint.add_electrode(/*type=*/2, /*Psi [V]=*/1.0);
joint.thole_a = 0.39;
joint.electroneutrality = true;
joint.joint_tol = 1e-10;
// default_bindings → md::Phase::pre_force @ -50 (before pair/kspace force)
sim.add(std::move(joint));
```

Parameters = ctor args + public members. Writes electrode \(Q\) to
`kAtom_charge`. Exposes diagnostics: per-electrode totals, last \(\mu\),
joint residual / iteration count, neutrality residual.

## Physics (foothold)

Block Gauss–Seidel / alternating BO:

1. **Elastance** \(A\) from Gaussian electrode sites (same kernels as
   `ConstPotBO`; optional 2D \(xy\) images via `xy_image_range`).
2. **RHS** \(b_i = \Phi_i^{\mathrm{perm}} + \Phi_i^{\mu}\) — permanent
   non-electrode charges (Gaussian–point as `ConstPotBO`) **plus** potential
   from electrolyte induced dipoles
   \(\Phi_i^{\mu} = k_e\sum_j \boldsymbol{\mu}_j\cdot\mathbf{r}_{ij}/r_{ij}^3\)
   (Thole-scaled when \(\alpha_i\alpha_j\) active; electrode sites \(\alpha=0\)).
3. Solve \(A Q = \Psi - b\) (+ optional neutrality projection).
4. Write \(Q\) to Frame electrode rows.
5. **PIM SCC** on polarizable non-electrode atoms: permanent field includes
   electrode \(Q\) (point multipole kernels as E4 bulk PIM) + mutual \(T\cdot\mu\).
6. Outer loop until \(\max|\Delta Q|\) and \(\max|\Delta\mu|\) below `joint_tol`.

**Reductions:** \(\alpha=0\) → pure BO \(Q\) (matches `ConstPotBO` on same geom);
no polarizable sites → electrode-only.

**Out of scope (this slug):** MaZe charge dynamics, ASPC warm-start on joint
variables, production neighbor-list multipole path, full MW 2D reciprocal
electrode matrix, AIM/DAIM, forces on frozen graphene walls, CUDA.

## Geometry (regression)

Minimal **IL + graphene parallel-plate** toy (not a full molten-salt run):

| Site | type | role | \(q\) [e] | \(\alpha\) [Å³] | \((x,y,z)\) [Å] |
| ---- | ---- | ---- | --------- | --------------- | --------------- |
| bot wall | 1 | electrode \(\Psi=0\) | solved | 0 | (0, 0, 2) |
| top wall | 2 | electrode \(\Psi=1\) V | solved | 0 | (0, 0, 18) |
| cation | 3 | electrolyte | +1 | 1.0 | (0, 0, 8) |
| anion | 4 | electrolyte | −1 | 1.5 | (0, 0, 12) |

Box orthorhombic \(L_x=L_y=10\) Å, \(L_z=20\) Å; \(\eta=2\) Å⁻¹;
`xy_image_range=0` (min-image); `thole_a=0.39`; electroneutrality on.

## Offline pin

- **Oracle:** metalwalls (form of Coretti et al. JCP 2022 / MW multipole
  + Gaussian electrode BO; no live MW binary required at ctest).
- **Pin path:** `tools/metalwalls_oracle/inputs/metalwalls/dipole-electrode-joint/pin_joint.py`
  and hard-coded values in `regressions/md/dipole-electrode-joint/REFERENCE.md`.
- **Observables:** \(Q_{\mathrm{bot}}\), \(Q_{\mathrm{top}}\), \(\mu\) of cation/anion,
  joint residual, and \(\alpha=0\) reduction vs `ConstPotBO`.

## Build (binding)

MD-only: `ATV_BUILD_MD=ON`, `ATV_BUILD_KSPACE=ON`, tests+regressions ON;
**HF/PW/GPW/CUDA/Torch OFF**.

## Non-goals

- New Driver category / `ElectrochemistryDriver`
- Public `*Config` bags or `make_*` / `build_*` factories
- Absorbing joint solve into `KSpace*`
- Live MetalWalls or LAMMPS at ctest time

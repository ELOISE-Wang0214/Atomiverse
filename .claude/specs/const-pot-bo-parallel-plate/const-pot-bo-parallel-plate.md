# Spec: const-pot-bo-parallel-plate

**Status:** open  
**Branch:** `wf/const-pot-bo-parallel-plate`  
**Map slice:** E1 (constant-potential v1, map §1 P0)  
**Method:** Constant-potential Born–Oppenheimer (dense matrix inverse) electrode charge solve  
**Geometry:** Parallel-plate capacitor foothold — two electrode sites + two electrolyte ions (LAMMPS ELECTRODE `madelung` twin); also supports multi-atom planar walls by type  

## Goal

Land a **Component** on `cpu::md::MD` that, each force evaluation, solves electrode
charges \(Q\) at fixed applied potentials \(\Psi\) under the Gaussian BO
elastance matrix \(A\), with optional electroneutrality, and writes
`frame.span(kAtom_charge)`.

Not in scope: CG large electrodes (E2), Gaussian-vs-Ewald full kspace electrode
matrix (E2), freeze integrator (sibling slug), charging-work energy channel
(sibling), Thomas–Fermi, MaZe, PIM coupling.

## Architecture (binding)

| Rule | Choice |
| --- | --- |
| Driver | **None new** — Component on existing `cpu::md::MD` |
| Kind | `Component` (not `Potential`, not inside `KSpace*`) |
| Phase | `md::Phase::pre_force` (early order) so pair/kspace see updated \(Q\) |
| Parameters | ctor args + public members — **no** `*Config` bag, **no** `make_*` |
| Workspace | private \(A\), \(Q\), electrode index maps |
| I/O | Frame only: read \(x\), electrolyte \(q\), types; write electrode \(q\) |
| Units | \(\Psi\) [V], \(Q\) [e], \(r\) [Å], energy channel later in eV |

## Physics (v1)

Gaussian electrode charges of width \(\eta\) [1/Å] (Siepmann–Sprik / MetalWalls /
LAMMPS ELECTRODE form):

\[
A_{ii} = k\,\eta\sqrt{2/\pi},\qquad
A_{ij} = k\,\frac{\mathrm{erf}(\eta_{ij} r_{ij})}{r_{ij}}\ (i\neq j),
\quad
\eta_{ij}=\frac{\eta_i\eta_j}{\sqrt{\eta_i^2+\eta_j^2}},
\quad k=\mathtt{units::coulomb}.
\]

Electrolyte potential at electrode site \(i\):

\[
b_i = \sum_{j\notin\mathrm{elec}} k\, q_j\,\frac{\mathrm{erf}(\eta_i r_{ij})}{r_{ij}}.
\]

BO solution with optional electroneutrality (LAMMPS `symm on`):

\[
Q = A_{\mathrm{sym}}^{-1}(\Psi - b),\quad
A_{\mathrm{sym}}^{-1} = A^{-1} - \frac{A^{-1}ee^\top A^{-1}}{e^\top A^{-1}e}.
\]

Optional **2D lattice images** in \(xy\) (planar electrochemical cell, matching
the ELECTRODE `madelung` / `plate_cap.py` reference sum). Full 3D Ewald
electrode matrix is deferred to E2 (`slab-2d-ewald-foothold` / kspace electrode).

## Public API (target)

```cpp
cpu::md::ConstPotBO conp(/*eta=*/2.0);
conp.add_electrode(/*type=*/1, /*Psi [V]=*/0.0);
conp.add_electrode(/*type=*/2, /*Psi [V]=*/1.0);
conp.electroneutrality = true;
conp.xy_image_range = 200;   // 0 = primary cell only

cpu::md::MD sim;
sim.on(md::Phase::pre_force).add(std::move(conp));
// … neighbor + coulomb/kspace as needed …
```

Accessors after launch: `total_charge(electrode_index)`, `delta_psi()`,
`neutrality_residual()`, `solve_residual()`.

## Tasks

1. Spec + acceptance (this file + `.acceptance.md`).
2. Offline pin: LAMMPS ELECTRODE madelung log + literature plate_cap Madelung sum;
   `REFERENCE.md` with `oracle:` / `oracle_version:`.
3. Implement `cpu::md::ConstPotBO` (header + cpp), register in CMake MD sources.
4. Unit tests: neutrality, residual, small no-image Q.
5. E2e regression `regressions/md/const-pot-bo-parallel-plate/`.
6. MD-only build + run regression to pass.
7. Perf gate: no ELECTRODE-enabled twin binary → `perf_exempt` with reason.
8. Commit on `wf/const-pot-bo-parallel-plate`.

## Dual-gate

| Gate | Criterion |
| --- | --- |
| Correctness | Regression vs offline pin (see acceptance) |
| Speed | Fair wall vs same oracle binary if available; else exempt |

## Oracle decision

- Algorithm **exists in LAMMPS** (`fix electrode/conp`, package ELECTRODE).
- Local `lmp` **lacks** ELECTRODE; offline pin from shipped example log
  `examples/PACKAGES/electrode/madelung/log.1Dec2022.madelung-ewald-ew3dc.g++.1`
  plus literature Madelung sum (`plate_cap.py`).
- **oracle = lammps** (offline pin); MetalWalls not required for this twin.

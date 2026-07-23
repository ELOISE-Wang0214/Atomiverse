# Spec: const-pot-cg-large-electrode

**Status:** open  
**Branch:** `wf/const-pot-cg-large-electrode`  
**Map slice:** E2 (electrostatic rigor — CG large electrodes)  
**Method:** Conjugate-gradient electrode charge solve (vs dense BO inverse)  
**Geometry:** Larger planar wall — many electrode sites per plate + electrolyte ions  

## Goal

Land a **Component** on `cpu::md::MD` that solves the constant-potential electrode
linear system \(A Q = \Psi - b\) with **preconditioned conjugate gradient**
instead of dense Cholesky / full inverse, so large \(N_{\mathrm{elec}}\) walls
scale as \(O(N_{\mathrm{iter}} N_{\mathrm{elec}}^2)\) (or better with future
matrix-free \(A\)) rather than \(O(N_{\mathrm{elec}}^3)\).

Primary correctness gate: **site charges \(Q\) match the dense BO path**
(`cpu::md::ConstPotBO`) on the same geometry, \(A\), \(\eta\), \(\Psi\), and
electroneutrality projector, within a tight numerical tolerance.

Not in scope: ASPC charge predictor (sibling slug), Thomas–Fermi, MaZe, PIM,
charge solve inside `KSpace*`, CUDA.

## Architecture (binding)

| Rule | Choice |
| --- | --- |
| Driver | **None new** — Component on existing `cpu::md::MD` |
| Kind | `Component` (not `Potential`, not inside `KSpace*`) |
| Phase | `md::Phase::pre_force` (early; same role as `ConstPotBO`) |
| Parameters | ctor args + public members — **no** `*Config` bag, **no** `make_*` |
| Workspace | private \(A\), \(Q\), CG vectors; electrode index maps |
| I/O | Frame only: read \(x\), electrolyte \(q\), types; write electrode \(q\) |
| Units | \(\Psi\) [V], \(Q\) [e], \(r\) [Å], \(\eta\) [1/Å] |

## Physics

Same Gaussian elastance \(A\) and electrolyte RHS \(b\) as E1 `ConstPotBO`
(Siepmann–Sprik / MetalWalls / LAMMPS ELECTRODE form). Unconstrained solve:

\[
A Q = \Psi - b.
\]

With **electroneutrality** (MetalWalls / LAMMPS `symm on`):

\[
Q = v - \frac{e^{\mathsf T} v}{e^{\mathsf T} u}\, u,\quad
A v = \Psi-b,\quad A u = e,
\]

where each linear solve uses **Jacobi-preconditioned CG** to residual
\(\lVert A x - rhs\rVert_2 / \lVert rhs\rVert_2 \le\) `cg_tol` (or absolute
when \(\lVert rhs\rVert\) is tiny).

Optional 2D \(xy\) lattice images via `xy_image_range` (same semantics as BO).

## Public API (target)

```cpp
cpu::md::ConstPotCG conp(/*eta [1/Å]=*/2.0);
conp.add_electrode(/*type=*/1, /*Psi [V]=*/0.0);
conp.add_electrode(/*type=*/2, /*Psi [V]=*/1.0);
conp.electroneutrality = true;
conp.xy_image_range = 0;   // primary / min-image for large multi-site walls
conp.cg_tol = 1e-12;
conp.max_cg_iter = 500;

cpu::md::MD sim;
sim.on(md::Phase::pre_force).add(std::move(conp));
// or sim.add_ref(conp);
```

Accessors after launch: `total_charge(ie)`, `delta_psi()`,
`neutrality_residual()`, `solve_residual()`, `n_electrode()`, `cg_iterations()`
(last CG iter count for the primary solve).

## Tasks

1. Spec + acceptance (this file + `acceptance.md`).
2. Offline dual pin: dense BO twin on multi-site plate (primary) + document
   oracle=metalwalls (MW-only CG path; no LAMMPS CG electrode twin required);
   `REFERENCE.md` with `oracle:` / `oracle_version:`.
3. Implement `cpu::md::ConstPotCG` (header + cpp), register in CMake MD sources.
4. Unit tests: CG matches dense BO on small madelung; neutrality; larger plate.
5. E2e regression `regressions/md/const-pot-cg-large-electrode/`.
6. MD-only build + run regression to pass.
7. Perf gate: no MetalWalls twin binary → `perf_exempt` with reason; optional
   wall vs dense BO on same geometry recorded in REFERENCE.
8. Commit on `wf/const-pot-cg-large-electrode`.

## Dual-gate

| Gate | Criterion |
| --- | --- |
| Correctness | Site \(Q\) vs dense `ConstPotBO` (see acceptance) |
| Speed | Fair wall vs MetalWalls CG if binary available; else exempt |

## Oracle decision

- Dense inverse / ELECTRODE `conp` exists in LAMMPS, but **iterative CG electrode
  solve** is the MetalWalls-style large-\(N\) path (map E2). Queue marks
  **oracle=metalwalls**.
- Acceptance dual-gate uses **in-tree dense BO** as the numeric twin (same \(A\))
  so ctest never needs a live MW binary.
- Offline pin path: `regressions/md/const-pot-cg-large-electrode/REFERENCE.md`.

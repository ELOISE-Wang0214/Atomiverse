# Spec — const-charge-electrode

**Slug:** `const-charge-electrode` (metalwalls-method-loop, wave 5)
**Method:** Constant-charge electrode control (fixed total \(Q\) per electrode,
potential \(\Psi\) is the response) — the LAMMPS ELECTRODE `fix electrode/conq`
/ MetalWalls constant-charge twin of the existing constant-potential BO solve.
**Map row:** §1 "Constant-charge electrodes (control)" — status `~` (foothold
was manual fixed charge + `FixFreeze`; no self-consistent const-\(Q\)
electrode type).

## Problem

`ConstPotBO` solves \(A Q = \Psi - b\) at **fixed applied potential**. The
dual control mode — prescribe each electrode's **total charge**
\(\sum_{i\in e} Q_i = Q_e\) and let the electrode potential float — is absent.
Fixed per-site charge + `FixFreeze` cannot distribute a total charge
self-consistently over the electrode sites, and gives no potential readout.

## Design (architecture-law compliant)

Extend the existing Component `cpu::md::ConstPotBO` — **no** new Driver
category, **no** charge solve inside `KSpace*`, **no** `*Config` bag, no
factories. Reuse the Gaussian elastance machinery (`build_matrix_`,
`build_rhs_`, dense Cholesky + CG paths) — zero duplication.

Math (LAMMPS `electrode/conq`): with indicator \(d_e\) of const-\(Q\) group
\(e\), unknown uniform group potentials \(\psi_e\),

\[
Q = A^{-1}\bigl(\Psi_{\text{base}} + \textstyle\sum_e \psi_e d_e - b\bigr),
\qquad d_e^{\mathsf T} Q = Q_e .
\]

Solve the small SPD system \(M\psi = r\) with
\(M_{ef} = d_e^{\mathsf T} A^{-1} d_f\),
\(r_e = Q_e - d_e^{\mathsf T} A^{-1}(\Psi_{\text{base}} - b)\); then
\(Q = Q_0 + \sum_e \psi_e u_e\) with \(u_e = A^{-1} d_e\). Works on both the
dense-Cholesky and CG paths.

### Public surface (members + methods on `ConstPotBO`, no new types)

- `ElectrodeGroup` gains `bool const_charge` + `double q_target` (defaults
  keep existing aggregate init valid).
- `add_electrode_charge(int type, double q_total)` (`&` and `&&` chaining,
  mirroring `add_electrode`) — registers a constant-charge group.
- `electrode_potential(std::size_t ie)` [V] — applied \(\Psi\) for const-Ψ
  groups; **solved** \(\Psi\) for const-Q groups after `launch`.
- `delta_psi()` now reports solved potentials (identical to applied for pure
  const-Ψ runs — existing behaviour preserved bit-for-bit).
- When any const-Q group is registered the `electroneutrality` projection is
  skipped: the group-total constraints are the charge control (targets summing
  to \(-q_{\text{elyt}}\) give global neutrality). Documented in the header.

## Acceptance criteria

1. **Exact charge control:** on the madelung parallel-plate twin, conq with
   targets \(\pm Q\) gives `total_charge(e)` equal to the targets to
   \(10^{-12}\) e (dense) — "Q fixed".
2. **V response vs LAMMPS oracle (inverse pin):** feeding the LAMMPS ELECTRODE
   madelung conp pin \(Q_{\text{top}} = +0.279214485147238\) e (log
   `log.1Dec2022.madelung-ewald-ew3dc.g++.1`, LAMMPS 3 Nov 2022) as the conq
   target must return \(\Delta\Psi = \Psi_{\text{top}}-\Psi_{\text{bot}} =
   1.000\) V within **1 %** (v1 2D-image Gaussian \(A\) vs full
   `ewald/electrode`, same gate as `const-pot-bo-parallel-plate`).
3. **conp/conq duality (multi-site):** on a multi-site-per-plate geometry with
   asymmetric electrolyte, conq fed the conp group totals reproduces the conp
   **per-site** charges and applied \(\Delta\Psi\) to \(10^{-8}\) (unit test —
   the distribution over sites is self-consistent, not uniform).
4. **CG path parity:** `use_cg = true` conq matches dense conq to \(10^{-8}\).
5. **Existing const-Ψ behaviour untouched:** all existing `ConstPotBO` /
   `ConstPotCG` unit tests and `reg_const-pot-bo-parallel-plate` still pass.
6. **Regression:** `regressions/md/const-charge-electrode/` out-of-library
   consumer + `REFERENCE.md` (oracle: lammps, version, cited numbers, units),
   targets `reg_const-charge-electrode` / `run_reg_const-charge-electrode`,
   ≤ 60 s.
7. **Unit tests** under `tests/cpu/dynamics/test_const_charge_electrode.cpp`,
   registered in `tests/CMakeLists.txt`, seconds-fast, single-function each.
8. **MD-only build** stays green (HF/PW/GPW/CUDA/TORCH OFF).

## Perf

Local `lmp` lacks `PKG_ELECTRODE` (verified 2026-07-27: installed packages =
EXTRA-*, KSPACE, MANYBODY, MISC, MOLECULE, OPENMP, QEQ, REPLICA, RIGID), so no
fair `fix electrode/conq` twin binary → `perf_exempt`, same as
`const-pot-bo-parallel-plate`.

## Out of scope

Full `ewald/electrode` reciprocal matrix (E2 foothold exists), per-electrode
\(\eta\) (separate wave-5 slug), MPI halo for solved charges (separate slug),
ramped/time-dependent \(Q(t)\) schedules.

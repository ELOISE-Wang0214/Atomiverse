# Spec: slab-2d-ewald-foothold

**Slug:** `slab-2d-ewald-foothold`  
**Map slice:** E2 — electrostatic rigor (2D / slab Ewald)  
**Branch:** `wf/slab-2d-ewald-foothold`  
**Oracle:** **LAMMPS** (`kspace_style ewald` + `kspace_modify slab`) — mandatory on overlap  
**Status:** open

## Goal

Ship a **LAMMPS-compatible quasi-2D (slab) Ewald** path for electrochemical
planar geometry: 3D Ewald on a z-extended cell plus the Yeh–Berkowitz
\(k=0\) dipole correction (`kspace_modify slab <volfactor>`).

This is **not** a new Driver category. It is a public member on the existing
reciprocal component `cpu::KSpaceEwald` (parameters = ctor args / members; no
`*Config` bag, no `make_*` factory). Charge solves stay outside `KSpace*`.

## Physics (binding)

Yeh & Berkowitz, *J. Chem. Phys.* **111**, 3155 (1999); non-neutral extension
Ballenegger et al., *J. Chem. Phys.* **131**, 094107 (2009). Matches LAMMPS
`Ewald::setup` + `Ewald::slabcorr` (develop / 30 Mar 2026 tree):

1. Extended period: \(L_z' = L_z \cdot f\), \(V' = L_x L_y L_z'\),
   \(k_z = 2\pi m_z / L_z'\).
2. Standard reciprocal + self + net-charge terms use \(V'\).
3. Dipole correction (neutral system \(Q=0\)):
   \[
   E_{\mathrm{slab}} = k_e \frac{2\pi}{V'} M_z^2,\quad
   M_z = \sum_i q_i z_i
   \]
   \[
   F_{z,i} = -k_e \frac{4\pi}{V'} q_i M_z
   \]
   with \(k_e =\) `units::coulomb` (eV·Å/e²). Non-neutral terms follow LAMMPS
   `slabcorr` when \(Q\neq 0\).

`slab_factor <= 1` (default `0`) disables the path (pure 3D Ewald).

## Public API (minimal)

```cpp
cpu::KSpaceEwald recip(/*alpha=*/0.5, /*accuracy=*/1e-6);
recip.slab_factor = 3.0;   // public member; LAMMPS kspace_modify slab 3.0
// pair with cpu::CoulLong{alpha, cutoff} as today
```

No public `SlabConfig`, no free `make_slab_ewald`.

## Tasks

1. [x] Spec + acceptance (this file + `acceptance.md`).
2. [x] Offline LAMMPS pin under `tools/metalwalls_oracle/inputs/lammps/` +
      numbers in `regressions/md/slab-2d-ewald-foothold/REFERENCE.md`
      (`oracle: lammps`, `oracle_version:`).
3. [x] Implement `KSpaceEwald::slab_factor` (volume scale + Yeh–Berkowitz
      energy/force); pure math may live in `longrange/` header if shared.
4. [x] Unit tests under `tests/cpu/longrange/` (neutral Mz correction; FD force
      with slab; 3D path unchanged when `slab_factor==0`).
5. [x] E2E regression `regressions/md/slab-2d-ewald-foothold/` vs offline pin
      (energy + \(F_z\) sign/magnitude).
6. [x] MD-only build (HF/PW/GPW/CUDA/Torch **OFF**); regression PASS.
7. [x] Perf gate vs same LAMMPS binary (or `perf_exempt` with reason).
8. [ ] Commit on `wf/slab-2d-ewald-foothold` only.

## Out of scope

- True 2D Ewald / MetalWalls native 2D (later; this foothold is 3D+slab).
- PME/PPPM slab (`kspace_modify slab` on PPPM) — follow-up.
- Constant-potential electrodes (E1), freeze, efield (other slugs).
- CUDA.

## Dual-gate

| Gate | Criterion |
|---|---|
| Correctness | Regression energy within stated tol of LAMMPS offline pin; unit FD forces |
| Speed | \(T_{\mathrm{atv}}/T_{\mathrm{lammps}} \le 10\) on fair twin, or exempt |

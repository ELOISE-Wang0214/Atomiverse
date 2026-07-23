# Offline pin notes — `piston-electrode-pz`

**oracle:** metalwalls (MW-only rigid electrode piston; no LAMMPS twin)

**oracle_version:** MetalWalls algorithm as published in
Marin-Laflèche et al., JOSS 5, 2373 (2020) DOI 10.21105/joss.02373 —
*constant pressure \(P_z\) by treating electrodes as rigid pistons*.

## Why not LAMMPS

LAMMPS `fix wall/piston` (SHOCK package) is a **shock-wave** driven wall,
not an electrode rigid body with \(F_{\mathrm{ext}}=\pm P A\) for
electrochemical constant-\(P_z\). Capability map + loop oracle rules mark
this slug **oracle=metalwalls**.

## Pin procedure (no MW binary on worker)

1. Build the parallel-plate LJ geometry in
   `regressions/md/piston-electrode-pz/reg_piston-electrode-pz.cpp`.
2. Probe step with \(P_{\mathrm{target}}=0\): record
   \(P_z^{\mathrm{probe}}=F_z^{\mathrm{atoms}}/A\).
3. Production with \(P_{\mathrm{target}}=P_z^{\mathrm{probe}}\) and COM damping:
   assert \(\lvert\langle P_z\rangle-P_{\mathrm{target}}\rvert/P_{\mathrm{target}}\le 0.10\).
4. Unit tests pin \(F_{\mathrm{net}}=0\) and \(P_z=F/A\) identities to 1e-12.

When a MetalWalls binary is available offline, re-run a comparable piston
example from `ampere2/metalwalls` `example/` and replace the self-consistent
probe pin with tabulated MW \(\langle P_z\rangle\), then fair-bench walls.

## Numbers live in

`regressions/md/piston-electrode-pz/REFERENCE.md`

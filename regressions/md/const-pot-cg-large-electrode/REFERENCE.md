# Reference — constant-potential CG solve (large multi-site electrode)

oracle: metalwalls  
oracle_version: MetalWalls-style Jacobi-PCG electrode charge solve (map E2); dual twin = in-tree `cpu::md::ConstPotBO` dense path on identical Gaussian elastance \(A\) (2026-07-24)

## System

Larger planar capacitor with many electrode sites (queue geom: larger wall /
many electrode sites):

| Species | Type | Geometry | initial q [e] |
| --- | --- | --- | --- |
| Bottom electrode | 1 | \(4\times 4\) grid, \(z=8\) Å | 0 (solved) |
| Top electrode | 2 | \(4\times 4\) grid, \(z=12\) Å | 0 (solved) |
| Electrolyte + | 3 | box center \(xy\), \(z=9\) Å | +0.5 |
| Electrolyte − | 3 | box center \(xy\), \(z=11\) Å | −0.5 |

- \(N_{\mathrm{elec}} = 32\) (16 + 16)
- Box: \(L_x=L_y=4\) Å, \(L_z=20\) Å (orthorhombic, fully periodic ATV)
- Applied potentials: \(\Psi_{\mathrm{bot}}=0\) V, \(\Psi_{\mathrm{top}}=1\) V
- Gaussian width \(\eta=2.0\) Å⁻¹
- Electroneutrality on (`symm`)
- Lattice images: primary / min-image (`xy_image_range = 0`) for multi-site wall

## Offline pin / dual twin

| Source | Role |
| --- | --- |
| **Dense twin** `cpu::md::ConstPotBO` | Binding numeric oracle for \(Q_i\) (same \(A\), \(\eta\), \(\Psi\), neutrality projector) |
| MetalWalls CG | Algorithm class (iterative electrode solve for large \(N\)); no live MW binary required at ctest |

Primary-cell madelung sanity (2 electrode atoms, same as E1 unit test):

| Observable | Value |
| --- | --- |
| \(Q_{\mathrm{top}}\) [e] (dense BO / CG, `xy_image_range=0`) | **0.27175356** |

## What Atomiverse asserts

Public API: `cpu::md::ConstPotCG` on `md::Phase::pre_force`, writes
`kAtom_charge`.

| Check | Tol |
| --- | --- |
| \(N_{\mathrm{elec}} \ge 18\) | exact |
| \(\lvert\sum Q_i\rvert\) | \(\le 10^{-10}\) e |
| \(\max_i \lvert Q_i^{\mathrm{CG}} - Q_i^{\mathrm{BO}}\rvert\) | \(\le 10^{-8}\) e |
| \(\lvert Q_{\mathrm{group}}^{\mathrm{CG}} - Q_{\mathrm{group}}^{\mathrm{BO}}\rvert\) | \(\le 10^{-8}\) e |
| Primary \(Q_{\mathrm{top}}\) vs 0.27175356 | abs \(\le 10^{-6}\) e |

## Perf

```
perf_exempt: true
perf_exempt_reason: no MetalWalls runtime binary on worker PATH for fair CG twin; dual correctness twin is in-tree ConstPotBO (not an external oracle wall)
ratio_atv_over_oracle: n/a
within_10x: exempt
```

Internal informational wall (CG vs dense BO on the 32-site plate) is printed by
the regression binary (`wall_ms CG` / `wall_ms BO`) and is **not** the 10×
oracle gate.

## Repro of dual twin (offline / local only — never ctest)

```bash
# Build MD-only, then:
./build/reg_const-pot-cg-large-electrode
# asserts CG vs BO site charges; no LAMMPS/MetalWalls launch
```

MetalWalls (when available offline): run an equivalent multi-site constant-potential
CG electrode demo and record group \(Q\) into this file; keep ctest free of the binary.

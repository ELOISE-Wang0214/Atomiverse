# surface-tension-coupling

## Goal

Ship **constant surface-tension (NPγT) ensemble coupling** on CPU classical MD:
semi-isotropic MTK barostat with **LAMMPS `fix npt … couple xy`** lockstep
lateral edges, plus a target surface tension \(\gamma\) that sets the lateral
pressure offset

\[
P_\perp^\mathrm{target}(t) = P_n - \frac{n_\mathrm{surf}\,\gamma}{L_z(t)},
\quad P_z^\mathrm{target} = P_n.
\]

Oracle: **LAMMPS** (`fix npt` surface-tension path). Geometry: LJ slab / bulk
with anisotropic targets. Observables: \(\gamma\) and/or \(L_z\) (and \(L_x=L_y\)).

## Non-goals

- CUDA / Torch barostat parity
- Full triclinic shear (`xy`/`xz`/`yz` tilt DOFs) — that is PR-NPT
- Live LAMMPS at ctest time
- SCF / GPW / PW

## Public API

```cpp
// Dedicated preset over NoseHooverNPT (no *Config bag, no make_*)
cpu::md::SurfaceTensionNPT st{
    /*T [K]=*/100.0,
    /*tau_T [fs]=*/50.0,
    /*P_normal [eV/Å³]=*/0.0,
    /*gamma [eV/Å²]=*/0.001,
    /*tau_P [fs]=*/500.0,
    /*dt [fs]=*/1.0,
    /*n_surfaces=*/2};

// Or configure on NoseHooverNPT:
cpu::md::NoseHooverNPT npt{T, tau_T, P_n, tau_P, dt};
npt.set_couple_xy(true);
npt.set_surface_tension(/*gamma=*/0.001, /*n_surfaces=*/2);
```

- `SurfaceTensionNPT` **is** an `Integrator` (via `NoseHooverNPT`); replaces
  VelocityVerlet / NoseHoover / Barostat on `cpu::md::MD`.
- Instantaneous planar surface tension (diagnostic, not a second driver):

\[
\gamma_\mathrm{inst} = \frac{L_z}{2}\Bigl(P_{zz}-\tfrac12(P_{xx}+P_{yy})\Bigr)
\]

  with \(P_{\alpha\alpha} = (\sigma_{\alpha\alpha}^\mathrm{vir}+2\,\mathrm{KE}_\alpha)/V\)
  (LAMMPS `compute pressure` sign: positive compressive).

## Tasks

1. Extend `cpu::md::NoseHooverNPT` with:
   - `set_couple_xy(bool)` — LAMMPS `couple xy` (average \(P_{xx},P_{yy}\);
     lockstep \(\dot\varepsilon_x=\dot\varepsilon_y\); \(z\) independent).
   - `set_surface_tension(gamma, n_surfaces=2)` — enables couple-xy and the
     \(P_\perp(L_z)\) offset above (`P_target` = normal pressure \(P_n\)).
2. Add public `cpu::md::SurfaceTensionNPT` constructor that wires the above
   (parameters as ctor args; no Config bag).
3. Unit tests under `tests/cpu/dynamics/` (API + couple-xy \(L_x=L_y\) short run).
4. Offline LAMMPS pin + `regressions/md/surface-tension-coupling/` e2e.
5. Fair wall vs LAMMPS twin; record ratio in REFERENCE.md.

## Design notes

- When `couple_xy` is on, **do not** apply the full-ISO `pdim==3` stress average;
  only average the lateral pair.
- `gamma == 0` + couple-xy ⇒ semi-isotropic isobaric NPT (same \(P\) on all
  axes, lockstep \(xy\)) — valid intermediate.
- Dynamic \(L_z\) update of \(P_\perp\) is the constant-\(\gamma\) ensemble;
  LAMMPS `fix npt` takes fixed \(P_{start}/P_{stop}\) only, so the offline pin
  uses **fixed** \(P_\perp = P_n - n\gamma/L_z(0)\) for a short twin (standard
  community path for LAMMPS surface-tension NPT).

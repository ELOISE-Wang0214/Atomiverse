# Spec — `cpmd-rattle-water`

**Slug:** `cpmd-rattle-water`  
**Domain:** dynamics (GPW / AIMD)  
**Status:** in progress  
**Branch:** `wf/cpmd-rattle-water`  
**Parent:** `dev` (after SHAKE-free `cpmd-water` foothold)

## Intent

Upgrade Car–Parrinello electronic propagation from the **Löwdin re-orthonormalization
foothold** (`cpmd-water`) to **full electronic RATTLE**: constraint-preserving
position and velocity stages for the occupied MO coefficients under
\(C^\top S C = I\). Same closed-shell Γ GPW water NVE envelope as the foothold,
with **tighter** \(E_{\mathrm{cons}}\) drift and explicit residual gates on the
constraint manifold.

## Literature

- Car & Parrinello, *Phys. Rev. Lett.* **55**, 2471 (1985)
- Andersen, *J. Comput. Phys.* **52**, 24 (1983) — RATTLE (velocity + position)
- Marx & Hutter, *Ab Initio Molecular Dynamics* (CUP) — CP Lagrangian,
  holonomic orthonormality, conserved
  \(E_{\mathrm{cons}} = T_{\mathrm{n}} + T_{\mathrm{e}} + E_{\mathrm{KS}}\)
- Twin family: `cpmd-water` (Löwdin foothold), `gpw-bomd-water`, `xl-bomd-water`

## Public API (no Config bags / factories)

| Symbol | Role |
| ------ | ---- |
| `cpu::gpw::Cpmd` | Electronic DOF propagator — **RATTLE** position + velocity (replaces Löwdin-only `advance`) |
| `GPWEvaluator::{cpmd,cpmd_mu,cpmd_dt}` | Unchanged engagement (members, not Config) |
| `Cpmd::last_position_iters()` | Diagnostics: RATTLE position iterations |
| `Cpmd::last_constraint_residual()` | \(\lVert C^\top S C - I\rVert_\infty\) after last position stage |

Compose with `cpu::md::MD` + `VelocityVerlet` at `md::Phase::force` (same as foothold).

## Algorithm

**Constraint** (holonomic, matrix-valued):

\[
\sigma(C) = C^\top S C - I = 0.
\]

**Electronic force** (unchanged, S-tangent OT gradient, closed-shell \(P=2CC^\top\)):

\[
F = -(I - S C C^\top)(4 H C),\qquad \mu\,\ddot C = F.
\]

**Velocity-Verlet + electronic RATTLE** (Andersen form on MO coefficients):

1. Half-step velocity: \(\dot C_{1/2} = \dot C + (\Delta t_{\mathrm{au}}/(2\mu)) F\)
2. Unconstrained position: \(C_{\mathrm{raw}} = C + \Delta t_{\mathrm{au}}\,\dot C_{1/2}\)
3. **Position RATTLE (SHAKE):** \(C' = C_{\mathrm{raw}} + C\,\Lambda\) with iterative
   \(\Lambda\) so \(C'^\top S C' = I\) (reference \(C\) is the pre-step orthonormal set).
4. Redefine half-step from constrained displacement:
   \(\dot C_{1/2} \leftarrow (C' - C)/\Delta t_{\mathrm{au}}\)
5. Complete velocity: \(\dot C_{\mathrm{unc}} = \dot C_{1/2} + (\Delta t_{\mathrm{au}}/(2\mu)) F\)
   (single force evaluation; matches fixed-density CPMD launch)
6. **Velocity RATTLE:** \(\dot C' = \dot C_{\mathrm{unc}} + C'\,\Gamma\) with
   \(\Gamma = -\mathrm{sym}(C'^\top S \dot C_{\mathrm{unc}})\) so
   \(C'^\top S \dot C'\) is antisymmetric (\(\dot\sigma = 0\)).
7. \(T_e = \tfrac12\mu\lVert\dot C'\rVert_F^2\) [Ha]; store \(C\leftarrow C'\), \(\dot C\leftarrow\dot C'\).

**Geometry change** (`align_to_overlap`): Löwdin transport of \(C\) under the new \(S\),
then velocity RATTLE projection (metric change is not a dynamics step).

**Nuclear path** (unchanged): cold SCF seed → fixed-density KS + analytic forces
→ electronic advance; production `iters==0`.

## Scope / non-goals

| In | Out |
| -- | --- |
| Closed-shell RKS Γ GPW, PBE, MT, DZVP water | UKS, OT, k-points |
| Electronic RATTLE (pos + vel) | Nuclear RATTLE (`FixRattle` is classical MD) |
| Tighter \(E_{\mathrm{cons}}\) envelope vs Löwdin foothold | Electronic thermostat |
| Unit residual gates on \(\sigma,\dot\sigma\) | Production multi-ps traj / live CP2K in ctest |
| Offline CP2K twin input (oracle pin notes) | Silent fallback to BOMD |

## Acceptance (summary)

See `acceptance.md`. Unit: RATTLE residual + velocity antisymmetry. E2E:
`regressions/dynamics/cpmd-rattle-water` — CP \(E_{\mathrm{cons}}\) drift
\(\le 0.05\), constraint residual, `iters==0`, molecule intact.

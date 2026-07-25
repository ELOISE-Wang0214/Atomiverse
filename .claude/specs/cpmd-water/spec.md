# Spec — `cpmd-water`

**Slug:** `cpmd-water`  
**Domain:** dynamics (GPW / AIMD)  
**Status:** in progress  
**Branch:** `wf/cpmd-water`

## Intent

Ship a **Car–Parrinello MD (CPMD) foothold** on closed-shell Γ-point GPW:
fictitious electronic degrees of freedom (occupied MO coefficients) are
propagated with a Verlet-like step and fictitious mass \(\mu\), coupled to
nuclear Velocity Verlet. This is **not** BOMD (no SCF each step after seed)
and **not** XL-BOMD (no auxiliary density + SCF).

## Literature

- Car & Parrinello, *Phys. Rev. Lett.* **55**, 2471 (1985) — extended
  Lagrangian with fictitious electronic kinetic energy.
- Marx & Hutter, *Ab Initio Molecular Dynamics* (CUP) — CP conserved energy
  \(E_{\mathrm{cons}} = T_{\mathrm{n}} + T_{\mathrm{e}} + E_{\mathrm{KS}}\).
- CP2K Quickstep AIMD / CPMD docs (electronic mass, short H₂O NVE).

## Public API (no Config bags / factories)

| Symbol | Role |
| ------ | ---- |
| `cpu::gpw::Cpmd` | Electronic DOF propagator: \(\mu\) [a.u. mass], Verlet + Löwdin reorth |
| `GPWEvaluator::cpmd` | Engage CPMD path (mutually exclusive with `xl_bomd`) |
| `GPWEvaluator::cpmd_mu` | Fictitious electron mass [a.u.] (default 400) |
| `GPWEvaluator::cpmd_dt` | Electronic timestep [fs] — **must match** nuclear `VelocityVerlet` dt |

Compose with `cpu::md::MD` + `VelocityVerlet` at `md::Phase::force`.

## Algorithm (foothold)

**Electronic Lagrangian (closed-shell RKS, occupied \(C \in \mathbb{R}^{n_b\times n_{\mathrm{occ}}}\)):**

\[
\mathcal{L} = \sum_I \tfrac12 M_I \dot R_I^2
            + \sum_i \tfrac12 \mu \lVert \dot c_i\rVert^2
            - E_{\mathrm{KS}}[C,R]
            - \mathrm{constraints}\ (C^\top S C = I).
\]

**EOM (unconstrained force matches OT band-energy gradient factor for \(P=2CC^\top\)):**

\[
\mu\,\ddot C = -4\,H[P]\,C \quad(+\;S\text{-orthonormality via Löwdin after Verlet}).
\]

**Position Verlet + Löwdin (SHAKE-free foothold):**

\[
C_{\mathrm{raw}} = 2C_n - C_{n-1} + \Delta t_{\mathrm{au}}^2\,\ddot C,\quad
C_{n+1} = C_{\mathrm{raw}}\,(C_{\mathrm{raw}}^\top S C_{\mathrm{raw}})^{-1/2}.
\]

**Nuclear forces** at fixed \(C\): analytic GPW force assembly with
\(P=2CC^\top\) and energy-weighted \(W = 2\,C\,(C^\top H C)\,C^\top\)
(constraint contribution off the BO surface).

**Per MD force launch:**

1. **Seed step (once):** full SCF → analytic BO forces → `Cpmd::seed(C_occ)`.
2. **Later steps:** fixed-density KS eval at \(P(C)\) (no SCF) → forces at
   that state → electronic Verlet advance → write \(E_{\mathrm{KS}}\) and
   \(F_I\) to Frame.

**Conserved quantity (eV):**

\[
E_{\mathrm{cons}} = E_{\mathrm{kin}}^{\mathrm{nuc}}
                  + E_{\mathrm{KS}}
                  + T_{\mathrm{e}},\quad
T_{\mathrm{e}} = \tfrac12\mu\sum_{i}\lVert\dot c_i\rVert^2.
\]

## Scope / non-goals

| In | Out |
| -- | --- |
| Closed-shell RKS Γ GPW, PBE, MT, DZVP water | UKS, OT, k-points |
| Fictitious \(\mu\), electronic KE, CP \(E_{\mathrm{cons}}\) | Electronic thermostat |
| Löwdin reorth (foothold) | Full RATTLE/SHAKE constraint forces |
| Short H₂O NVE envelope | Production long traj / CP2K energy twin |
| Refuse UKS/OT/xl_bomd loudly | Silent fallback to BOMD |

## Acceptance (summary)

See `acceptance.md`. Unit: propagator + evaluator guards. E2E:
`regressions/dynamics/cpmd-water` — CP \(E_{\mathrm{cons}}\) drift envelope,
finite trajectory, molecules intact, **not** full-SCF-per-step after seed.

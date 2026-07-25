# Car–Parrinello MD — H₂O NVE (GPW fictitious electronic DOF)

*(Technical-report draft; narrative house style TBD — rewrite later.)*

| Field | Value |
| ----- | ----- |
| Slug | `cpmd-water` |
| Regression | `regressions/dynamics/cpmd-water/` |
| Public API | `cpu::md::MD` + `VelocityVerlet` + `cpu::gpw::GPWEvaluator{…, compute_forces=true, cpmd=true, cpmd_mu, cpmd_dt}` at `md::Phase::force`; propagator `cpu::gpw::Cpmd` |
| Roadmap / queue | [cp2k-method-queue.md](../../developer/cp2k-method-queue.md) row `cpmd-water`; MD roadmap Car–Parrinello |
| Dual gate | Correctness A1–A15 of the spec acceptance; speed **`perf_exempt`** |
| Wave | CP2K capability loop — wave 3 |

## Summary

**Car–Parrinello molecular dynamics** (Car & Parrinello, *PRL* **55**, 2471
(1985)) treats the occupied Kohn–Sham orbitals as classical fictitious
degrees of freedom with mass \(\mu\). The extended Lagrangian is

\[
\mathcal{L} = \sum_I \tfrac12 M_I \dot R_I^2
            + \sum_i \tfrac12 \mu \lVert \dot c_i\rVert^2
            - E_{\mathrm{KS}}[C,R]
            - \text{constraints }(C^\top S C = I).
\]

Nuclei integrate with the ordinary MD `VelocityVerlet`. Electrons integrate
with a position-Verlet step on the occupied MO matrix \(C\)
(\(\mu\,\ddot C = -4 H C\) for closed-shell \(P = 2 C C^\top\)), followed by
Löwdin re-orthonormalization (SHAKE-free foothold). After one cold SCF seed,
**production steps do not re-converge SCF** — they evaluate the KS energy and
matrix once at the current density and take analytic nuclear forces at that
fixed electronic state.

This is **not** Born–Oppenheimer MD (BOMD: SCF every step) and **not**
XL-BOMD (auxiliary density + SCF). The literature conserved quantity is

\[
E_{\mathrm{cons}} = T_{\mathrm{n}} + T_{\mathrm{e}} + E_{\mathrm{KS}},
\]

with \(T_{\mathrm{e}} = \tfrac12\mu\sum_i\lVert\dot c_i\rVert^2\) exposed as
`last_electronic_ke_ev()` [eV] and \(E_{\mathrm{KS}}\) on
`frame.meta.get(kEnergy)`.

Scope: closed-shell RKS Γ-point GPW. UKS, OT, and simultaneous `xl_bomd` are
refused loudly at `build()`.

## Usage

```cpp
constexpr double dt = 0.1; // fs — nuclear and electronic share one step
cpu::gpw::GPWEvaluator eval(
    {O, H}, /*cutoff_ry*/ 200.0, /*max_iter*/ 60, /*tol*/ 1e-5,
    /*mix*/ 0.3, cpu::gpw::XcKind::PBE, /*nonperiodic*/ true,
    /*compute_forces*/ true);
eval.cpmd = true;
eval.cpmd_mu = 2000.0;  // a.u.; AO residual is stiffer than PW textbook ~400
eval.cpmd_dt = dt;      // MUST match VelocityVerlet dt

cpu::md::MD sim;
sim.set_integrator(cpu::md::VelocityVerlet{dt});
sim.on(md::Phase::force).add_ref(eval, /*order=*/0);
sim.build();
sim.run(frame, n_steps);

// CP conserved energy [eV]:
const double e_cons = frame.meta.get(kKinetic)
                    + frame.meta.get(kEnergy)
                    + eval.last_electronic_ke_ev();
```

The propagator primitive (`cpu/gpw/cpmd.h`):

```cpp
cpu::gpw::Cpmd cp(2000.0);  // throws if mu <= 0
cp.seed(C_occ);             // after cold SCF
cp.align_to_overlap(S);     // after nuclei move
cp.advance(H, S, dt_fs);    // projected Verlet + Löwdin
double te_ha = cp.kinetic_ha();
```

## Pin

- **Unit:** mass/seed/advance guards; `GPWEvaluator` refuses UKS / OT /
  `xl_bomd` / non-positive \(\mu,\Delta t\) (`test_cpu_gpw_cpmd`).
- **E2E:** `reg_cpmd-water` — H₂O NVE twin family (PBE / DZVP-MOLOPT-GTH /
  GTH-PBE / MT, 200 Ry, dt 0.1 fs, \(\mu=2000\)): relative \(E_{\mathrm{cons}}\)
  drift ~2×10⁻³ (tol 0.10) after Te is engaged; production `iters==0`;
  bounded O–H; wall ~32 s. See
  `regressions/dynamics/cpmd-water/REFERENCE.md`.

## Honest non-goals

- Full RATTLE constraint forces (Löwdin reorth only)
- Electronic thermostat / NVT-CP
- UKS, OT seed, k-points
- CP2K binary energy twin (literature envelope + structural gates;
  `perf_exempt`)
- Second-generation CP (`second-gen-cpmd-water`) — separate queue row

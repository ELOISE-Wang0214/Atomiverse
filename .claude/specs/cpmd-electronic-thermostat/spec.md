# Spec — `cpmd-electronic-thermostat`

**Slug:** `cpmd-electronic-thermostat`  
**Domain:** dynamics (GPW / AIMD)  
**Status:** in progress  
**Branch:** `wf/cpmd-electronic-thermostat`  
**Parent:** `dev` (after CPMD + RATTLE + UKS product; residual electronic \(T\) control)

## Intent

Add a **Nosé–Hoover electronic thermostat** on the fictitious Car–Parrinello
degrees of freedom so \(T_e = \tfrac12\mu\lVert\dot C\rVert_F^2\) is regulated
toward a target kinetic energy \(K_0\). This is the classical CP remedy for
energy flow into the electronic subsystem (adiabaticity restoration;
Blöchl & Parrinello; Sprik; Marx & Hutter §2.4). Nuclei remain NVE
`VelocityVerlet` in the textbook pin; only the electronic DOF are thermostatted.

## Literature

- Car & Parrinello, *Phys. Rev. Lett.* **55**, 2471 (1985)
- Andersen, *J. Comput. Phys.* **52**, 24 (1983) — electronic RATTLE (already product)
- Blöchl & Parrinello, *Phys. Rev. B* **45**, 9413 (1992) — separate ionic/electronic thermostats
- Sprik, *J. Phys. Chem.* **95**, 2283 (1991) — electronic thermostat for adiabaticity
- Marx & Hutter, *Ab Initio Molecular Dynamics* (CUP) — NH on electrons; extended
  conserved energy; metallic gap caveat (BOMD preferred for metals)
- Twin family: `cpmd-rattle-water`, `cpmd-water`, `cpmd-uks`

## Public API (no Config bags / factories)

| Symbol | Role |
| ------ | ---- |
| `cpu::gpw::Cpmd` | Electronic VV + RATTLE + optional NH thermostat (members) |
| `Cpmd::{thermostat, ekin_target_ha, tau_e_fs}` | Enable + target \(K_0\) [Ha] + \(\tau_e\) [fs] |
| `Cpmd::thermostat_energy_ha()` | \(\tfrac12 Q\xi^2 + 2 K_0\eta\) [Ha] for \(E_{\mathrm{ext}}\) |
| `Cpmd::xi()` / `eta()` | NH friction / coordinate diagnostics |
| `Cpmd::seed(C)` / `seed(C, Cdot)` | Rest seed or seeded electronic velocity |
| `GPWEvaluator::{cpmd_thermostat, cpmd_ekin_target_ha, cpmd_tau_e}` | Engagement (members) |
| `GPWEvaluator::last_electronic_thermostat_ev()` | Last \(E_{\mathrm{therm}}\) [eV] |
| `GPWEvaluator::last_cp_energy_ev()` | Still \(E_{\mathrm{KS}}+T_e\); caller adds \(E_{\mathrm{therm}}\) for NH |

Compose with `cpu::md::MD` + `VelocityVerlet` at `md::Phase::force` (same as
`cpmd-rattle-water`).

## Algorithm

**Electronic force + RATTLE** unchanged (see `cpmd-rattle-water`).

**Single-chain Nosé–Hoover** on \(\dot C\) (target kinetic energy \(K_0 > 0\)):

\[
Q = 2 K_0\,\tau_e^2,\qquad
\dot\xi = \frac{2 T_e - 2 K_0}{Q},\qquad
\dot\eta = \xi.
\]

Symmetric half-steps around mechanical VV+RATTLE (atomic units of time):

1. \(\xi \leftarrow \xi + (\Delta t/2)\,(2T_e-2K_0)/Q\)
2. \(\dot C \leftarrow \dot C\,e^{-\xi\,\Delta t/2}\); \(\eta \leftarrow \eta + \xi\,\Delta t/2\)
3. Mechanical VV half-vel → pos → position RATTLE → complete vel → velocity RATTLE
4. Recompute \(T_e\); repeat (1)–(2) for the second half-step
5. \(T_e = \tfrac12\mu\lVert\dot C\rVert_F^2\);
   \(E_{\mathrm{therm}} = \tfrac12 Q\xi^2 + 2 K_0\eta\)

**Extended conserved energy** (NVE nuclei + electronic NH):

\[
E_{\mathrm{ext}} = T_{\mathrm{n}} + T_{\mathrm{e}} + E_{\mathrm{KS}} + E_{\mathrm{therm}}.
\]

When `thermostat == false`, \(\xi=\eta=0\), \(E_{\mathrm{therm}}=0\), recovers NVE CP.

**UKS:** independent NH state per spin channel (`cpmd_` / `cpmd_beta_`); same
\(K_0,\tau_e\) from evaluator members; \(T_e\) and \(E_{\mathrm{therm}}\) sum both.

## Scope / non-goals

| In | Out |
| -- | --- |
| Closed-shell RKS Γ GPW water + NH electronic \(T\) | Nuclear NVT as the pin (allowed as compose, not required) |
| RATTLE + electronic NH | Langevin / CSVR electronic (later depth) |
| Extended \(E_{\mathrm{ext}}\) drift envelope | Live CP2K in ctest |
| Unit gates on regulation + NH energy | Multi-k CPMD, OT+CPMD, metallic CP |
| Offline CP2K twin input (BOMD-style NVT note) | Nuclear electronic joint thermostat chains |

## Acceptance (summary)

See `acceptance.md`. Unit: thermostat regulates \(T_e\), RATTLE still holds,
invalid params throw. E2E: `regressions/dynamics/cpmd-electronic-thermostat` —
\(E_{\mathrm{ext}}\) drift, \(T_e\) ceiling, `iters==0`, molecule intact.

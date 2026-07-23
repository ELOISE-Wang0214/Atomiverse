# nph_lj — isotropic NPH enthalpy conservation (LJ)

## What is checked

Thermostat-free isotropic NPH integrator on 64 LJ atoms. After warm-up,
the extended enthalpy

\[
H = E_\mathrm{kin} + E_\mathrm{pot} + P V + \tfrac12 W\dot{\varepsilon}^2
\]

(published on `frame.meta[kConserved]`) must keep relative drift
\(|\Delta H|/|H_0|\) below a hard-coded gate.

| Quantity | Gate | Rationale |
| -------- | ---- | --------- |
| \(\max\|\Delta H\|/\|H_0\|\) | \(\le 0.25\) | small-N Trotter oscillation after SC warm-up; fails on secular leak |
| volume | positive finite | barostat stability |

## Parameters

- Same LJ argon set as `lj_argon` / `pr_npt_lj`
- \(P_\mathrm{target}=0\), \(\tau_P=500\,\mathrm{fs}\), \(dt=2\,\mathrm{fs}\)
- Warm-up 200 steps before measuring \(H_0\)
- No temperature target (NPH)

## Literature

Martyna, Tobias, Klein, *J. Chem. Phys.* **101**, 4177 (1994) — MTK
barostat without atom thermostat yields conserved enthalpy.

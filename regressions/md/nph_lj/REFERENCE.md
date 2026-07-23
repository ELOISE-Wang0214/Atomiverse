# NPH LJ argon — literature / gate basis

## What this regression pins

Public thermostat-free `cpu::md::NphIntegrator` (isotropic NPH):

1. Extended enthalpy
   \(H = E_\mathrm{kin}+E_\mathrm{pot}+P V+\tfrac12 W\dot\varepsilon^2\)
   written to `kConserved` stays within a relative drift gate over production.
2. Volume remains positive and non-runaway.

## Physical setup

| Quantity | Value | Note |
| -------- | ----- | ---- |
| Potential | LJ argon-like | Same \(\varepsilon,\sigma\) as PR regression |
| Lattice | Compressed FCC | Dense start |
| Ensemble | NPH (no atom thermostat) | |
| \(P\) | 0 eV/Å³ | Hydrostatic |
| \(\tau_P\) | 1000 fs | Barostat mass only |
| \(dt\) | 1 fs | |

## References

- Martyna, Tobias, Klein, *J. Chem. Phys.* **101**, 4177 (1994) — MTK barostat;
  NPH is the thermostat-free subset (conserved extended enthalpy).

Relative drift bound is an \(O(dt^2)\) engineering envelope for short production.

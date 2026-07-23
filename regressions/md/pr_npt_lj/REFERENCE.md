# pr_npt_lj — Parrinello–Rahman NPT on LJ argon

## What is checked

Short flexible-cell NPT run of a 64-atom LJ fluid (Rahman-like
parameters) with hydrostatic \(P_\mathrm{target}=0\).

Hard-coded gates (stability / finite-sample windows, not EOS pins):

| Quantity | Gate | Rationale |
| -------- | ---- | --------- |
| mean density | 0.30–1.50 × start | SC lattice expands under P=0 (same as MTK); collapse/explosion fails |
| \(\lvert\langle P\rangle\rvert\) | \(\le 5\times 10^{-3}\) eV/Å³ | soft bound ~ few × \(\rho k_B T\) |
| cell shear / edge | < 0.05 | cubic hydrostatic → off-diagonal \(h\) stays small |
| volume / conserved | positive + finite | integrator stability |

## Parameters

- \(\varepsilon/k_B = 120\,\mathrm{K}\), \(\sigma = 3.4\,\mathrm{Å}\), \(r_c=2.5\sigma\)
- \(\rho_0 = 0.02071\,\mathrm{Å}^{-3}\), \(T=94.4\,\mathrm{K}\)
- \(\tau_T=50\,\mathrm{fs}\), \(\tau_P=500\,\mathrm{fs}\), \(dt=2\,\mathrm{fs}\)

## Literature context

Parrinello & Rahman, *J. Appl. Phys.* **52**, 7182 (1981);
Martyna, Tobias, Klein, *J. Chem. Phys.* **101**, 4177 (1994).
Gates document stable public Integrator dynamics, not a tabulated EOS point.

# Spec — `xl-bomd-linearized-forces`

**Slug:** `xl-bomd-linearized-forces`  
**Loop:** CP2K capability (GPW / periodic DFT)  
**Branch:** `wf/xl-bomd-linearized-forces`  
**Method:** Zero-SCF linearized-functional XL-BOMD forces (Niklasson shadow
Born–Oppenheimer potential)

## Problem

SCF-converged XL-BOMD is product (`xl-bomd-water`, `xl-bomd-uks`, multi-k
variants): the auxiliary density \(P\) is only an SCF **guess**, and every MD
step still iterates to `tol`. The Niklasson residual is the **zero-SCF** form
that evaluates forces of the **linearized / shadow** energy functional around
the auxiliary density — one KS construction + one diagonalization per step,
no iterative SCF.

## Science

- A. M. N. Niklasson, *Phys. Rev. Lett.* **100**, 123004 (2008) — XL-BOMD.
- A. M. N. Niklasson et al., *J. Chem. Phys.* **130**, 214109 (2009) —
  dissipative Table-I coefficients (already unit-pinned).
- A. M. N. Niklasson & M. J. Cawkwell, *Phys. Rev. B* **86**, 174308 (2012);
  Niklasson, *J. Chem. Phys.* **147**, 054103 (2017) — shadow / linearized
  free-energy functional; zero iterative SCF after seed.
- Shadow energy (closed-shell KS, linearized around auxiliary \(n\)):

  \[
  \mathcal{E}[D,n] = \mathrm{Tr}(H[n]\,D) - \mathrm{Tr}(H[n]\,n)
                   + E_{\mathrm{KS}}[n]
  \]

  with \(D = q[n]\) the idempotent density from **one** diagonalization of
  \(H[n]\). At the fixed point \(D = n\), \(\mathcal{E} = E_{\mathrm{KS}}\).
- Nuclear forces: Hellmann–Feynman of \(\mathcal{U}(R,n)\) at fixed \(n\),
  assembled with density \(D\) and KS matrix \(H[n]\)
  (`gpw_forces_at_density`).

## Public API (no Config bag)

```cpp
cpu::gpw::GPWEvaluator eval({O, H}, /*cutoff_ry*/200, /*max_iter*/60,
                            /*tol*/1e-5, /*mix*/0.3, XcKind::PBE,
                            /*nonperiodic*/true, /*compute_forces*/true);
eval.xl_bomd = true;
eval.xl_bomd_k = 5;
eval.xl_bomd_linearized = true;   // NEW — zero-SCF shadow forces after seed
```

### Semantics

| Flag combo | Behaviour |
| ---------- | --------- |
| `xl_bomd` only | Existing SCF-converged XL (guess = aux, full SCF). |
| `xl_bomd` + `xl_bomd_linearized` | Cold seed: full SCF. Warm steps: \(H[n]\) + one diag → \(D\); shadow energy + forces at \((D,H[n])\); XL advance(\(D\)). `last_result().iters == 0`. |
| `xl_bomd_linearized` without `xl_bomd` | `build()` throws. |

### Scope (foothold)

- Closed-shell RKS, Γ-point, diagonalization path.
- Requires `compute_forces = true`.
- Refuses loudly (no silent degrade): UKS, multi-k, `use_ot`, `cpmd`,
  `second_gen_cp`.

## Out of scope

- Kernel / Krylov preconditioned electronic EOM (Niklasson 2020+).
- UKS / multi-k linearized XL.
- OT + XL.
- Live CP2K in ctest.

## Deliverables

1. Spec (this file) + `acceptance.md`.
2. `GPWEvaluator::xl_bomd_linearized` + launch path.
3. Unit tests in `test_gpw_xl_bomd.cpp` (A1–A4).
4. Regression `regressions/dynamics/xl-bomd-linearized-forces/` +
   `REFERENCE.md` (literature / offline CP2K twin when available).
5. Oracle sketch `tools/cp2k_oracle/inputs/xl-bomd-linearized-forces/`.

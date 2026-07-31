# Spec — `xl-bomd-linearized-uks`

**Slug:** `xl-bomd-linearized-uks`  
**Loop:** CP2K capability (GPW / periodic DFT)  
**Branch:** `wf/xl-bomd-linearized-uks`  
**Method:** Zero-SCF linearized-functional XL-BOMD forces on Γ-UKS
(dual-spin shadow Born–Oppenheimer potential)

## Problem

Zero-SCF linearized XL-BOMD is product on the Γ RKS diagonalization path
(`xl-bomd-linearized-forces`): after the cold full-SCF seed, each warm step
builds \(H[n]\) once from the auxiliary density, diagonalizes once, and
evaluates shadow forces with **zero iterative SCF**. The evaluator throws
for UKS at `gpw_evaluator.cpp` (`"Γ RKS only (UKS linearized XL is a
follow-on)"`) — the declared wave-14 follow-on. SCF-converged XL on Γ-UKS
is already product (`xl-bomd-uks`, dual α/β `XlBomd` propagators), so the
propagator plumbing exists; only the linearized warm step is RKS-only.

## Science

Same references as `xl-bomd-linearized-forces`, spin-resolved:

- A. M. N. Niklasson, *Phys. Rev. Lett.* **100**, 123004 (2008) — XL-BOMD.
- A. M. N. Niklasson et al., *J. Chem. Phys.* **130**, 214109 (2009) —
  dissipative Table-I coefficients (unit-pinned; identical per spin channel).
- A. M. N. Niklasson & M. J. Cawkwell, *Phys. Rev. B* **86**, 174308 (2012);
  Niklasson, *J. Chem. Phys.* **147**, 054103 (2017) — shadow / linearized
  functional. The linearization is spin-diagonal for collinear UKS: each
  channel σ ∈ {α, β} propagates its own auxiliary density \(n^\sigma\), and

  \[
  \mathcal{E}[D^\alpha, D^\beta, n^\alpha, n^\beta]
    = \sum_\sigma \left[\mathrm{Tr}(H^\sigma[n]\,D^\sigma)
                        - \mathrm{Tr}(H^\sigma[n]\,n^\sigma)\right]
      + E_{\mathrm{KS}}[n^\alpha, n^\beta]
  \]

  with \(D^\sigma = q^\sigma[n]\) from **one** generalized diagonalization
  of \(H^\sigma[n]\) each (occupancy 1). At the fixed point
  \(D^\sigma = n^\sigma\), \(\mathcal{E} = E_{\mathrm{KS}}\).
- Forces: Hellmann–Feynman of \(\mathcal{U}(R, n^\alpha, n^\beta)\) at fixed
  \(n^\sigma\), assembled with densities \(D^\sigma\) and KS matrices
  \(H^\sigma[n]\) — the existing dual-spin fixed-density overload of
  `gpw_forces_at_density` (CPMD-UKS path).

## Implementation

All pieces exist; the work is wiring the UKS variant of the linearized warm
branch in `GPWEvaluator::launch` and dropping the build() UKS throw:

1. `build()`: remove the `n_alpha + n_beta > 0` refusal under
   `xl_bomd_linearized` (keep multi-k / OT / CPMD / SGCP refusals). Dual
   `xl_` / `xl_beta_` propagators are already emplaced for UKS.
2. `launch()` linearized branch: for UKS require both propagators ready,
   run the fixed-density dual-KS eval (`max_iter = 0`,
   `guess_P = n^α`, `guess_P_beta = n^β` — existing solver path), then:
   one generalized diagonalization per spin → \(D^\sigma\) (occ 1),
   dual-spin shadow energy as above, Löwdin ⟨S²⟩ at \(D^\sigma\),
   `iters = 0` / `converged = true`, `advance(D^σ)` on both propagators,
   forces via the dual-spin `gpw_forces_at_density` overload.
3. Cold seed unchanged: full SCF falls through to the SCF-XL path which
   already advances both spin propagators from \(P^\sigma_{\mathrm{SCF}}\).

## Public API (no Config bag — unchanged shape)

```cpp
cpu::gpw::GPWEvaluator gpw({O}, /*cutoff_ry*/200, /*max_iter*/80,
                           /*tol*/1e-5, /*mix*/0.4, XcKind::PADE,
                           /*nonperiodic*/true, /*compute_forces*/true);
gpw.n_alpha = 7;   // O2 triplet (12 valence e-)
gpw.n_beta = 5;
cpu::electronic::XlBomdForce eval{std::move(gpw), 5,
                                  cpu::electronic::XlBomdLinearized{}};
```

### Regression system (decided during impl)

The queue sketched an OH-radical "water-family" twin; measured on the worker
host the OH doublet cold seed is a pathological DIIS case (converges only at
~196 iterations regardless of `mix`/`tol`, PBE; PADE worse at 300). The
regression instead uses the **identical twin of `xl-bomd-uks`** (O₂ triplet,
PADE, seed ≈ 23 iterations) with only the `XlBomdLinearized{}` tag added —
a direct SCF-XL vs zero-SCF linearized-XL A/B on the same open-shell system.
The OH finding is recorded in the regression REFERENCE.

### Semantics

| Case | Behaviour |
| ---- | --------- |
| Γ RKS + linearized | unchanged (existing product) |
| Γ UKS + linearized | **NEW** — cold seed full SCF; warm steps dual-spin zero-SCF shadow forces, `last_result().iters == 0` |
| multi-k + linearized | `build()` throws (follow-on `xl-bomd-linearized-kpoints`) |
| OT / CPMD / SGCP + linearized | `build()` throws (unchanged) |

## Out of scope

- Multi-k linearized XL (`xl-bomd-linearized-kpoints`, held sibling row).
- Kernel / Krylov preconditioned electronic EOM (Niklasson 2020+).
- OT + XL.
- Live CP2K in ctest.

## Deliverables

1. Spec (this file) + `acceptance.md`.
2. UKS linearized warm branch in `src/cpu/gpw/gpw_evaluator.cpp` +
   guard/doc updates (`gpw_evaluator.h`, `xl_bomd_force.h`).
3. Unit tests in `tests/cpu/gpw/test_gpw_xl_bomd.cpp` (A1–A3).
4. Regression `regressions/dynamics/xl-bomd-linearized-uks/` +
   `REFERENCE.md` (literature envelope; offline CP2K twin when a binary
   is available).
5. Oracle sketch `tools/cp2k_oracle/inputs/xl-bomd-linearized-uks/`.

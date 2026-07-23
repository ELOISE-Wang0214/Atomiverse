# Reference — electrode charge predictor (ASPC, E2)

oracle: metalwalls  
oracle_version: algorithm MetalWalls const-pot BO/CG + Kolafa ASPC (J. Comput. Chem. 25, 335 (2004); Coretti et al., J. Chem. Phys. 157, 184801 (2022); Marin-Laflèche et al., JOSS MetalWalls)

## Physics

Constant-potential electrode charges \(Q(t)\) evolve smoothly as the electrolyte
moves. MetalWalls warm-starts the iterative (CG) charge solve from a short
history of prior \(Q\) via Kolafa’s Always-Stable Predictor-Corrector (ASPC).
LAMMPS ELECTRODE has **no** ASPC twin → primary offline oracle is **metalwalls**
(algorithm pin; no live MW binary in ctest).

## ASPC pin (Kolafa \(K=2\))

Lagrange coefficients for predicting at step \(n+1\) from history
\(\{Q_{n-2}, Q_{n-1}, Q_n\}\) (oldest → newest):

\[
\omega = \{1,\ -3,\ 3\}
\]

so \(Q^{\mathrm{pred}}_{n+1} = Q_{n-2} - 3 Q_{n-1} + 3 Q_n\).

Linear history \(Q_k = k\) predicts \(Q_{n+1} = n+1\) exactly (order-2 exact for
polynomials of degree ≤ 2).

## Geometry

### CG vs dense pin (2 electrode sites)

Madelung-like 4-atom primary cell (same family as
`const-pot-bo-parallel-plate`, `xy_image_range=0`):

| Site | Type | role | \(z\) [Å] |
| --- | --- | --- | --- |
| 0 | 1 | bot electrode | 8 |
| 1 | 2 | top electrode | 12 |
| 2 | 3 | electrolyte + | 9 |
| 3 | 3 | electrolyte − | 11 |

Box \(L_x=L_y=1\), \(L_z=20\) Å. \(\eta=2\) Å⁻¹, \(\Psi=(0,1)\) V.

### CG iteration comparison (multi-site)

3×3 atoms per wall (18 electrode sites) + 4 electrolyte charges in
\(L_{xy}=6\) Å, \(L_z=20\) Å. Electrolyte \(z\) drifts each launch so
\(b=\Phi_{\mathrm{elec}}\) changes. Primary cell only.

## What Atomiverse asserts

Public API:

- `cpu::electrodes::ChargeASPC` (vector ASPC)
- `cpu::md::ConstPotBO` with `use_cg`, `predictor_order`, `last_iters()`

| Check | Tol / gate |
| --- | --- |
| ASPC \(K=2\) coeffs \(\{1,-3,3\}\) | exact |
| CG \(Q\) vs dense Cholesky on fixed RHS | abs ≤ \(10^{-8}\) e |
| Mean CG iters with `predictor_order=2` **&lt;** cold (`predictor_order=0`) after history fill | strict inequality (observed cold≈6, warm≈5 on 18-site plate) |
| \(\lvert\sum Q_i\rvert\) | ≤ \(10^{-8}\) e |

## Perf

```
perf_exempt: true
perf_exempt_reason: no MetalWalls twin binary on PATH for fair wall; ASPC is MW-only algorithm pin
ratio_atv_over_oracle: n/a
within_10x: exempt
```

## Repro of pin (offline only — never ctest)

Kolafa coefficients are analytic (Eqs. 7–11 of Kolafa 2004). MetalWalls
documents charge extrapolation for CG/BO electrode solves; re-run a MW
const-pot demo with/without predictor offline if a binary is available and
update iteration counts here — not required for this foothold.

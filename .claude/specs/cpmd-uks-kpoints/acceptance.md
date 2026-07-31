# Acceptance — `cpmd-uks-kpoints`

## Correctness

| ID | Gate | Where |
| -- | ---- | ----- |
| A1 | `cpmd` + multi-k + UKS (`n_alpha`/`n_beta`) builds (no throw); multi-k CPMD **thermostat** still throws at `build` | unit |
| A2 | Multi-k UKS fixed-density eval (`max_iter == 0`, dual k-guess) exports `H_ks_k` + `H_ks_k_beta` + `S_k` + `P_k` + `P_k_beta` (size `n_kpts`), `iters == 0`, energy ≈ cold SCF at the same density; missing β guess throws | unit |
| A3 | Γ-limit identity: k-path dual-spin fixed-density eval at one explicit Γ point reproduces the Γ dual-spin fixed-density eval (energy and \(H^\sigma\)) | unit |
| A4 | E2E periodic H₂ NVE, mesh 1×1×2 unreduced, UKS closed-shell fillings (α+β Bloch CP): cold seed converged; baseline + production `iters == 0`, `n_kpts >= 2`; \(T_e > 0\) once propagated; relative \(E_{\mathrm{cons}}\) drift ≤ 0.10; UKS-vs-RKS closed-shell single-point identity at `tol = 1e-7` ≤ 3e-4 eV; finite coords; 0.4 < r_HH < 1.5 Å | `reg_cpmd-uks-kpoints` |

## Perf

| Field | Value |
| ----- | ----- |
| `perf_exempt` | **true** — oracle family lock: CP2K 2026.1 (present and working on host) has no fictitious-mass CP propagator; the CP conserved quantity with fictitious \(T_e\) is not a number it can print (family policy of `cpmd-kpoints` / `cpmd-uks`, queue-assigned) |
| Else | \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\) |

## Dual gate

Green only when A1–A4 hold and (`within_10x` or `perf_exempt`).

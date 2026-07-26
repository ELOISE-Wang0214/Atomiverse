# Acceptance — ot-uks-dftu-smearing

## Dual gate

| Gate | Requirement |
| ---- | ----------- |
| Unit | `test_cpu_gpw_ot` — multi-k UKS OT+U+FD ≡ DIIS+U+FD free energy |
| Regression | `reg_ot-uks-dftu-smearing` / `run_reg_ot-uks-dftu-smearing` |

## Regression gates

| ID | Check | Tol |
| -- | ----- | ---: |
| A0 | `build` accepts multi-k UKS OT + U + FD | — |
| A1 | OT+U+FD and DIIS+U+FD converge | — |
| A2 | \(\lvert F_{\mathrm{OT}}-F_{\mathrm{DIIS}}\rvert\) | \(10^{-4}\) Ha |
| A3 | \(E_U\) non-trivial; OT vs DIIS \(E_U\) | \(10^{-3}\) Ha |
| A4 | \(N_e = 2\) (valence) | \(10^{-5}\) |
| A5 | \(\lvert F(U)-F(0)\rvert\) activation under OT+FD | \(\ge 10^{-4}\) Ha |
| A6 | absolute vs frozen first-green free energy | \(10^{-3}\) Ha |

## System (regression)

Equal-spin H₂ (DZV-GTH-PADE / GTH-PADE-q1), 8 Bohr cubic, PADE, 100 Ry,
explicit κ = (¼,0,0), `n_alpha=1`/`n_beta=1` seed (shared free-energy
minimum under FD), H \(L=0\) \(U-J=4\) eV, Fermi–Dirac `degauss=0.02` Ry,
`extra_bands=2`, `use_ot=true`.

## Performance

CP2K 2026.1 has no OT+k product twin → **perf_exempt** (internal OT≡DIIS
identity is decisive; same class as `ot-dftu-smearing`).

## Done when

- `acceptance_met = true`
- Spec deleted on merge close (harness)
- Residual prose: multi-k UKS OT+U+smear closed in roadmap/queue

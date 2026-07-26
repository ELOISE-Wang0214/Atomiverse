# Acceptance — `gpw-uks-smearing`

## Gates

| ID | Gate | Pass criterion |
| -- | ---- | -------------- |
| U1 | Unit `test_cpu_gpw_smearing` | UKS+FD `build` accepts (no throw); RKS reject OT+smear still; spinor `find_fermi` recovers Ne |
| R1 | SCF converges | multi-k UKS+FD `converged == true` |
| R2 | Electron count | \(\lvert n_{\mathrm{elec\_occ}} - N_e\rvert < 10^{-5}\) |
| R3 | Fermi level finite | \(\lvert e_{\mathrm{fermi}}\rvert < 10^3\) Ha |
| R4 | Free energy pin | \(\lvert F - F_{\mathrm{pin}}\rvert\) within REFERENCE tol |
| R5 | Closed-shell identity | UKS \(n_\alpha=n_\beta\) + smear ≡ RKS \(n_{\mathrm{occ}}\) + smear within \(5\times10^{-6}\) Ha on H₂ toy (or Na closed shell) |
| P1 | Perf | \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) **or** `perf_exempt` with reason |

## System (regression primary pin)

- BCC Na conventional cubic, 2 atoms, \(a = 4.225\) Å  
- `GTH-PBE-q9` + `DZVP-MOLOPT-SR-GTH`  
- PADE, cutoff 150 Ry, MP \(2\times2\times2\) TR-folded  
- UKS seed `n_alpha=10`, `n_beta=8` (\(N_e=18\); ferro seed)  
- `smearing = fermi_dirac`, `degauss = 0.02` Ry  

Fast identity (optional same binary): high-spin or closed-shell H₂ toy multi-k.

## Offline oracle

Document CP2K UKS+SMEAR input under `tools/cp2k_oracle/inputs/gpw-uks-smearing/`.
Absolute free energy is **ATV freeze** when CP2K wall/grid differs (same class
as `gpw-metallic-smearing`). Run offline CP2K when available for algorithm twin.

## Definition of done

- Spec + unit + regression green under scoped build/test  
- Commit on `wf/gpw-uks-smearing`  
- `acceptance_met = true` when U1 + R1–R5 pass and P1 within_10x or exempt  

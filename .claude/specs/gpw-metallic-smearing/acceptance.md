# Acceptance — `gpw-metallic-smearing`

## Gates

| ID | Gate | Pass criterion |
| -- | ---- | -------------- |
| U1 | Unit `test_cpu_gpw_smearing` | FD closed form + `find_fermi` recovers Ne; public `GpwDft` members compile without Config bag |
| R1 | SCF converges | `GpwResult::converged == true` within `max_iter` |
| R2 | Electron count | \(\lvert n_{\mathrm{elec\_occ}} - N_e\rvert < 10^{-5}\) |
| R3 | Fermi level finite | \(\lvert e_{\mathrm{fermi}}\rvert < 10^3\) Ha |
| R4 | Free energy pin | \(\lvert F - F_{\mathrm{pin}}\rvert\) within REFERENCE tol |
| R5 | Smearing activation | Fixed-occ run (if offered) differs or FD entropy term is nonzero finite when metal |
| P1 | Perf | \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) **or** `perf_exempt` with reason |

## System (regression)

- BCC Na conventional cubic, 2 atoms, \(a = 4.225\) Å  
- `GTH-PBE-q9` + `DZVP-MOLOPT-SR-GTH` (in-tree tables)  
- PADE, cutoff 150 Ry, MP \(2\times2\times2\) TR-folded  
- `smearing = fermi_dirac`, `degauss = 0.02` Ry, `n_occ = 9` (\(N_e=18\))

## Offline oracle

CP2K 2026.1 twin input
`tools/cp2k_oracle/inputs/gpw-metallic-smearing/na_bcc_fermi_dirac.inp`
documents SMEAR + multi-k; absolute free energy stays frozen ATV (grid
details differ). Dual-gate wall: ATV ≈ 1.0 s, CP2K ≈ 89.5 s → within_10x.

## Definition of done

- Spec + regression + unit green under scoped build/test  
- Commit on `wf/gpw-metallic-smearing`  
- `acceptance_met = true` when R1–R4 + U1 pass and P1 within_10x or exempt  

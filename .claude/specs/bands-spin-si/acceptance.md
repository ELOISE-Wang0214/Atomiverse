# Acceptance — bands-spin-si

## Correctness (regression `reg_bands-spin-si`)

System: diamond Si, a = 5.43 Å, 2-atom fcc primitive, ONCVPSP LDA Si UPF
(`tests/data/pw/Si_lda_nc.upf`), ecut = 15 Ry, SCF k-mesh 2×2×2,
`XCType::LSDA`, `spin_polarized = true`, `multiplicity = 1` (\(M=0\)),
`n_bands = 8`, diamond path Γ–X–W–L–Γ–K (`n_seg=4`).

| # | Gate |
| - | ---- |
| 1 | SCF converges (`br.scf_converged`) |
| 2 | `br.nspin == 2`; α and β each have 8 sorted bands per path k |
| 3 | `\|br.magnetization\| ≤ 1e-6` (fixed occ \(N_\uparrow=N_\downarrow=4\)) |
| 4 | Every local same-\(k\) gap (α and β, \(n_{\mathrm{occ}}=4\)) \(> 0\) |
| 5 | Path includes Γ and X; \(\varepsilon_c^\alpha(X) < \varepsilon_c^\alpha(\Gamma)\) |
| 6 | Γ→X gap (α): \(0.05 \le E_{\Gamma\to X} \le 2.50\) eV (Yin–Cohen LDA ~0.45 eV window; same twin settings as closed-shell bands) |
| 7 | α≈β at Γ VBM and X CBM: \(\lvert\varepsilon^\alpha-\varepsilon^\beta\rvert \le 5\times 10^{-4}\) Ry |
| 8 | Public energy still via SCF result [Ry]/[Ha]; Frame not required for bands post-process |

## Dual gate

1. **Correctness:** offline pins above (REFERENCE.md). Never live CP2K in ctest.
2. **Perf:** \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) when a fair twin wall is
   available; else `perf_exempt` + reason in REFERENCE.

## Done when

- Spec files present on `wf/bands-spin-si`
- `Bands::compute` no longer throws on `spin_polarized`
- Unit test binary `test_cpu_pw_bands_dos` green (scoped)
- `reg_bands-spin-si` PASS
- Commit on branch; build wiped; worker JSON reports
  `acceptance_met=true` when gates hold

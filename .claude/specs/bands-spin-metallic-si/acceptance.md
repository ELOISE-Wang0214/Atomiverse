# Acceptance — bands-spin-metallic-si

## Correctness (required)

1. `Bands::compute` with `spin_polarized && smearing != none` does **not**
   throw (prior reject removed).
2. SCF converges on the Si twin; `br.nspin == 2`.
3. α and β path eigenvalues: 8 bands/k, sorted ascending; path includes
   Γ and X.
4. Local same-\(k\) gaps positive for both spins at integer-fill
   \(n_\mathrm{occ}=4\) (Si valence).
5. `gap_ry > 0`; \(E_F\) finite and within \(\pm 4\,\mathrm{degauss}\) of
   \([E_\mathrm{vbm}, E_\mathrm{cbm}]\).
6. \(|M| \le 10^{-3}\) on \(M=0\) seed Si (numerics under FD).
7. α≈β at Γ VBM and X CBM: \(\lvert\varepsilon^\alpha-\varepsilon^\beta\rvert
   \le 10^{-3}\) Ry.
8. Γ→X gap (α) in eV window \([0.05, 2.50]\) (Yin–Cohen LDA class);
   \(\varepsilon_c^\alpha(X) < \varepsilon_c^\alpha(\Gamma)\).
9. Offline pin + REFERENCE.md + CP2K note; no live oracle in ctest.
10. Unit test: spin+smear is accepted (guard-level) or prior reject test
    inverted.

## Performance

- Prefer \(T_\mathrm{atv}/T_\mathrm{cp2k} \le 10\) when a CP2K binary is
  available.
- Else `perf_exempt=true` with reason (no CP2K module on host) — same
  policy as `bands-spin-si` / `bands-metallic-si`.

## Done when

`reg_bands-spin-metallic-si` PASS + unit gate green + commit on
`wf/bands-spin-metallic-si` with `acceptance_met=true`.

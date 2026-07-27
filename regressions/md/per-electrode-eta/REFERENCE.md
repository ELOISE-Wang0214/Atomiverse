# Reference — per-electrode Gaussian width η_e (ConstPotBO multi-wall)

oracle: literature
oracle_version: Reed–Lanning–Madden JCP 126, 084704 (2007) DOI 10.1063/1.2464084; Scalfi et al. JCP 153, 174704 (2020) DOI 10.1063/5.0028232; MetalWalls JOSS 5(53), 2373 (2020) DOI 10.21105/joss.02373

Loop row `per-electrode-eta` is **oracle=metalwalls**; no MetalWalls binary is
installed on the host, so per the `tools/metalwalls_oracle/README.md` contract
the pin is the **literature closed form** with full citation
(`perf_exempt=true`). Would-be MW twin input preserved at
`tools/metalwalls_oracle/inputs/metalwalls/per-electrode-eta/runtime.inpt`;
derivation in `tools/metalwalls_oracle/pin_notes_per-electrode-eta.md`.

## Physics pinned

MetalWalls-style per-electrode Gaussian charge widths on the
constant-potential BO solve (public surface:
`ConstPotBO::add_electrode(type, Psi, eta)`), with the published elastance
elements (units: eV, Å, e; \(k = 14.399645\) eV·Å/e²):

- \(A_{ii} = k\sqrt{2/\pi}\,\eta_i\)
- \(A_{ij} = k\,\mathrm{erf}(\eta_{ij} r_{ij})/r_{ij}\),
  \(\eta_{ij}=\eta_i\eta_j/\sqrt{\eta_i^2+\eta_j^2}\)

## Case 1 — two-site vacuum capacitor, dissimilar η (exact literature pin)

\(\eta_1 = 1.979\) Å⁻¹ (LAMMPS ELECTRODE default, Ahrens-Iwers et al. JCP 157,
084801 (2022)), \(\eta_2 = 0.955\) Å⁻¹ (Merlet/Salanne graphite width used in
MW examples), \(d = 4\) Å, \(\Delta\Psi = 1\) V, electroneutral:

\[
Q = \frac{\Delta\Psi}{A_{11}+A_{22}-2A_{12}} = \mathbf{0.037722100829787}\ e
\]

(\(A_{11}=22.737234510126967\), \(A_{22}=10.972237977347778\),
\(A_{12}=3.599907138217925\), \(\eta_{12}=0.860091344088193\) Å⁻¹).

| Check | Tol |
| --- | --- |
| \(Q_{\mathrm{top}}\) vs 0.037722100829787 e | rel ≤ 1e-9 |
| \(\lvert\sum Q\rvert\) | ≤ 1e-12 e |

## Case 2 — multi-atom two-wall, mixed η, 2D xy Gaussian images

2×2-atom square walls (a = 2 Å, \(L_x=L_y=4\) Å), wall separation 5 Å,
\(\eta_{\mathrm{bot}}=1.979\), \(\eta_{\mathrm{top}}=0.955\) Å⁻¹,
\(\Delta\Psi=1\) V, `xy_image_range = 100`.

Self-consistency ordering (smaller η ⇒ smaller self-elastance ⇒ larger
capacitance, monotone in each diagonal element):

| Check | Tol |
| --- | --- |
| \(\lvert\sum Q\rvert\) | ≤ 1e-10 e |
| \(Q_{\mathrm{top}} > 0\) (with \(\Psi_{\mathrm{top}}>\Psi_{\mathrm{bot}}\)) | sign |
| \(Q^{\mathrm{shared}\ 1.979} < Q^{\mathrm{mixed}} < Q^{\mathrm{shared}\ 0.955}\) | strict |
| 3-arg `add_electrode` at equal η ≡ 2-arg shared path | ≤ 1e-14 e/site |

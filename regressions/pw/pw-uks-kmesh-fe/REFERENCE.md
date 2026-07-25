# Reference — k-point collinear UKS (O₂ multi-k)

## Method

Collinear spin-polarized plane-wave KS (`cpu::pw::PW` with
`spin_polarized = true`, `multiplicity = 3`) on a **Monkhorst–Pack
`k_mesh = 2`** mesh (not Γ-only). Fixed occupations N↑ = 7, N↓ = 5 for 12
valence electrons (oxygen NC-PP Z_val = 6 × 2). Closes the wave-2 ladder
step after Γ UKS (`pw-spin-polarized-o2`) and restricted metallic smearing
(`pw-smearing-al`).

**Geometry choice:** Queue allows Fe or O₂ cell. Bulk Fe metal needs
spin **+** metallic smearing (not co-supported yet on `PW::run`). This pin
uses the O₂ cubic cell so multi-k collinear UKS is gated offline without a
new public API.

## System

| Quantity | Value |
| -------- | ----- |
| Molecule | O₂ |
| R(O–O) | 1.2075 Å |
| Cell | cubic L = 5.5 Å |
| PP | SG15 ONCV PBE oxygen `tests/data/pw/O_ONCV_PBE-1.2.upf` |
| E_cut | 5 Ry |
| k-mesh | Monkhorst–Pack 2×2×2 |
| Multiplicity | 3 (triplet ³Σg⁻) |
| SCF tol (density L2) | 2×10⁻³ |
| max_iter | 20 |
| mix_beta | 0.5 |

## Pins

### Magnetization (literature + fixed occupations)

O₂ ground state is the triplet ³Σg⁻ (Herzberg, *Molecular Spectra and
Molecular Structure* I). With fixed integer occupations and multiplicity 3:

\[
M = N_\uparrow - N_\downarrow = 2.
\]

**Gate:** `|M − 2| ≤ 1×10⁻⁹` on `PWResult::magnetization`.

### Total energy (offline method pin)

\[
E_{\mathrm{tot}} = -25.1087178036~\mathrm{Ha}
\quad (| \Delta E | \le 5\times 10^{-3}~\mathrm{Ha}).
\]

Recorded 2026-07-25 on Arrhenius (x86_64) with the settings above
(converged in 5 SCF iters).

**Re-pin 2026-07-25:** previous frozen value −24.2395602747 Ha (recorded
2026-07-24) was produced by a Davidson subspace step that truncated the
(complex Hermitian) subspace Hamiltonian to its real part — first-order
wrong Ritz vectors at every non-TRIM k — and by G-sphere enumeration
bounds that clipped sphere-edge plane waves. Both fixed with the
`berry-polarization` case (see its REFERENCE.md, "Root-cause fixes"); the
corrected solver lowers this fixed-settings pin. Magnetization (the
literature hard gate) is unchanged.

**Oracle note:** Host CP2K module not used as an absolute-energy twin
(GTH-GPW ≠ NC-UPF PW; same class of limitation as `pw-spin-polarized-o2`).
Absolute e_tot is a **fixed-settings offline method pin** for regression
stability. Magnetization remains the literature hard gate. Units: public
energy on `PWResult` is also reported in Ry (`energy = 2 × energy_ha`).

Optional CP2K documentation twin (algorithm only, not energy-compared):
Quickstep UKS + `KPOINTS` Monkhorst–Pack 2×2×2 on the same O₂ cell — not
required at ctest.

## Perf

| Field | Value |
| ----- | ----- |
| perf_exempt | true |
| reason | CP2K binary not used as fair single-rank PW/NC-UPF twin wall; no loadable comparable twin on pin host |
| T_atv (this host, OMP=8, single rank) | **~31 s** |
| ratio T_atv/T_cp2k | n/a |

## Dual gate

1. Correctness: mag + e_tot gates above; SCF must converge with `k_mesh = 2`.
2. Perf: exempt (no CP2K twin wall).

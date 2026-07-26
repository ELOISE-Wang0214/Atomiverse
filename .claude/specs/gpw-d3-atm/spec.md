# Spec: gpw-d3-atm

## Goal

Add **DFT-D3 three-body Axilrod–Teller–Muto (ATM / C9)** dispersion on the
GPW path, matching offline CP2K `TYPE DFTD3` + `CALCULATE_C9_TERM .TRUE.`.

## Science

Grimme D3 (J. Chem. Phys. 132, 154104, 2010) three-body term as in CP2K
`qs_dispersion_d3.F`:

\[
C_9^{ABC} = -\sqrt{C_6^{AB}\,C_6^{BC}\,C_6^{CA}}
\]
\[
E_9 = s_9\,f_d\,C_9\,\frac{3\cos\theta_a\cos\theta_b\cos\theta_c + 1}{(R_{AB}R_{BC}R_{CA})^3}
\]

with \(s_9 = s_6\), damping \(f_d\) via `damping_d3(rabc, r0, 4/3, 16, rcut)`,
\(r_{abc}=(R_{AB}^2 R_{BC}^2 R_{CA}^2)^{1/6}\),
\(r_0=(R0_{AB}R0_{BC}R0_{CA})^{1/3}\). Angle form uses
`ang = 0.375 s1 s2 s3 / r2abc + 1` (CP2K). CN-dependent \(C_6\) (not
`REFERENCE_C9_TERM`).

## API

- `dftd3_pbe(..., bool calculate_c9 = false)` / `dftd3_pbe_bj(..., bool calculate_c9 = false)`
- `GpwDft::d3_atm` / `GPWEvaluator::d3_atm` — when true and `d3` or `d3bj`,
  include ATM in \(E_\mathrm{disp}\) and analytic forces.

Default remains ATM **off** (two-body product unchanged).

## Pin (offline CP2K 2026.1)

Isolated water dimer, PERIODIC NONE, 12 Å box, PBE-D3 zero + C9:

- **E_disp = −0.00081414118193 Ha**
- Input: `tools/cp2k_oracle/inputs/gpw-d3-atm/water_dimer_d3_atm.inp`
- Two-body-only baseline on same geometry: −0.00081399213205 Ha (ATM Δ ≈ −1.49×10⁻⁷ Ha)

## Acceptance

1. Unit: `dftd3_pbe(..., calculate_c9=true)` matches pin to 1e-9 Ha.
2. Unit: analytic ATM forces match FD to 1e-7 Ha/Bohr.
3. Regression `reg_gpw-d3-atm`: `GpwDft{d3=true,d3_atm=true}` `e_disp` within 1e-9 Ha.
4. Two-body default (`d3_atm=false`) still matches existing D3 pins.
5. Perf: geometry-only E_disp → **perf_exempt** (same as other D3/D4 E_disp slugs).

## Out of scope

- D4 ATM (`gpw-d4-atm` sibling)
- `REFERENCE_C9_TERM` (constant CN C9)
- Long-range C9 LRC
- Molecular `cpu::scf::D3BJ` path

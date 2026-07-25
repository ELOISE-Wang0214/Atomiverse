# Acceptance — `pw-spin-gga`

## Gates

| ID | Check | Tol |
|----|-------|-----|
| A1 | O₂ triplet SCF converges with `spin_polarized` + `xc=PBE` true spin-GGA | — |
| A2 | Magnetization identity M = multiplicity − 1 = 2 | ≤ 1×10⁻⁹ |
| A3 | Finite total energy; \|E − offline method pin\| | ≤ 1×10⁻³ Ha |
| A4 | GGA activation: \|E_PBE_spin − E_LSDA_spin\| > 1×10⁻⁴ Ha (same geometry/settings) | > 1×10⁻⁴ Ha |
| A5 | Unit: closed-shell collapse (ρ↑=ρ↓) spin-GGA vs restricted `step` on nonuniform density | ≤ 1×10⁻⁸ (rel) on E_xc / v |

## Unit

- Extend `test_cpu_pw_spin` or add XC-builder unit proving true spin-GGA is not
  σ=0 for nonuniform spin densities; magnetization identity stays cheap.
- Do **not** put full O₂ e2e in `tests/` (regressions only).

## Perf

`perf_exempt` if no fair CP2K/QE NC-UPF wall twin on pin host (same policy as
`pw-spin-polarized-o2` / `pw-pbe-gga-si`). Target ratio ≤ 10 when twin available.

## Done when

- `XCBuilder::step_spin(ws, xc)` uses spectral ∇ρ_σ for GGA (no production σ=0
  PBE path for spin)
- `ctest -R test_cpu_pw_spin` green
- `run_reg_pw-spin-gga` PASS
- Commit on `wf/pw-spin-gga`

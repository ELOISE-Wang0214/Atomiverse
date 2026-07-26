# Acceptance — `dfpt-metallic-si`

## Checklist

- [x] Spec present under `.claude/specs/dfpt-metallic-si/`
- [x] `DfptPhonon::compute` + `dfpt_gamma_force_constants` accept restricted metallic smearing
- [x] Spin + smear, finite-q + smear still throw
- [x] Unit tests updated (`test_cpu_pw_dfpt`)
- [x] Regression `regressions/response/dfpt-metallic-si/` with REFERENCE.md offline pins
- [x] Scoped build: `atomiverse test_cpu_pw_dfpt reg_dfpt-metallic-si`
- [x] Scoped test green; dual gate (correctness + perf_exempt)
- [x] Commit on `wf/dfpt-metallic-si`
- [ ] Builds wiped (exit cleanup)

## Pins

| Quantity | Pin | Tol |
| -------- | --- | --- |
| Free-energy FD \(C[(0x),(0x)]\) | 17.2393 eV/Å² | 0.30 |
| DFPT \(C_{00}\) | 17.3416 eV/Å² | vs FD pin |
| Optical Γ triplet | 579.45–579.60 cm⁻¹ | [500, 650] |
| Acoustic | ≈ −2 cm⁻¹ | ≤ 30 |

## Perf

`perf_exempt = true` — no CP2K module on host (same as sibling DFPT regs).
ATV wall ≈ 16.5 s.

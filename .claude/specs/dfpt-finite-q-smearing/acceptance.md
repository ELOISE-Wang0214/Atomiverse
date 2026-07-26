# Acceptance — `dfpt-finite-q-smearing`

## Checklist

- [x] Spec present under `.claude/specs/dfpt-finite-q-smearing/`
- [x] `DfptPhonon::compute` accepts restricted metallic smear at commensurate
      finite \(\mathbf{q}\) and direct finite-\(q\) Sternheimer
- [x] `dfpt_q_force_constants` no longer throws on restricted smear
- [x] Spin + smear still throws
- [x] Unit tests updated (`test_cpu_pw_dfpt`)
- [x] Regression `regressions/response/dfpt-finite-q-smearing/` + REFERENCE.md
- [x] Scoped build: `atomiverse test_cpu_pw_dfpt reg_dfpt-finite-q-smearing`
- [x] Scoped test green; dual gate (correctness + perf_exempt)
- [ ] Commit on `wf/dfpt-finite-q-smearing`
- [ ] Builds wiped (exit cleanup)

## Pins (product twin — X + FD)

| Quantity | Pin / window | Tol |
| -------- | ------------ | --- |
| System | diamond Si, a=5.43 Å, ecut 10 Ry, k=1 on SC, FD 0.02 Ry | — |
| \(\mathbf{q}\) | \(X=(1/2,0,0)\) commensurate | — |
| \(\lvert\omega\rvert\) low 3 | 130, 130, 230 | [50, 400] cm⁻¹ |
| \(\lvert\omega\rvert\) high 3 | 448, 673, 673 | [300, 850] cm⁻¹ |
| n_imag | 2 | ≤ 2 |
| asymmetry | 8.4e-8 eV/Å² | ≤ 1e-3 |

## Perf

`perf_exempt = true` — no CP2K module on host (same as sibling DFPT regs).
ATV wall ≈ 40 s.

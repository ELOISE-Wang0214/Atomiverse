# Acceptance — `gpw-fd-phonon-si`

| ID | Gate | Status (2026-07-26) |
| -- | ---- | ------------------- |
| A1 | \(3N = 6\) modes | blocked (reg not green) |
| A2 | ≥3 modes with \(\|\omega\| < 80\) cm⁻¹ (Γ acoustic residual) | blocked |
| A3 | Three non-acoustic modes within **50 cm⁻¹** of CP2K pin (sorted) | blocked |
| A4 | All \(\omega\) finite | blocked |
| A5 | Public `Phonon` + `GPWEvaluator{compute_forces=true}` only | **PASS** (unit wiring) |
| A6 | \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) or `perf_exempt` | n/a |
| A7 | Regression wall ≤ 60 s (Release, 4 OMP) | n/a |

**Block reason:** GPW **analytic forces on Si** (nonperiodic MT) disagree with
energy finite-difference (~271×) and offline CP2K (~40×) while total energy
matches to ~0.014 Ha at L=12 Å. Periodic Γ MOLOPT exceeds the e2e budget.
Water analytic forces remain green. Fix force assembly for Si GTH, then
re-open A1–A4 against the recorded CP2K vib pin.

**Offline CP2K pin (2026-07-26):** CP2K 2026.1, compact diamond fragment
vib optical **−263.481, −194.722, 383.783** cm⁻¹ (wall 10.1 s); L=12
ENERGY_FORCE \(E=-7.585334\) Ha, \(F_x=\pm 3.370\times 10^{-3}\) Ha/Bohr.

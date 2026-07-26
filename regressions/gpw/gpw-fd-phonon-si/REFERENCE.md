# Reference — FD Γ phonons with GPW analytic forces (Si₂)

## Method

Finite-displacement (central difference, \(h = 10^{-3}\) Å) Cartesian force
constants at the **Γ** point of a fixed cell, mass-weighted dynamical matrix,
full \(3N\) spectrum in cm⁻¹. Force engine is **GPW analytic forces** via

```text
cpu::vib::Phonon(1,1,1)
  → vib::Phase::force += cpu::gpw::GPWEvaluator{…, compute_forces=true}
```

No `NeighborList` (GPW is not a `NeighborConsumer`). No public `*Config` /
`make_*`.

## Status (2026-07-26 worker)

| Gate | Result |
| ---- | ------ |
| Unit wiring (`test_cpu_vib_phonon_gpw_wiring`) | **PASS** — `Phonon` accepts `GPWEvaluator` without NeighborList |
| Offline CP2K pin | **recorded** (see below) |
| ATV optical vs CP2K | **BLOCKED** — see diagnosis |
| Perf within 10× | n/a until correctness green |

### Diagnosis (blocking)

1. **Periodic Γ DZVP-MOLOPT** Si₂: single SCF does not finish in ≤ 45 s on
   epyc login (same product constraint as `reg_gpw-pbesol-si`). 12× FD is
   outside the ≤ 60 s e2e budget.
2. **Nonperiodic MT** (ATV `nonperiodic=true`, CP2K `PERIODIC NONE` +
   `PSOLVER MT`), L = 12 Å, Si–Si = 2.351 Å, PADE, CUTOFF 150 Ry, OT:
   - Total energy agrees: ATV \(E \approx -7.571\) Ha vs CP2K
     \(E = -7.585334\) Ha (\(\Delta E \approx 0.014\) Ha).
   - **Analytic force is wrong:** ATV \(F_x \approx 0.131\) Ha/Bohr vs
     CP2K \(F_x = 3.370\times 10^{-3}\) Ha/Bohr (~40×).
   - **Energy finite-difference force** on the same ATV surface:
     \(F_x^{\mathrm{FD}} \approx 4.8\times 10^{-4}\) Ha/Bohr — analytic is
     **~271×** the energy-FD force. Water analytic forces remain green
     (`reg_gpw-forces-water`).
3. Consequence: FD Hessian from analytic forces fails Γ ASR and cannot match
   CP2K `VIBRATIONAL_ANALYSIS` frequencies within 50 cm⁻¹.

**Required follow-up (out of Phonon wiring scope):** fix GPW analytic forces
for Si (nonperiodic MT / multi-valence GTH) so analytic ≡ energy-FD within
~1e-3 Ha/Bohr; then re-pin CP2K vib and re-open this regression gate.

## System (intended twin)

| Setting | Value |
| ------- | ----- |
| Structure | Si₂, bond along \(x\), mid-cell, cubic \(L = 12\) Å preferred |
| Alternate compact | diamond fragment in \(a = 5.43\) Å (tight vacuum; MT quality poor) |
| BC | nonperiodic MT |
| Basis | DZVP-MOLOPT-GTH |
| PP | GTH-PBE-q4 |
| XC | PADE |
| CUTOFF | 150 Ry |
| FD step | \(10^{-3}\) Å |
| Atoms / modes | \(N=2\), \(3N=6\) |

## Offline CP2K pin (recorded)

| Field | Value |
| ----- | ----- |
| Engine | **CP2K 2026.1** (Development `85b8a9b49c`) |
| Host | Arrhenius epyc, offline only |
| Input | `tools/cp2k_oracle/inputs/gpw-fd-phonon-si.inp` |
| RUN_TYPE | `VIBRATIONAL_ANALYSIS` (compact diamond fragment) + ENERGY_FORCE (L=12) |
| Compact vib optical (cm⁻¹) | **−263.481, −194.722, 383.783** (FULLY_PERIODIC T, wall **10.1 s**) |
| L=12 ENERGY_FORCE | \(E=-7.585334\) Ha, \(F_x=\pm 3.370\times 10^{-3}\) Ha/Bohr |

## Acceptance (when unblocked)

| Check | Tol |
| ----- | --- |
| Mode count | exact 6 |
| Acoustic residual | ≥3 modes with \(\|\omega\| < 80\) cm⁻¹ |
| Optical max \(\|\Delta\omega\|\) vs CP2K | **50 cm⁻¹** |
| Wall | ≤ 60 s |
| Perf | \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) |

## Notes

- Sibling SW product: `regressions/pw/fd-phonon-supercell`.
- Sibling PW DFPT: `regressions/response/dfpt-phonon-si`.
- Water force product remains the GPW force foothold (`reg_gpw-forces-water`).

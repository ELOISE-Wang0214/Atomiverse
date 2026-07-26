# Spec — `gpw-d4-atm`

**Slug:** `gpw-d4-atm`  
**Loop:** CP2K method (GPW / periodic DFT)  
**Oracle:** offline **dftd4 3.7.0** (Arrhenius; CP2K 2026.1 lacks DFTD4)  
**Branch:** `wf/gpw-d4-atm`

## Goal

Ship **DFT-D4 three-body Axilrod–Teller–Muto (ATM)** on the GPW product path
with scale \(s_9 = 1\) (dftd4 `pbe-D4-ATM` / `--mbdscale 1.0`), folded into
`GpwResult::e_disp` alongside the existing two-body BJ term.

## Current state (pre)

- Two-body PBE-D4 on GPW is product (`gpw-d4-water`, forces, stress).
- `cpu::scf::dftd4_energy` / `dftd4` document **ATM off** (`s9 = 0`).
- Molecular `cpu::scf::D4` is two-body only.

## Deliverable

1. **ATM energy** in `src/cpu/scf/dftd4.cpp` matching simple-dftd4
   `get_atm_dispersion_energy` (Chai–Head-Gordon zero damping, `alp=16`,
   critical radii from BJ \(a_1,a_2\), C9 from \(\sqrt{|C_6^{ij}C_6^{ik}C_6^{jk}|}\),
   triple-scale for image-self terms). **C6 for ATM use neutral charges**
   (`q = 0` reweight; same as dftd4 `disp.f90`).
2. Public free-function args: `s9` (default **0**, preserves two-body pins)
   and `disp3_cut` (default **40** Bohr). Forces via existing central-diff
   of `dftd4_energy` automatically include ATM when `s9 ≠ 0`.
3. GPW flag `GpwDft::d4_atm` / `GPWEvaluator::d4_atm` (requires `d4`):
   when true, pass `s9 = 1.0` into energy / forces / stress D4 calls.
4. Unit tests: isolated + periodic water-dimer E_disp vs dftd4 pins with
   `s9=1`; two-body regression pins stay green at default `s9=0`.
5. Regression `regressions/gpw/gpw-d4-atm/` hard gate on periodic ATM pin.
6. Offline oracle inputs under `tools/cp2k_oracle/inputs/gpw-d4-atm/`.

## Hard pin (acceptance)

Periodic water dimer (same geometry as `gpw-d4-water`), cubic 12 Å cell:

| Quantity | Value |
| -------- | ----- |
| E_disp (PBE-D4 + ATM, \(s_9=1\)) | **−9.9554475834240288×10⁻⁴ Ha** |
| E_disp (two-body only, \(s_9=0\)) | −9.9458842745846774×10⁻⁴ Ha (context) |
| Source | dftd4 3.7.0, `-f pbe` (default ATM) / `--mbdscale 1.0` |
| Tolerance | **1×10⁻⁸ Ha** |

Isolated (unit test): −9.8819741533958208×10⁻⁴ Ha (`s9=1`).

## Non-goals

- CP2K live QS twin (no DFTD4 on Arrhenius CP2K build).
- Changing default two-body-only product (`d4` without `d4_atm`).
- Molecular post-HF / PySCF-loop work.
- Analytic closed-form ATM forces (FD of energy is sufficient, same as two-body).

## Perf

Geometry-only O(N³ × images) on 6 atoms; full GPW SCF dominates.
**perf_exempt = true** (same policy as `gpw-d4-water`).

## Acceptance checklist

- [x] Spec + offline pin JSON under `tools/cp2k_oracle/inputs/gpw-d4-atm/`
- [x] ATM in `dftd4_energy` when `s9 ≠ 0`; default `s9=0` keeps old pins
- [x] `d4_atm` wired through solver / forces / stress / evaluator
- [x] Unit tests ATM energy (iso + periodic)
- [x] `reg_gpw-d4-atm` PASS within 1e-8 Ha
- [x] Scoped build+test green; commit on `wf/gpw-d4-atm`

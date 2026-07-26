# Spec — gpw-d4-forces

## Goal

Ship **analytic ionic forces** for Γ-point (and multi-k) **GPW + DFT-D4**
(Becke–Johnson, two-body, PBE params; ATM off). Energy (`gpw-d4-water`) and
fixed-P cell stress (`gpw-d4-stress`) already fold geometry-only
\(E_{\mathrm{disp}}^{\mathrm{D4}}\); force assembly currently **skips** D4
(energy-only policy in `gpw_forces.cpp`).

## Method

1. Extend `cpu::scf::dftd4` with a force twin:
   - `D4Result { energy, forces }` where \(F_A = -\partial E_{\mathrm{disp}}/\partial R_A\)
     [Hartree/Bohr].
   - Same energy path as `dftd4_energy` (EN-weighted erf CN, EEQ2019, BJ pair
     sum; periodic Ewald EEQ when requested).
   - Gradient of the geometry-only energy: product path uses a central
     finite-difference of `dftd4_energy` (exact \(\nabla E_{\mathrm{disp}}\) to
     \(O(h^2)\); \(h=10^{-4}\) Bohr). Dual-gates against offline **dftd4 3.7.0
     `-g`** (gradients) and internal energy-FD self-consistency.
2. Wire into `gpw_forces` assembly when `cfg.d4` (mutually exclusive with other
   dispersion flags), same add-on pattern as D3/TS/XDM.
3. No public `*Config` bag; no `make_*` factories. Parameters stay on existing
   `GpwDft::d4` / `GPWEvaluator::d4` members.

## Out of scope

- Three-body ATM / MBD scaling (`s9 ≠ 0`)
- CP2K QS ENERGY_FORCE with `TYPE DFTD4` (Arrhenius CP2K 2026.1 built without
  libdftd4 — same as `gpw-d4-water` / `gpw-d4-stress`)
- Molecular Gaussian post-HF; classical electrochem MD; CUDA

## Acceptance (see acceptance.md)

- Unit: periodic + isolated water-dimer D4 forces match offline dftd4 3.7.0
  gradients (sign: force \(=-\) gradient) within 1e-7 Ha/Bohr; self-FD residual
  within 1e-7.
- Product: `GPWEvaluator` / `gpw_forces_*` with `d4=true` includes D4 contribution
  (delta vs bare electronic forces matches standalone D4 forces).
- Regression dual-gate: analytic product forces vs multi-SCF `gpw_forces_fd`
  with `d4=true` on the periodic water-dimer twin (or documented cheaper
  dual-gate); E_disp pin parity retained.
- perf: \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\) **or** `perf_exempt`
  (CP2K lacks DFTD4 on Arrhenius → dual-gate dominated by multi-SCF FD twin).

## Files (expected)

- `include/atomiverse/cpu/scf/dftd4.h`, `src/cpu/scf/dftd4.cpp`
- `src/cpu/gpw/internal/gpw_forces.cpp`
- `tests/cpu/gpw/test_cpu_gpw_dftd4.cpp` (+ optional forces dual-gate)
- `regressions/gpw/gpw-d4-forces/`
- `tools/cp2k_oracle/inputs/gpw-d4-forces/` (pin notes; reuses d4-water POSCAR)

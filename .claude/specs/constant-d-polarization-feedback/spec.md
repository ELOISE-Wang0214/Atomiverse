# constant-d-polarization-feedback

Wave 5, metalwalls-method-loop. Branch `wf/constant-d-polarization-feedback`.

## Goal

Upgrade the displacement-field foothold `cpu::md::FixDfield` (vacuum
E ≡ D) to the self-consistent **constant-D ensemble** of Zhang & Sprik:
the Maxwell field responds to the instantaneous itinerant polarization,

    E_z(t) = D_z − 4π k_e P_z(t),      P_z = (1/V) Σ_i q_i z_i^unwrapped

with the constant-D electric enthalpy term

    U_D = V/(8π k_e) (D_z − 4π k_e P_z)²   [eV]

so that F_{i,z} = −∂U_D/∂z_i = q_i E_z exactly. `k_e = units::coulomb =
14.399645 eV·Å/e²` (1/(4πε₀) in the eV/Å/e unit system); the 4π k_e factor
converts P [e/Å²] into a field [V/Å].

Scope per queue row: bulk electrolyte, constant **D_z** (feedback along z
only — MetalWalls applies D normal to the cell; x,y keep the vacuum
foothold coupling). Observable: P_z(t), E_z response vs the
vacuum-foothold twin. No new Driver category, no charge solve inside
KSpace — this stays a Bias Component on `cpu::md::MD`.

## Design

- Reuse `FixDfield` (never duplicate): add public member
  `bool polarization_feedback = false`. Default `false` preserves the
  wave-4 vacuum foothold bit-for-bit (existing
  `reg_displacement-dfield-md` must stay green untouched).
- `polarization_feedback = true`:
  - pass 1 computes itinerant P_z from unwrapped coordinates (image
    flags), pass 2 applies F_z += q_i E_z with E_z = D_z − 4π k_e P_z;
    x,y stay F += q D (foothold).
  - energy: U = U_D(z) − Σ q_i (D_x r_x + D_y r_y) added to
    `frame.meta.get(kEnergy)`.
  - throws on triclinic boxes (unchanged) and on V ≤ 0 (field
    ill-defined; no silent fallback).
- New observable accessor `last_Ez()` [V/Å]; `last_Pz()` [e/Å²] kept.

## Acceptance criteria

1. **Foothold lock**: `polarization_feedback = false` (default) reproduces
   the wave-4 behavior exactly; `reg_displacement-dfield-md` green
   without modification; `last_Ez() == Dz` on the foothold path.
2. **Field**: with feedback on, `last_Ez() == Dz − 4π·14.399645·Pz` and
   forces are `q_i · last_Ez()` (z), `q_i · D_{x,y}` (x,y).
3. **Energy/force consistency**: the applied z-force equals
   −∂U/∂z_i of the reported energy (central finite difference, ≤ 1e-8
   relative), i.e. constant-D dynamics conserves KE + U_D.
4. **Itinerant polarization**: P_z uses unwrapped z (image flags shift
   P_z by q·n·L_z/V).
5. **Unit tests** `tests/cpu/dynamics/test_fix_dfield.cpp` registered as
   `test_cpu_fix_dfield`, single-function, seconds.
6. **Regression** `regressions/md/constant-d-polarization-feedback/`
   (out-of-library, public headers only, ≤ 60 s, `reg_` +
   `run_reg_` targets):
   - Gate A (static): ±1 e pair, d = 2 Å, L = 10 Å, D_z = 0.5 V/Å →
     P_z = 2e-3 e/Å², E_z = 1.3809744842959204e-1 V/Å,
     U_D = 5.2696244279044818e-2 eV (hard-coded literature-formula pins);
     foothold twin at the same D keeps E_z = 0.5 (screening response).
   - Gate B (dynamics): D_z = 0, ±1 e pair at rest, m = 1 amu, d0 = 2 Å
     → itinerant-polarization harmonic well U = 2π k_e q² d²/V, plasma-like
     oscillation ω = sqrt(4π k_e q²/(μV) · force_to_accel) =
     5.9091697154893236e-2 fs⁻¹ (T = 106.33 fs); pin d(t) at quarter /
     half period vs d0·cos(ωt), energy conservation |KE + U − U0| small.
   - Gate C (vacuum twin): same initial state, feedback off, D = 0 →
     zero force, atoms never move; E_z stays D_z. This is the queue's
     "E_z response vs vacuum-foothold twin" observable.
7. **Oracle**: row says oracle=metalwalls; no MetalWalls binary is
   installed and mainline LAMMPS has no constant-D ensemble → literature
   pin with full citation (Zhang & Sprik, Phys. Rev. B 93, 144201 (2016);
   Cox & Sprik, J. Chem. Phys. 151, 064506 (2019)), `perf_exempt = true`.
   Pin derivation saved under `tools/metalwalls_oracle/`
   (`inputs/literature/constant_d_pins.py` +
   `pin_notes_constant-d-polarization-feedback.md`); REFERENCE.md states
   oracle/oracle_version/citations/units.
8. **MD-only build** honest: HF/PW/GPW/CUDA/TORCH stay OFF; scoped
   targets only; dual gate (ctest unit filter + run_reg) green.

---
slug: gpw-d4-atm-forces
status: approved
created: 2026-07-31
domain: gpw / dispersion
queue: cp2k-method-loop (wave 1 2026-07-31)
depends_on: []
sibling: gpw-d3-atm-forces (parallel worker; shares only the ATM triple-derivative
  pattern conceptually — no shared code file is touched by both)
---

# Spec: `gpw-d4-atm-forces` — analytic D4 ATM (s9) ionic forces

**Queue row:** `gpw-d4-atm-forces` | Analytic ionic forces with DFT-D4
three-body ATM (s9) | gpw | GapCensus wave 14, re-verified wave 15.

## Gap

`gpw-d4-atm` shipped the ATM **energy** (dftd4 3.7.0 pin, neutral-charge C6
reweight). Forces with `d4_atm = true` work today, but the ATM leg of
`cpu::scf::dftd4()` is a **central finite difference** of the ATM energy
(`grad_h` parameter, `src/cpu/scf/dftd4.cpp`), re-running CN + neutral-C6 +
triple lattice sums 6N times. There is no oracle pin on the ATM gradient
anywhere in the tree — only the two-body `-g` pin (`gpw-d4-forces`).

## Deliverable

1. **Analytic ATM gradient** inside `cpu::scf::dftd4()` (the free function GPW
   forces / stress and molecular `D4::forces` all delegate to):
   - explicit triple-geometry derivatives of `E_ATM = Σ −rr·c9·triple` w.r.t.
     `r²ij, r²ik, r²jk` (Axilrod–Teller angle polynomial + zero-damping
     `fdmp`, alp = 16);
   - C9 chain: `∂E/∂c6 = ½ E_triple/c6`, folded through the **neutral-charge**
     weight derivative `∂gw0/∂CN_d4 = dgwk_dcn · ζ(q=0)` and the existing
     `dcn_d4` neighbour lists. No EEQ adjoint on the ATM leg — `q = 0` is
     geometry-independent by construction (dftd4 `disp.f90` neutral reweight).
   - Drop the now-dead `grad_h` parameter from `dftd4()` (no BC shims).
2. **Offline pins** (dftd4 3.7.0, produced 2026-07-31, stored under
   `tools/cp2k_oracle/inputs/gpw-d4-atm/`):
   - `pin_periodic_atm_grad.json` — periodic water dimer, `-f pbe --grad`
     (s9 = 1): E = −9.9554475834240310e−04 Ha + 6×3 gradient;
   - `pin_isolated_atm_grad.json` — same nuclei isolated, ATV atom order:
     E = −9.8819741533958229e−04 Ha + gradient.
3. **Unit tests** (`tests/cpu/gpw/test_cpu_gpw_dftd4.cpp`,
   `tests/cpu/scf/test_dispersion_gradient.cpp`): isolated + periodic
   s9 = 1 forces vs −gradient pins; analytic-vs-FD of the full s9 = 1
   functional tightened to the analytic level; ATM-increment
   (s9=1 − s9=0) vs oracle pin delta.
4. **Regression** `regressions/gpw/gpw-d4-atm-forces/` (+ REFERENCE.md):
   GPW PBE-D4+ATM product analytic forces on the periodic water dimer
   (12 Å, 200 Ry, DZVP-MOLOPT-SR-GTH), dual gate:
   - hard: standalone `dftd4(s9=1)` forces vs −gradient pin;
   - hard: `e_disp` vs the ATM energy pin (existing);
   - wiring: product(d4_atm) − product(plain) ≈ standalone D4+ATM forces.

## Oracle

dftd4 3.7.0 (`/software/sse2/el9_epyc9005/easybuild/prefix/software/dftd4/3.7.0-gfbf-2025b`),
the same libdftd4 CP2K links when built with DFTD4 (Arrhenius CP2K 2026.1
lacks it — same policy as `gpw-d4-atm` / `gpw-d4-forces`). Pins are numbers
the binary printed; re-run commands in REFERENCE.md.

## Tolerances

| Check | Tol |
| ----- | --- |
| standalone max\|F + g_pin\| (s9=1, periodic, direct) | 5e−8 Ha/Bohr |
| conservative parity \|F_x(O0) − FD(oracle E±)\| | 1e−9 Ha/Bohr |
| ATM increment max\|ΔF − Δpin\| | 5e−8 Ha/Bohr (non-vacuity max\|ΔF\| > 1e−8) |
| standalone max\|F + g_pin\| (s9=1, isolated, unit) | 1e−9 Ha/Bohr |
| \|E_disp − pin\| | 1e−8 Ha |
| wiring max\|(F_on − F_plain) − F_standalone\| | 1e−7 Ha/Bohr |
| unit analytic-vs-FD (s9=1 full functional) | 1e−7 Ha/Bohr |

**Oracle finding (2026-07-31, drove the gate split):** dftd4 3.7.0's
*periodic* ATM gradient is non-conservative by ~2.4e−8 Ha/Bohr against its
own printed energy (isolated is conservative to 2.6e−10; FD-verified
offline, commands in REFERENCE.md). ATV's s9=1 energy equals dftd4's at
base and displaced geometries to a constant 3.8e−12 Ha, so the direct pin
gate sits at 5e−8 (absorbs the defect, still 8× below the 4.05e−7 ATM
increment ⇒ fails if the ATM gradient is dropped), and the sharp 1e−9 gate
is against the FD of the oracle's own printed displaced energies
(`pin_periodic_atm_dispO0x_{plus,minus}.json`).

## Perf

`perf_exempt = true` — same policy and reason as `gpw-d4-forces` /
`gpw-d4-atm`: no CP2K QS ENERGY_FORCE twin on this host (CP2K built without
DFTD4); the hard deliverable is gradient-pin parity. The analytic path
replaces 6N ATM energy evaluations with one pass (strictly faster than
shipped FD).

## Out of scope

- D3 ATM forces (`gpw-d3-atm-forces`, parallel worker).
- ATM contribution to the **stress** tensor beyond what `gpw-d4-stress`
  already covers (energy-route FD stress unaffected).
- CUDA.

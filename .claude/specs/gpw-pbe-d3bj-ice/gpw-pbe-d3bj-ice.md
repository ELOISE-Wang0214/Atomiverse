---
slug: gpw-pbe-d3bj-ice
status: approved
created: 2026-07-23
track: cp2k-method-loop
chapter: gpw
grilled: false
---

# gpw-pbe-d3bj-ice

## Summary

Wire **Grimme DFT-D3 with Becke–Johnson (rational) damping** into the **GPW**
path for **periodic** water (ice-relevant H-bonded water dimer cell). Close the
gap left after zero-damping D3 on GPW (`d3`) and molecular BJ (`cpu::scf::D3BJ`).

Public surface: member `d3bj` on `cpu::gpw::GpwDft` and `cpu::gpw::GPWEvaluator`
(no `*Config` bag). Dispersion energy is exposed on `GpwResult::e_disp` [Ha]
and folded into `e_total`. Analytic D3(BJ) forces participate when
`gpw_forces_analytic` runs with `d3bj`.

## Domain basis

Grimme D3 + BJ damping (Grimme et al., *J. Comput. Chem.* **32**, 1456 (2011)):

\[
E_{\mathrm{disp}} = -\sum_{A<B,\,T}\left[
  \frac{s_6\,C_6^{AB}}{R_{AB,T}^6 + f^6}
  + \frac{s_8\,C_8^{AB}}{R_{AB,T}^8 + f^8}
\right],
\quad
f = a_1\sqrt{C_8/C_6}+a_2.
\]

PBE-D3(BJ) parameters (CP2K `REFERENCE_FUNCTIONAL PBE`):  
\(s_6=1\), \(s_8=0.7875\), \(a_1=0.4289\), \(a_2=4.4407\) Bohr.

Periodic orthorhombic Γ lattice sum (same image enumeration as zero-damping
GPW D3). Three-body ATM **off** (CP2K default with `CALCULATE_C9_TERM .FALSE.`).

**CN form for GPW / CP2K parity:** use CP2K's `cnparam_d3` (includes the
long-range \(\frac12(1-\tanh(R-2 r_{\mathrm{cov}}))\) factor) for the GPW
BJ path, so E_disp matches offline CP2K `TYPE DFTD3(BJ)`. The molecular
`dftd3_pbe_bj` / `cpu::scf::D3BJ` path keeps the original Grimme CN (PySCF /
simple-dftd3 twin) — do not break that pin.

## Design

| Symbol | Role |
| ------ | ---- |
| `GpwDft::d3bj` | bool member; enables BJ dispersion on SCF energy |
| `GpwDft::d3` | existing zero-damping flag; mutually exclusive with `d3bj` (throw if both) |
| `GPWEvaluator::d3bj` | threads into solver at `build` |
| `GpwResult::e_disp` | dispersion energy [Hartree] (0 when off) |
| `dftd3_pbe_bj_cp2k` (internal) or flag on BJ engine | CP2K-CN BJ for GPW |

Reuse in-tree `D3Params` / `dftd3_tables.h`. No public Config. No factories.

### Reuse decision

- `reuse` `dftd3_pbe` infrastructure (getc6, image_counts, force CN-chain).
- `reuse` molecular BJ formula; **change CN for GPW CP2K path only**.
- `pattern` zero-damping GPW `d3` wiring in `gpw_solver.cpp` / `gpw_forces.cpp`.
- `new` `d3bj` members + `e_disp` on result; unit + e2e regression.

## Geometry (pin)

Water dimer [Å] (same as zero-damping CP2K D3 unit test), cubic cell
**12 Å**, **PERIODIC XYZ** (ice-relevant H-bond pair in a periodic box; not a
full ice Ih supercell — keeps e2e ≤ 60 s):

| Atom | x | y | z |
| ---- | - | - | - |
| O | 0.000000 | 0.000000 | 0.000000 |
| H | 0.758602 | 0.000000 | 0.504284 |
| H | 0.758602 | 0.000000 | −0.504284 |
| O | 3.000000 | 0.000000 | 0.000000 |
| H | 3.758602 | 0.000000 | 0.504284 |
| H | 3.758602 | 0.000000 | −0.504284 |

Basis for optional full-energy twin: `DZVP-MOLOPT-SR-GTH` + `GTH-PBE-q*`.

## Offline CP2K pin

- **Engine:** CP2K **2026.1** (`git:5e54ba2`), module tree
  `CP2K/2026.1-gmpflf-2026.03` on Arrhenius.
- **Input:** `tools/cp2k_oracle/inputs/gpw-pbe-d3bj-ice/water_dimer_d3bj.inp`
- **Settings:** GPW PBE, `TYPE DFTD3(BJ)`, `REFERENCE_FUNCTIONAL PBE`,
  `R_CUTOFF 20.0` Bohr, no C9, no LRC, `CUTOFF 300` Ry.
- **E_disp:** **−0.00137627769481 Ha**
- **E_tot (CP2K QS):** −34.33962825563683 Ha (context only; absolute GPW
  e_tot is **not** the hard gate — SCF/OT/grid differences).
- Recorded in `regressions/gpw/gpw-pbe-d3bj-ice/REFERENCE.md`.

## Tasks

- [x] Spec + acceptance (this file + `.acceptance.md`)
- [x] CP2K offline pin → REFERENCE.md under `regressions/gpw/gpw-pbe-d3bj-ice/`
- [x] CP2K-CN BJ path for GPW (`dftd3_pbe_bj` with CP2K CN, or dedicated entry)
- [x] Wire `GpwDft::d3bj`, `GpwResult::e_disp`; mutually exclusive with `d3`
- [x] Wire forces (`gpw_forces_analytic`) and `GPWEvaluator::d3bj`
- [x] Unit tests: isolated/periodic E_disp vs pin + FD forces (`tests/cpu/gpw/`)
- [x] E2e regression `regressions/gpw/gpw-pbe-d3bj-ice/` (public API only)
- [x] CMake: `atv_add_regression(gpw gpw-pbe-d3bj-ice)` under `ATV_BUILD_GPW`
- [x] Build + unit + regression green; perf note vs CP2K wall (~13.4 s)

## Testing strategy

- **Unit:** water-dimer periodic E_disp vs CP2K pin (1e-8 Ha); analytic
  forces vs FD; `d3`+`d3bj` both true throws; `e_total` includes `e_disp`.
- **Regression:** public `GPWEvaluator` or direct `GpwDft` path reports
  `e_disp` within tol of pin; integration `Δe_total = e_disp`.

## Out of scope

- Full ice Ih supercell / phonons
- DFT-D4 / ATM three-body
- Molecular `D3BJ` / PySCF pin changes
- k-points, UKS, stress
- CUDA GPW

---

**Librarian note:** `cpu/gpw/`; confidence high (engine exists; wire + CP2K CN).

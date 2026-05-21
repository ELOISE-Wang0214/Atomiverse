---
title: Phase 4 Vibrations 03 — Harmonic thermochemistry + output components
status: approved
created: 2026-05-21
---

# Phase 4 Vibrations 03 — Harmonic thermochemistry + output components

## Summary

This is the final link of the three-spec vibrational-analysis chain. Given the normal-mode frequencies produced by link 2 and stored in `HessianResult`, it completes the analysis by computing harmonic-oscillator / rigid-rotor / translation thermochemistry — zero-point energy, internal energy U, enthalpy H, entropy S, and Gibbs free energy G — and by exposing two user-facing output components. After a Hessian run, a user obtains a parseable TSV summary of frequencies, imaginary-mode count, ZPE, and G, plus optional per-mode XYZ animations suitable for visualization. The thermochemistry is a faithful port of PySCF `pyscf/hessian/thermo.py::thermo()`, so results match PySCF to the published reference values within the acceptance tolerances.

## Domain basis

All thermodynamic quantities are reported in Hartree (Eh) and Eh/K, consistent with PySCF `thermo()`. Constants come from the NIST-mirrored module introduced in link 2 (CODATA-2014, matching PySCF): `PLANCK = 6.626070040e-34` J·s, `BOLTZMANN = 1.38064852e-23` J/K, `AVOGADRO = 6.022140857e23` /mol, `LIGHT_SPEED_SI = 299792458` m/s, plus `ATOMIC_MASS` and `HARTREE2J`. Per-element masses follow PySCF `elements.MASSES` (H 1.008, C 12.011, N 14.007, O 15.999).

Reference (binding, port verbatim): PySCF `pyscf/hessian/thermo.py` — `thermo()`, `harmonic_analysis()`, `rotation_const()`, `_get_rotor_type()`. PySCF: Sun et al., *J. Chem. Phys.* 153, 024109 (2020), DOI 10.1063/5.0006074. Statistical-thermodynamics formulas: McQuarrie, *Statistical Mechanics* (2000); Sackur-Tetrode equation.

Gas constant in atomic units (per molecule): `R_Eh = BOLTZMANN / HARTREE2J` [Eh/K].

Per-mode vibrational temperature and helpers, for real positive frequencies only (`idx = freq.real > 0`; imaginary modes dropped):
- `theta_i = h * nu_i / kB` where `nu_i` is the link-2 mode frequency converted to Hz (`freq_au_i * au2hz`); `theta_i` in K.
- `rt_i = theta_i / T`; `x_i = exp(-theta_i / T)`.

Zero-point energy:
- `ZPE = R_Eh * 0.5 * sum_i theta_i` [Eh].

Vibrational partition-function contributions:
- `S_vib = R_Eh * sum_i ( rt_i * x_i / (1 - x_i) - ln(1 - x_i) )`
- `E_vib = H_vib = ZPE + R_Eh * T * sum_i rt_i * x_i / (1 - x_i)`
- `Cv_vib = R_Eh * sum_i x_i * rt_i^2 / (1 - x_i)^2`

Translational contributions (M = total molecular mass in kg):
- `q_trans = (2*pi*M*kB*T / h^2)^1.5 * (kB*T / P)`
- `S_trans = R_Eh * (2.5 + ln q_trans)` (Sackur-Tetrode)
- `E_trans = 1.5 * R_Eh * T`; `H_trans = 2.5 * R_Eh * T` (the `+kB*T` enthalpy term lives in the translational component).

Rotational contributions — rotor type from `_get_rotor_type`; principal moments of inertia from diagonalizing the inertia tensor about the center of mass; rotational constants `B_A,B_B,B_C` from the inertia eigenvalues (`rotation_const()`):
- Nonlinear: `q_rot = (kB*T/h)^1.5 * sqrt(pi) / (sigma * sqrt(B_A B_B B_C))`; `S_rot = R_Eh * (1.5 + ln q_rot)`; `E_rot = H_rot = 1.5 * R_Eh * T`.
- Linear: `q_rot = kB*T / (sigma * h * B)`; `S_rot = R_Eh * (1 + ln q_rot)`; `E_rot = H_rot = R_Eh * T`.
- Single atom: all rotational terms 0.

Electronic contributions:
- `S_elec = R_Eh * ln(multiplicity)`; `E_elec = H_elec = E0` (SCF total energy). Closed-shell singlets (multiplicity 1) give `S_elec = 0`.

Totals:
- `S_tot = S_elec + S_trans + S_rot + S_vib`
- `E_tot = E_elec + ZPE + (E_vib - ZPE) + E_trans + E_rot`
- `H_tot = H_elec + H_vib + H_trans + H_rot`
- `G_tot = H_tot - T * S_tot`

The acceptance compares the **thermal correction to Gibbs**, `G_tot - E0`, which is independent of the electronic energy and isolates the thermochemistry machinery; the simplest sub-check is `E_0K - E0 = ZPE`.

Symmetry number sigma enters `q_rot` as a `1/sigma` factor. A wrong sigma shifts G by `kB*T*ln(sigma)` (~2.4 mEh for sigma=12 at 298 K). Point-group / sigma auto-detection is explicitly out of scope for this link: sigma is **supplied per molecule** as a `Thermochem` config input. Test-molecule symmetry numbers: H2O = 2, NH3 = 3, CH4 = 12, benzene = 12, formaldehyde = 2.

## Design

Three new components in the `ATV::cpu::vib` namespace, all consuming the `HessianResult` struct populated by link 2 (which already carries the real/imaginary normal-mode frequencies). None of the three introduce a new driver phase enum — they bind to the existing `vib::Phase` slots defined in links 1-2.

- `cpu::vib::Thermochem` — `cpu::ComponentBase`, `role() == ComponentRole::Observer`, bound to `vib::Phase::thermochem`. Holds a `Thermochem::Config { double temperature = 298.15; double pressure = 101325.0; }` and consumes per-molecule structural inputs (atomic masses, Cartesian coordinates, multiplicity, symmetry number sigma) from the frame / molecule. On `launch()` it performs the faithful PySCF `thermo()` port: it computes the electronic, translational, rotational, and vibrational components separately (mirroring PySCF's component-wise `dict` structure and `_sum` aggregation) and writes ZPE, U, H, S, G, and `n_imaginary` into the `HessianResult` struct. Rotor type is determined via a `_get_rotor_type` port from the inertia eigenvalues; `rotation_const()` is ported to derive `B` from those eigenvalues. No allocation occurs in `launch()` — scratch buffers for the inertia tensor and per-mode arrays are sized in `build()`.
- `cpu::vib::FrequencyReport` — `cpu::ComponentBase`, `role() == ComponentRole::Recorder`, bound to `vib::Phase::end_of_step`. Constructed with an output path (default `"frequencies.tsv"`). It gates on the final displacement pass (`ctx.step == 6N`) so it fires exactly once after `thermochem` has populated `HessianResult`. It writes a tab-separated file: one row per frequency (index, frequency value, real/imaginary flag) plus a trailer summarizing imaginary-mode count, ZPE, and G(T). Writing only — it computes nothing.
- `cpu::vib::NormalModeWriter` — `cpu::ComponentBase`, `role() == ComponentRole::Recorder`, bound to `vib::Phase::end_of_step`. Optional; constructed with an output-path prefix and an amplitude scale. It also gates on `ctx.step == 6N`. For each real mode it writes a multi-frame XYZ file animating the displacement of the equilibrium geometry along that normal-mode eigenvector (a sinusoidal sweep), reusing the existing `cpu::recorder` XYZ writer primitives for atom-symbol formatting. Writing only.

Lifecycle for all three: `build()` (one-time scratch alloc / open-file handles if applicable) → `launch()` (per-call compute or I/O) → `destroy()` (flush and close). Single responsibility per CLAUDE.md: `Thermochem` computes derived quantities (Observer); `FrequencyReport` and `NormalModeWriter` only perform I/O (Recorder).

## Files to create or modify

- `include/atomiverse/cpu/vib/thermochem.h` (new) — `Thermochem` component + `Thermochem::Config`.
- `src/cpu/vib/thermochem.cpp` (new) — PySCF `thermo()` port: elec / trans / rot / vib component functions, rotor-type and inertia helpers, total aggregation.
- `include/atomiverse/cpu/vib/frequency_report.h` (new) — `FrequencyReport` Recorder.
- `src/cpu/vib/frequency_report.cpp` (new) — TSV writer.
- `include/atomiverse/cpu/vib/normal_mode_writer.h` (new) — `NormalModeWriter` Recorder.
- `src/cpu/vib/normal_mode_writer.cpp` (new) — per-mode XYZ animation writer.
- `include/atomiverse/cpu/vib/hessian_result.h` — extend the link-2 `HessianResult` struct with thermochemistry fields (`zpe`, `u_tot`, `h_tot`, `s_tot`, `g_tot`, `g_thermal`, `n_imaginary`).
- `tests/cpu/vib/test_thermochem_vs_pyscf.cpp` (new) — ZPE / G-correction / per-component checks vs PySCF references.
- `tests/cpu/vib/test_vib_output_components.cpp` (new) — TSV and XYZ output round-trip checks.
- `tests/CMakeLists.txt` — register the two new test targets.
- `benchmarks/qc/reference/thermo/pyscf_thermo_references.json` (new) — checked-in PySCF `thermo()` reference values with `pyscf_version` provenance.
- `benchmarks/qc/reference/generate_references.py` — add a generator stanza emitting the PySCF `thermo()` references for H2O / NH3 / CH4 (and benzene / formaldehyde if tractable).

## Tasks

- [ ] Extend HessianResult with thermochemistry fields in include/atomiverse/cpu/vib/hessian_result.h
- [ ] Generate and check in PySCF thermo() reference data (benchmarks/qc/reference/thermo/pyscf_thermo_references.json, generate_references.py)
- [ ] Write failing tests for Thermochem vs PySCF (tests/cpu/vib/test_thermochem_vs_pyscf.cpp)
- [ ] Implement cpu::vib::Thermochem component porting PySCF thermo() (include/atomiverse/cpu/vib/thermochem.h, src/cpu/vib/thermochem.cpp)
- [ ] Write failing tests for FrequencyReport and NormalModeWriter output (tests/cpu/vib/test_vib_output_components.cpp)
- [ ] Implement cpu::vib::FrequencyReport TSV Recorder (include/atomiverse/cpu/vib/frequency_report.h, src/cpu/vib/frequency_report.cpp)
- [ ] Implement cpu::vib::NormalModeWriter XYZ-animation Recorder (include/atomiverse/cpu/vib/normal_mode_writer.h, src/cpu/vib/normal_mode_writer.cpp)
- [ ] Add Doxygen docstrings per doc.style with units to all three component headers and register test targets in tests/CMakeLists.txt
- [ ] Run full check + test suite

## Testing strategy

Happy path:
- `Thermochem` on H2O / NH3 / CH4 produces ZPE matching the checked-in PySCF reference to <= 0.001 eV (~0.1 kJ/mol).
- `Thermochem` produces a thermal Gibbs correction `G_tot - E0` at 298.15 K, 1 atm matching PySCF to <= 0.5 kJ/mol on the same molecule set (benzene / formaldehyde added if runtime allows).
- `FrequencyReport` writes a TSV that parses cleanly and contains the correct frequency count plus the ZPE value.
- `NormalModeWriter` writes a valid multi-frame XYZ that loads with the correct atom count and atom symbols matching the input molecule.

Edge cases:
- Imaginary modes are dropped from all vibrational sums (`idx = freq.real > 0`); `n_imaginary` is reported separately and surfaced in the TSV trailer.
- Linear molecule rotor branch (2 rotational DOF, `S_rot = R_Eh*(1 + ln q_rot)`) versus nonlinear branch is selected correctly by the `_get_rotor_type` port; single-atom case yields zero rotational terms.
- Closed-shell singlet (multiplicity 1) yields `S_elec = 0`.
- Supplied sigma is honored verbatim — a deliberately wrong sigma in a fixture shifts G by `kB*T*ln(sigma)`, confirming sigma is an input not auto-detected.

Domain validation:
- Component-resolved check on at least H2O: `S_trans` (Sackur-Tetrode), `S_rot`, and `S_vib` each match PySCF's per-component `thermo()` values to tight tolerance, so a failing total can be traced to the specific partition function at fault.
- Sub-check `E_0K - E0 == ZPE` holds (electronic-energy-independent identity).
- Component roles verified: `Thermochem.role() == Observer`; `FrequencyReport.role() == NormalModeWriter.role() == Recorder`.

## Out of scope

- IR / Raman intensities (require dipole derivatives).
- Anharmonic corrections (VPT2, vibrational SCF).
- Quasi-RRHO low-frequency-mode entropy treatment (Grimme damping).
- Point-group detection and symmetry-number auto-detection — sigma is supplied per molecule as a `Thermochem` config input.
- CUDA backend port of the vibrational thermochemistry pipeline.
- Analytic Hessian and the numerical Hessian / normal-mode machinery themselves (delivered by links 1 and 2 of this chain).

---
title: Phase 4 Vibrations — Normal-mode frequencies (mass-weighting, TR projection, diagonalization)
status: approved
created: 2026-05-21
---

# Phase 4 Vibrations — Normal-mode frequencies

## Summary

This spec turns the raw numerical Hessian produced by link 1 (`phase4-vibrations-01-numerical-hessian`) into physical vibrational frequencies. It adds two new components to the `cpu::vib` domain — `MassWeighting` and `cpu::vib::EigenSolver` — both bound to `vib::Phase::diagonalize`. Together they mass-weight the Hessian, project out the 3 translational and 2–3 rotational modes, diagonalize the projected mass-weighted Hessian, and convert the eigenvalues to wavenumbers (cm⁻¹). The resulting frequencies and normal-mode eigenvectors are written back into the `HessianResult` struct that link 1 defined. After this spec lands, a user running the `cpu::vib::Hessian` driver at an optimized geometry obtains exactly 3N−6 (nonlinear) or 3N−5 (linear) genuine vibrational frequencies, with imaginary modes reported as negative wavenumbers, matching PySCF `harmonic_analysis` to within 5 cm⁻¹ for modes above ~300 cm⁻¹.

## Domain basis

All formulas below are ported verbatim from PySCF `pyscf/hessian/thermo.py::harmonic_analysis` and the constants from `pyscf/data/nist.py` (CODATA-2014). The reference-port-fidelity rule applies — port faithfully, do not paraphrase the numerics.

**Mass-weighted Hessian.** Given the raw numerical Hessian `H` from link 1 in Hartree/Bohr²:

```
H^mw_{iα,jβ} = H_{iα,jβ} / sqrt(m_i * m_j)        [Hartree / Bohr² / amu]
```

with `H` in Hartree/Bohr² and atomic masses `m_i` in amu. (Link 1 records the Hessian in eV/Å²; the unit conversion to Hartree/Bohr² is applied here, using the constants module below.)

**Eigenproblem.** `H^mw` is real-symmetric, giving real eigenvalues `λ` (force constants in atomic units) and orthonormal eigenvectors `v`:

```
H^mw v = λ v
```

**Frequency conversion.** Frequency `ω = sqrt(λ)`. Wavenumber conversion reproduces PySCF's exact two-factor form — no hard-coded `5140.x` literal:

```
au2hz           = (HARTREE2J / (ATOMIC_MASS * BOHR_SI^2))^0.5 / (2*pi)
freq_wavenumber = freq_au * au2hz / LIGHT_SPEED_SI * 1e-2
```

with CODATA-2014 constants (mirrored verbatim from `pyscf/data/nist.py`):

```
HBAR           = 6.626070040e-34 / (2*pi)
E_MASS         = 9.10938356e-31
BOHR_SI        = 0.52917721092e-10
HARTREE2J      = HBAR^2 / (E_MASS * BOHR_SI^2)
ATOMIC_MASS    = 1e-3 / 6.022140857e23
LIGHT_SPEED_SI = 299792458
```

The combined constant is ≈ 2.194746e5 cm⁻¹ per unit `freq_au`. Rounding drift in any constant shifts frequencies by ~0.1 cm⁻¹, so a `nist`-mirrored constants module is mandatory.

**Imaginary frequencies.** When `λ < 0`, report the wavenumber as **negative** (`-sqrt(|λ|) * au2hz * ...`), the `imaginary_freq=False` convention universal in chemistry software. The count of imaginary (negative) modes is a reported diagnostic.

**Translation/rotation projection.** At a stationary point the 3N modes split into 3N−6 genuine vibrations (3N−5 for linear molecules) plus 3 translations and 3 rotations (2 for linear). PySCF's `harmonic_analysis` projects these out before diagonalizing:

1. Build the 6 (or 5) mass-weighted translation/rotation vectors (`_get_TR`): translations are `sqrt(m)*[unit axis]`; rotations are derived from the atomic coordinates expressed in the principal-axis frame.
2. QR-orthonormalize the TR vectors → `Q`.
3. Form the projector `P = I − Q Qᵀ`.
4. Diagonalize `H^mw` restricted to the complement subspace (`bvec`), yielding exactly 3N−6 / 3N−5 finite frequencies.

Linear-vs-nonlinear is decided by `_get_rotor_type`: `rot_const[0] > 1e8` and `rot_const[1] − rot_const[2] < 1e-3` ⇒ LINEAR.

**Atomic masses.** Masses are PySCF `elements.MASSES` (IUPAC-2013 standard atomic weights, isotope-averaged): H 1.008, C 12.011, N 14.007, O 15.999. A mass mismatch shifts every frequency (`ν ∝ m^(−1/2)`); using 12.000 for C costs ~0.7 cm⁻¹ at 1500 cm⁻¹.

**Expected accuracy.** Numerical Hessian (central difference, `h = 1e-3 Å`, tight SCF) vs PySCF analytic-Hessian frequencies: ≤ 5 cm⁻¹ for modes above ~300 cm⁻¹. Low modes are worst-conditioned (`δν ∝ 1/sqrt(λ)`); the 5 cm⁻¹ cap is applied over the projected 3N−6 / 3N−5 modes.

References:
- PySCF `pyscf/hessian/thermo.py` — `harmonic_analysis`, `_get_TR`, `_get_rotor_type`. PySCF: Sun et al., *WIREs Comput. Mol. Sci.* 8 (2018) e1340, DOI 10.1002/wcms.1340.
- CODATA-2014 fundamental constants — Mohr, Newell, Taylor, *Rev. Mod. Phys.* 88 (2016) 035009, DOI 10.1103/RevModPhys.88.035009.

## Design

**Constants module.** A new `cpu/vib/nist_constants.h` mirrors `pyscf/data/nist.py` verbatim: the CODATA-2014 raw constants listed above and the derived `HARTREE2J`, `au2hz`, and the `freq_au → cm⁻¹` factor, all as `constexpr double`. A new `cpu/vib/atomic_masses.h` holds the `elements.MASSES` IUPAC-2013 table keyed by atomic number, exposing a `mass_amu(int Z)` lookup that throws on an unsupported element (no fallback per project convention).

**`MassWeighting` component** (`cpu::vib::MassWeighting`, `role()==None`, `capabilities()` empty — it is a structural transform, not a frame-field writer). Bound to `vib::Phase::diagonalize` immediately after link 1's `Symmetrize` (order −1) at the default order. It reads atomic numbers from the frame to fetch per-atom masses, then divides the symmetrized Hessian `H` (held in `HessianResult` / the driver workspace from link 1) by `sqrt(m_i m_j)` element-wise, producing `H^mw` into the same workspace. No allocation in `launch()` — the workspace is sized in `build()`.

**`EigenSolver` component** (`cpu::vib::EigenSolver`, namespace-scoped — distinct from `cpu::scf::internal::EigenSolver`, no clash; `role()==None`, `capabilities()` empty). Bound to `vib::Phase::diagonalize` at order +1 (last). It: (a) determines linear vs nonlinear via `_get_rotor_type` on the principal-axis rotational constants; (b) builds the 5/6 mass-weighted TR vectors via `_get_TR`, QR-orthonormalizes them, forms `P = I − QQᵀ`; (c) restricts `H^mw` to the TR-complement subspace and diagonalizes it with the real-symmetric eigensolver (xtensor-blas `eigh`, LAPACK `dsyevd`, already linked via the SCF path); (d) converts eigenvalues to wavenumbers using `nist_constants.h` (negative for `λ < 0`); (e) writes `frequencies` (length 3N−6 / 3N−5, cm⁻¹), `normal_modes` (eigenvectors lifted back to the full 3N Cartesian basis), `n_imaginary`, and `is_linear` into the `HessianResult` struct link 1 defined.

**`HessianResult` extension.** Link 1 declared `HessianResult` with the raw `H`. This spec adds the fields `frequencies`, `normal_modes`, `n_imaginary`, and `is_linear`. The struct lives in link 1's `cpu/vib/hessian_result.h`; this spec edits that header to add the new members (link 1 is a hard dependency and assumed landed).

**Component ordering on `vib::Phase::diagonalize`:** `Symmetrize` (link 1, order −1) → `MassWeighting` (this spec, default order 0) → `EigenSolver` (this spec, order +1). The `cpu::vib::Hessian` driver's `default_bindings()` (from link 1) is extended to register `MassWeighting` and `EigenSolver`.

**Reference data decision.** The parity test compares Atomiverse numerical-Hessian frequencies against **PySCF's analytic-Hessian** `harmonic_analysis` frequencies, treated as the physical ground truth. The 5 cm⁻¹ tolerance absorbs the numeric-vs-analytic gap. Reference frequencies are generated by a genuine PySCF run and checked in as a JSON file (`tests/cpu/vib/reference/pyscf_frequencies.json`) carrying `pyscf_version`, method/basis, geometry, and per-molecule frequency arrays for provenance.

## Files to create or modify

- `include/atomiverse/cpu/vib/nist_constants.h` (new) — CODATA-2014 constants mirrored from `pyscf/data/nist.py`; derived `au2hz` and cm⁻¹ factor.
- `include/atomiverse/cpu/vib/atomic_masses.h` (new) — IUPAC-2013 `elements.MASSES` table + `mass_amu(int Z)` throwing lookup.
- `include/atomiverse/cpu/vib/mass_weighting.h` (new) — `cpu::vib::MassWeighting` component.
- `include/atomiverse/cpu/vib/eigen_solver.h` (new) — `cpu::vib::EigenSolver` component (TR projection + diagonalization + cm⁻¹ conversion).
- `src/cpu/vib/mass_weighting.cpp` (new) — `MassWeighting` implementation.
- `src/cpu/vib/eigen_solver.cpp` (new) — `EigenSolver` implementation: `_get_TR`, `_get_rotor_type`, QR projector, `dsyevd`, conversion.
- `include/atomiverse/cpu/vib/hessian_result.h` — add `frequencies`, `normal_modes`, `n_imaginary`, `is_linear` members (file owned by link 1).
- `src/cpu/vib/hessian.cpp` — extend `default_bindings()` to register `MassWeighting` and `EigenSolver` on `vib::Phase::diagonalize` (file owned by link 1).
- `tests/cpu/vib/test_mass_weighting.cpp` (new) — unit test for the mass-weighting transform.
- `tests/cpu/vib/test_tr_projection.cpp` (new) — unit test for `_get_TR` / `_get_rotor_type` / projector and genuine mode count.
- `tests/cpu/vib/test_frequencies_vs_pyscf.cpp` (new) — frequency parity test against checked-in PySCF references.
- `tests/cpu/vib/reference/pyscf_frequencies.json` (new) — PySCF `harmonic_analysis` reference frequencies with `pyscf_version` provenance.
- `tests/cpu/vib/reference/generate_pyscf_frequencies.py` (new) — PySCF script that regenerates the reference JSON.
- `tests/CMakeLists.txt` — register the three new `test_cpu_vib_*` targets.

## Tasks

- [ ] Write failing tests for nist constants and atomic masses (tests/cpu/vib/test_mass_weighting.cpp)
- [ ] Implement nist_constants.h and atomic_masses.h (include/atomiverse/cpu/vib/nist_constants.h, include/atomiverse/cpu/vib/atomic_masses.h)
- [ ] Write failing tests for MassWeighting transform (tests/cpu/vib/test_mass_weighting.cpp)
- [ ] Implement MassWeighting component and extend HessianResult (include/atomiverse/cpu/vib/mass_weighting.h, src/cpu/vib/mass_weighting.cpp, include/atomiverse/cpu/vib/hessian_result.h)
- [ ] Write failing tests for TR projection, rotor-type, and genuine mode count (tests/cpu/vib/test_tr_projection.cpp)
- [ ] Implement EigenSolver component with _get_TR, _get_rotor_type, QR projector, dsyevd, cm-1 conversion (include/atomiverse/cpu/vib/eigen_solver.h, src/cpu/vib/eigen_solver.cpp)
- [ ] Register MassWeighting and EigenSolver in Hessian default_bindings on vib::Phase::diagonalize (src/cpu/vib/hessian.cpp)
- [ ] Generate PySCF reference frequencies and write failing parity test (tests/cpu/vib/reference/generate_pyscf_frequencies.py, tests/cpu/vib/reference/pyscf_frequencies.json, tests/cpu/vib/test_frequencies_vs_pyscf.cpp)
- [ ] Register the three new vib test targets in tests/CMakeLists.txt and add doxygen with units to all new public symbols
- [ ] Run full check + test suite

## Testing strategy

- **Happy path — mass-weighting.** `test_mass_weighting.cpp`: a hand-built 6×6 Hessian on a 2-atom system with known masses; assert each `H^mw_{iα,jβ} == H_{iα,jβ}/sqrt(m_i m_j)` to machine precision. Assert `mass_amu(Z)` returns 1.008/12.011/14.007/15.999 for H/C/N/O and throws for an unsupported `Z`.
- **Happy path — TR projection.** `test_tr_projection.cpp`: on an optimized H₂O geometry, build the TR vectors, QR-orthonormalize, project; assert the projector annihilates the 6 TR vectors and that diagonalizing the projected `H^mw` yields exactly 3N−6 = 3 finite modes. Repeat for a linear molecule (CO₂ or HCN) and assert `is_linear == true` and 3N−5 modes.
- **Edge case — imaginary modes.** `test_frequencies_vs_pyscf.cpp` or `test_tr_projection.cpp`: feed a deliberately non-stationary geometry (or an NH₃ inversion saddle if a geometry is available) and assert exactly one mode reports a negative wavenumber and `n_imaginary == 1`.
- **Edge case — constants integrity.** `test_mass_weighting.cpp`: assert the derived `au2hz`-based cm⁻¹ factor equals ≈ 2.194746e5 within tight tolerance, confirming the `nist`-mirrored constants compose correctly.
- **Domain validation.** `test_frequencies_vs_pyscf.cpp`: run the full `cpu::vib::Hessian` driver on H₂O, NH₃, CH₄ at a tractable level (RHF or PBE / 6-31G); assert max frequency deviation over the 3N−6 / 3N−5 projected modes is ≤ 5 cm⁻¹ vs the checked-in PySCF analytic-Hessian references. Include benzene or formaldehyde if wall-time permits, otherwise note runtime.
- **Regression.** Full `ctest` suite stays green after the new components and `default_bindings()` extension.

## Out of scope

- Thermochemistry (ZPE, U/H/S/G, partition functions) — deferred to link 3 of the chain.
- IR / Raman intensities — requires dipole derivatives, not addressed here.
- Analytic Hessian (CPHF + second ERI derivatives) — link 1 produces a numerical Hessian; analytic is a separate future spec.
- CUDA backend — this domain is CPU-only for the chain.
- The numerical Hessian itself, the `cpu::vib::Hessian` driver, the `vib::Phase` enum, and the base `HessianResult` struct — all delivered by link 1 and assumed landed; this spec only consumes and extends them.

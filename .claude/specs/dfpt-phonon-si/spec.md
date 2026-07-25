---
slug: dfpt-phonon-si
status: approved
created: 2026-07-25
revised: 2026-07-25
domain: response / plane-wave DFT
queue: cp2k-method-loop
deferred: false
depends_on: [fd-phonon-supercell]
unblocks: [phonon-dispersion-si]
grilled: false
---

# Spec: `dfpt-phonon-si` — Γ-point DFPT phonons for the PW method (Si, NC-UPF)

**Slug:** `dfpt-phonon-si`
**Domain:** response (plane-wave DFT linear response)
**Queue row:** `docs/developer/cp2k-method-queue.md` `[~] dfpt-phonon-si` (wave 3)
**Dep:** `fd-phonon-supercell` done (FD Γ phonon machinery, `cpu::vib::Phonon`)
**Unblocks:** `phonon-dispersion-si` (queue: "prefer after `dfpt-phonon-si`")

## Summary

Deliver **density-functional perturbation theory (DFPT) phonons at Γ** for the
concrete plane-wave method `cpu::pw::PW` (NC-UPF / Kleinman–Bylander, LDA):
analytic dynamical matrix from self-consistent first-order response
(Sternheimer equation), instead of `6N` finite-displacement SCF solves.
Target system: diamond Si, 2-atom fcc primitive cell, `tests/data/pw/Si_lda_nc.upf`
(the same NC pseudopotential every PW Si regression uses).

Roadmap: `docs/roadmap/dft.md` §11 DFPT ○ (Γ FD phonons ✓ via
`fd-phonon-supercell`). CP2K twin: `VIBRATIONAL_ANALYSIS` / DFPT; literature
anchor Baroni–de Gironcoli–Dal Corso–Giannozzi, *Rev. Mod. Phys.* **73**, 515
(2001) and Giannozzi et al., *Phys. Rev. B* **43**, 7231 (1991) (Si Γ optical
LDA ≈ 517 cm⁻¹, exp 517.3 cm⁻¹).

## Method (q = 0, restricted insulator, fixed occupations)

Force constants `C[Iα][Jβ] = ∂²E/∂R_Iα∂R_Jβ` split as:

1. **Ion–ion:** analytic second derivative of the Ewald energy
   (real + reciprocal parts; self term is R-independent).
2. **Electronic expectation (same atom only):**
   `Σ_k w_k Σ_i f ⟨ψ_i|∂²V_ext^{(I)}/∂R_α∂R_β|ψ_i⟩` with
   `V_ext = V_loc + V_NL` (local term evaluated in G-space from ρ(G) and the
   per-atom radial form factor v(G); KB term from projector derivative sums).
3. **Response:** `Σ_k w_k Σ_i f · 2 Re⟨Δψ_i^{(Jβ)}|∂V_ext^{(Iα)}|ψ_i⟩`, where
   `Δψ` solves the **self-consistent Sternheimer equation** in the
   conduction subspace:
   `P_c (H_k − ε_i) P_c |Δψ_i⟩ = −P_c ΔV_KS |ψ_i⟩`,
   `ΔV_KS = ΔV_ext + ΔV_H[Δn] + f_xc Δn`, with `ΔV_H(G) = 4π Δρ(G)/G²`
   (G = 0 term zero) and the **LDA kernel** `f_xc = ∂V_xc/∂n` evaluated
   pointwise from the same `evaluate_xc(XCType::LSDA, …)` the SCF uses
   (numerical pointwise derivative; documented, LSDA only).
   Linear mixing on Δn; projected preconditioned CG per (k, band).

Mass weighting and ω [cm⁻¹] conversion reuse the `cpu::vib` constants
(`nist::EV_PER_ANG2_TO_HARTREE_PER_BOHR2`, `nist::au2wavenumber()`,
`vib::mass_amu`) so FD and DFPT spectra are unit-identical by construction.

## Scope guards

- **LSDA only** (`method.xc == XCType::LSDA`); GGA response kernel is out of
  scope — throw `std::invalid_argument` otherwise.
- **UPF only** (`upf_path` non-empty); toy V_loc path throws.
- **Fixed occupations, restricted** (no smearing, no spin) — metallic DFPT
  (Fermi-level shifts) out of scope; throw.
- **q = 0 only.** Full dispersions are the `phonon-dispersion-si` follow-on.
- No ASR imposition: acoustic modes must come out near zero on their own
  (that *is* part of the validation).

## Public surface (architecture-conform)

`include/atomiverse/cpu/planewave/dfpt_phonon.h` — mirrors the `Bands`
precedent (`bands_dos.h`): a concrete class, parameters as members, no
`*Config` bag, no factory:

```cpp
struct DfptPhononResult {
    xt::xtensor<double, 1> frequencies;       // cm^-1, ascending, imag < 0
    xt::xtensor<double, 2> force_constants;   // eV/Å², [3N][3N], Cartesian
    bool scf_converged;                       // ground-state SCF
    bool response_converged;                  // all 3N response SCFs
    int n_atoms;
};

class DfptPhonon {
  public:
    double response_tol = 1e-8;   // Δn L2 convergence per perturbation
    int response_max_iter = 80;
    double cg_tol = 1e-10;        // Sternheimer CG residual
    int cg_max_iter = 200;
    double mix_beta = 0.5;        // response linear mixing

    DfptPhononResult compute(const PW& method, const cpu::Frame& frame) const;
};
```

Internals in `cpu/planewave/internal/dfpt.{h,cpp}` (private workspace; never
on Frame). Shared pieces promoted on their **second caller** (per the
extract-on-second-use rule):

- the minimal SCF driver already private to `pw_forces.cpp` (`run_scf`)
  becomes a declared internal function, and — found defect, fixed here —
  it must **honor `method.xc`** (today it silently evaluates the default PBE
  functional regardless of the requested `xc`; forces/stress users with the
  default were unaffected, but DFPT/LSDA would be silently inconsistent);
- `upf_vloc_to_gspace` (radial v(G) form factor, private to `pw_setup.cpp`)
  gets an internal declaration (DFPT needs the per-atom form factor —
  dividing `V_loc_G` by the structure factor is ill-defined on diamond's
  S(G) = 0 shells).

## Validation strategy

Two-sided, per the `fd-phonon-supercell` precedent:

1. **Internal consistency (sharp):** DFPT diagonal force constant vs the
   central **second difference of the SCF total energy** at identical
   settings (the true discrete Hessian), pinned offline in REFERENCE.md.
   *Revised during implementation:* FD of `pw_forces_analytic` was the
   original plan, but measurement showed the v1 analytic forces deviate
   from the exact discrete-energy gradient by ≈ 2 % (egg-box/aliasing;
   documented in REFERENCE.md) — FD-of-forces is not a valid oracle.
   DFPT matches the energy second difference to 0.07 %.
2. **Physics window (literature):** Si Γ optical triplet inside a cited
   window around the DFPT-LDA/experimental 517 cm⁻¹; acoustic near zero
   with no ASR imposition; optical triple degeneracy; force-constant
   row sums ≈ 0 (translational invariance).

Offline pin recorded in `regressions/response/dfpt-phonon-si/REFERENCE.md`
per `tools/cp2k_oracle/README.md` (no live oracle at build/regression time;
CP2K binary availability on this host recorded there).

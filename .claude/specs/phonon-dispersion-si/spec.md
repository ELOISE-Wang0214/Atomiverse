---
slug: phonon-dispersion-si
status: approved
created: 2026-07-25
revised: 2026-07-25
domain: response / phonons
queue: cp2k-method-loop
deferred: false
depends_on: [fd-phonon-supercell, dfpt-phonon-si]
unblocks: []
grilled: false
---

# Spec: `phonon-dispersion-si` — phonon dispersion along a q-path (Si)

**Slug:** `phonon-dispersion-si`  
**Domain:** response (periodic phonons)  
**Queue row:** `docs/developer/cp2k-method-queue.md` `[~] phonon-dispersion-si`  
**Deps:** `fd-phonon-supercell` ✓ (supercell FD IFCs), `dfpt-phonon-si` ✓ (Γ DFPT)  
**Roadmap:** `docs/roadmap/dft.md` §11 — full phonon dispersions (q ≠ 0) ○

## Summary

Deliver **phonon dispersion ω(q)** along a high-symmetry q-path for a crystal
unit cell, by **Fourier interpolation of supercell interatomic force constants
(IFCs)** obtained from the existing FD path (`cpu::vib::Phonon`). Target
system for the regression twin: diamond Si + Stillinger–Weber (same classical
force twin as `fd-phonon-supercell`), so the case stays ≤ 60 s and needs no
PW/DFPT SCF wall.

This closes the product gap left by Γ-only FD (`Phonon`) and Γ-only DFPT
(`cpu::pw::DfptPhonon`): a public ω(q) path, not merely zone-folded supercell Γ
modes.

Oracle posture (offline only, per `tools/cp2k_oracle/README.md`):

1. **Sharp, same-Hamiltonian:** ω(Γ) from the Fourier path must match the
   acoustic/optical physics of the supercell IFCs (acoustic |ω| ≈ 0 at Γ;
   optical band in the SW Si window).
2. **Physics:** acoustic branches rise from Γ; high-symmetry X and L points
   land inside literature/SW windows; no live CP2K at regression time.
3. **CP2K twin (when binary available):** `VIBRATIONAL_ANALYSIS` + force
   constants / band structure on the same SW or GTH setup — record wall and
   optional ratio in REFERENCE.md. Host this shift: no CP2K module →
   `perf_exempt`.

## Method

Given Cartesian supercell force constants Φ [eV/Å²] from
`cpu::vib::Phonon::run` on an \(n_a×n_b×n_c\) expansion of the **unit cell**
(atom order = `expand_supercell`), form the mass-weighted dynamical matrix of
the primitive cell at fractional wavevector **q**:

\[
D_{\alpha\beta}(\kappa\kappa';\mathbf{q})
  = \frac{1}{\sqrt{m_\kappa m_{\kappa'}}}
    \sum_{\ell'}
    \Phi_{\alpha\beta}(0\kappa,\,\ell'\kappa')
    \,e^{i\mathbf{q}\cdot\mathbf{R}_{\ell'}} ,
\]

convert Φ to Hartree/Bohr² with the same NIST factors as
`MassWeighting` / `PhononDiagonalizer`, diagonalize the Hermitian \(D(\mathbf{q})\),
and report \(\omega = \mathrm{sign}(\lambda)\sqrt{|\lambda|}\) in **cm⁻¹**.

Default path (FCC/diamond, fractional of the primitive reciprocal lattice):
**Γ–X–W–L–Γ–K**, same waypoints as `cpu::pw::Bands::set_diamond_path`.

## Public surface (architecture-conform)

`include/atomiverse/cpu/vibration/phonon_dispersion.h` — concrete class,
parameters as members, **no** `*Config` bag, **no** factory:

```cpp
struct PhononDispersionResult {
    std::vector<std::array<double, 3>> qfrac;
    std::vector<std::vector<double>> frequencies_cm; // [iq][branch], ascending
    int n_atoms_unit = 0;
    int n_branches = 0;
    int na = 1, nb = 1, nc = 1;
};

class PhononDispersion {
  public:
    std::vector<std::array<double, 3>> qfrac; // empty → set_diamond_path(4)

    void set_diamond_path(int n_seg = 4);

    /// @param hessian_sc Cartesian Φ [eV/Å²], shape [3 N_sc, 3 N_sc]
    /// @param unit       Unit-cell Frame (types → masses; simbox = lattice)
    /// @param na,nb,nc   Supercell factors used to build `hessian_sc`
    PhononDispersionResult
    compute(const xt::xtensor<double, 2>& hessian_sc, const cpu::Frame& unit,
            int na, int nb, int nc) const;
};
```

Composition (caller's job — no all-in-one façade):

```cpp
cpu::vib::Phonon phonon(/*na=*/2, /*nb=*/2, /*nc=*/2);
phonon.set_neighbor(...);
phonon.on(vib::Phase::force).add(cpu::StillingerWeber{pot});
phonon.build();
auto pr = phonon.run(unit);

cpu::vib::PhononDispersion path;
path.set_diamond_path(3);
auto disp = path.compute(pr.hessian, unit, pr.na, pr.nb, pr.nc);
// disp.frequencies_cm[iq][ibranch]  [cm^-1]
```

Internals: Hermitian eigensolve via the real 2n embedding of \(D\) (reuse
`xt::linalg::eigh` — no new complex LAPACK surface). No Frame pollution;
spectra live on the result type only.

## Scope guards

- **In:** Fourier IFC path from any Cartesian supercell Hessian that shares
  `expand_supercell` ordering (classical SW twin for the regression; PW/DFPT
  supercell IFCs are a composition follow-on, not this slug's wall budget).
- **Out:** finite-q DFPT Sternheimer (true DFPT@q≠0); GPW path this wave;
  ASR projection as a required step; public Config bags / factories; CUDA.

## Validation

| Gate | What |
| ---- | ---- |
| Unit | D(q) size / Hermiticity; Γ acoustic ~0 on a tiny analytic IFC; diamond path length |
| Regression | SW Si 2×2×2 primitive supercell; Γ acoustic ≤ 30 cm⁻¹; Γ optical in SW window; X and L points finite and ordered; hard-coded pins in REFERENCE.md |
| Perf | `perf_exempt` if no CP2K wall twin; else \(T_\mathrm{ATV}/T_\mathrm{CP2K} ≤ 10\) |

## Acceptance

See `acceptance.md`. Green dual gate ⇒ `acceptance_met=true`.

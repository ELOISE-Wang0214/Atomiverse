---
slug: pw-smearing-al
status: approved
created: 2026-07-23
track: cp2k-method-loop
domain: pw
grilled: false
---

# pw-smearing-al

## Summary

Ship **metallic electronic smearing** on the concrete plane-wave method
`cpu::pw::PW` so odd-electron / metallic systems (textbook **FCC Al**,
NC-UPF) can run a stable SCF with fractional occupations. Expose Fermi–Dirac
and Methfessel–Paxton (plus Gaussian) as **members** of `PW` — no public
`*Config` bag, no `make_*` factory. Pin occupations and total energy against
an **offline** CP2K twin (algorithm + documented settings) and unit-test the
occupation formulas against closed-form references.

## Domain basis

- Fermi–Dirac occupations: Mermin, *Phys. Rev.* **137**, A1441 (1965).
- Methfessel–Paxton: *Phys. Rev. B* **40**, 3616 (1989); QE `wgauss` /
  CP2K `&SMEAR` METHFESSEL-PAXTON.
- Spin-unpolarized closed-shell convention: band occupation \(f_b \in [0,2]\),
  \(n_e = \sum_{k,b} w_k f_b(\varepsilon_{kb},\mu,\sigma)\).
- Fermi level \(\mu\) chosen so \(n_e\) matches the valence electron count.
- Mermin free energy \(F = E_{\mathrm{KS}} - TS\) (entropy term from the
  chosen smearing) is the variational total energy when smearing is on.
- Units: public energy surface remains eV via Frame meta; `PWResult.energy`
  stays **Ry** (QE-comparable); `degauss` is **Ry** (QE/CP2K electronic
  temperature width convention).

## Design

### Public API (`cpu::pw::PW` members)

```cpp
enum class Smearing : int {
    none = 0,              // fixed integer occupations (default, legacy)
    fermi_dirac = 1,       // Fermi–Dirac
    methfessel_paxton = 2, // MP order mp_order (default 1)
    gaussian = 3,          // Gaussian (MP order 0)
};

Smearing smearing = Smearing::none;
double degauss = 0.02;   // smearing width [Ry]; required > 0 when smearing != none
int mp_order = 1;        // Methfessel–Paxton order N (≥1 when kind is MP)
int extra_bands = 4;     // bands beyond ceil(n_elec/2) when smearing is on
```

`PWResult` gains:

```cpp
double e_fermi = 0;      // Fermi level [Ry]
double n_elec_occ = 0;   // Σ w f after last SCF iteration
double e_entropy = 0;    // −TS contribution included in energy [Ry]
```

### Internal

- `internal::smearing.h` / `.cpp`: spinless \(f(x)\in[0,1]\), closed-shell
  \(2f\), Fermi bisection, entropy weight.
- `DensityBuilder::accumulate_kpoint` takes per-band occupations
  (no hard-coded `*2`).
- SCF loop: Davidson → collect eigenvalues (MPI Allgather when active) →
  find \(\mu\) → occupations → density + band energy with \(f_b\) → free
  energy with \(−TS\).
- Default `Smearing::none` **must** preserve prior Si/H regressions
  (integer occupations, `n_bands = max(n_elec/2,1)`).

### System for e2e

- Geometry: FCC Al primitive cell, \(a = 4.05\) Å, 1 atom at origin.
- PP: `tests/data/pw/Al_pz_nc.upf` (NC, LDA/PZ, \(Z_{\mathrm{val}}=3\)).
- Settings (regression): `ecut_ry=20`, `k_mesh=2`, `smearing=fermi_dirac`,
  `degauss=0.02` Ry, `extra_bands=4`, `tol=1e-5`, `max_iter=60`.

### Offline CP2K twin

- CP2K 2026.1 Quickstep Γ/k-mesh Al with `&SMEAR METHOD FERMI_DIRAC` (and
  optional METHFESSEL-PAXTON note) under `tools/cp2k_oracle/inputs/`.
- Absolute \(E_{\mathrm{tot}}\) may differ (GTH-GPW vs NC-UPF PW); pin
  **occupation algorithm** (unit tests) hard; e2e pins **electron count**,
  **fractional occupations**, **finite Fermi level**, and ATV \(E_{\mathrm{tot}}\)
  against the offline-recorded value for the **identical** ATV settings
  after the first successful offline run is cross-checked for order of
  magnitude vs CP2K free energy / literature metallic Al SCF (documented in
  `REFERENCE.md`). Prefer matching ATV free energy to the recorded pin with
  tight tol once the implementation is frozen.

## Files to create or modify

- `.claude/specs/pw-smearing-al/` (this + acceptance)
- `include/atomiverse/cpu/planewave/pw_solver.h` — members + `PWResult` fields
- `include/atomiverse/cpu/planewave/internal/smearing.h` (new)
- `src/cpu/planewave/internal/smearing.cpp` (new)
- `include|src/.../density_builder.*` — occupation-weighted ρ
- `src/cpu/planewave/pw_solver.cpp` — SCF wiring
- `src/cpu/planewave/internal/pw_forces.cpp` — honor method smearing when set
- `CMakeLists.txt` — add `smearing.cpp`
- `tests/cpu/planewave/test_pw_smearing.cpp` (unit)
- `tests/CMakeLists.txt` — register unit test
- `regressions/pw/pw-smearing-al/` — e2e + `REFERENCE.md`
- `regressions/CMakeLists.txt` — register regression
- `tests/data/pw/Al_pz_nc.upf` — NC Al UPF
- `tools/cp2k_oracle/inputs/al_fermi_dirac.inp` — offline twin

## Tasks

- [ ] Write unit tests for FD / Gaussian / MP-1 occupations + Fermi search
- [ ] Implement `internal::smearing` + wire density / band energy / free energy
- [ ] Extend `PW` / `PWResult` public members (no Config bag)
- [ ] Offline CP2K twin input + `REFERENCE.md` pin notes
- [ ] E2E regression `regressions/pw/pw-smearing-al`
- [ ] Build + unit + regression green; record perf vs CP2K if wall available

## Testing strategy

- **Unit** (`tests/` only): closed-form \(f(x)\) at known \(x\); Fermi bisection
  recovers \(N_e\) on a synthetic spectrum; default `none` path unchanged API.
- **Regression**: Al FCC NC-UPF FD smear; assert convergence, \(|N_{\mathrm{occ}}-3|<\mathrm{tol}\),
  \(E\) within pin tol, occupations not all integers near \(E_F\).

## Out of scope

- Spin-polarized metallic UKS (`pw-spin-polarized-o2`, `pw-uks-kmesh-fe`)
- Cold smearing / Marzari–Vanderbilt
- GPW smearing (separate if needed)
- CUDA PW smearing
- Analytic forces correctness under smear (forces may inherit path; not pinned)

---

**Librarian note:** `cpu/planewave/`; roadmap “Smearing for metals” ○ → ✓;
confidence high on FD/MP math, medium on absolute E pin vs CP2K GTH twin.

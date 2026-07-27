# const-pot-il-edl-box — planar const-pot + non-polarizable IL EDL production box

Queue row (metalwalls-method-loop, wave 5):
`{"slug":"const-pot-il-edl-box","method":"Planar const-pot + non-polarizable IL
EDL production box","geom":"graphene-like walls + IL","observable":"Q, Psi(z),
short NVT","oracle":"metalwalls"}`

Map slice: **E9 scenarios** — §9 "Planar graphene + const-pot + IL" is `~`
(E1 Components ✓, production EDL box open; aqueous twin
`reg_const-pot-tip4p-edl` landed wave-4). Intentionally out: nanoporous
geometry (own row), polarizable IL / PIM-NPT (own row), 2D-PBC product stack,
CUDA, PLUMED.

## What lands

1. **Production EDL regression** `regressions/md/const-pot-il-edl-box/` —
   out-of-library consumer of public headers only:
   - Honeycomb ("graphene-like") single-layer walls, 4-atom rectangular cell,
     a_CC = 1.42 Å, 4×7 cells → 112 C atoms per wall, walls at z = 12 Å
     (type 2, Ψ = −1 V) and z = 28 Å (type 3, Ψ = +1 V), frozen
     (`FixFreeze`). Interior IL gap 16 Å < periodic return path 24 Å so the
     EDL response dominates.
   - Non-polarizable single-site coarse-grained IL: 12 cation (+1 e) +
     12 anion (−1 e) sites, IL-family LJ parameters (BMIM+/PF6−-like sizes
     and masses), per-type `PairLJCut` (LB ion–ion; explicit ion–wall
     entries; no intra-wall LJ), deterministic checkerboard start with a
     nascent-EDL polarization (one charge swap per outer layer).
   - Full E1/E2 pipeline: `cpu::md::ConstPotBO` (η = 1.805132 Å⁻¹ MW
     graphite width, CG solve + ASPC order 2, electrode-only neutrality) at
     `pre_force`, `CoulLong` + `KSpacePME` forces, `NeighborList`,
     `VelocityVerlet` + `NoseHoover` short NVT at 400 K, `set_threads(4)`.
   - Observables: per-electrode Q (`ElectrodeChargeObserver`), Ψ(z)
     (`PotentialProfile`), per-species ρ_s(z) (new observer below).
2. **New Component** `cpu::SpeciesDensityProfile`
   (`include/atomiverse/cpu/observe/species_density_profile.h` +
   `src/cpu/observe/species_density_profile.cpp`) — Observer kind,
   Frame-immutable; bins owned atoms of each requested `kAtom_type` along z
   into instantaneous and run-mean number-density profiles [Å⁻³].
   Ctor args only (`types, n_bins, z_min, z_max`) — no Config bag. This is
   the EDL layering observable of the production box (map §8 has no
   per-species profile today; `PotentialProfile` is charge-only).
3. **Unit tests** `tests/cpu/observe/test_species_density_profile.cpp`
   (single-function tests, seconds).
4. **Offline oracle pin** — no MetalWalls binary on this host, so the
   MW-form algebra is pinned by an independent numpy implementation:
   `tools/metalwalls_oracle/inputs/metalwalls/const-pot-il-edl-box/pin_il_edl.py`
   solves the same Gaussian-electrode BO system (electrode-only neutrality
   projection, `ConstPotBO` kernels, min-image) for the
   exact deterministic initial geometry and prints Q_left/Q_right for
   (a) the empty capacitor (IL charges zeroed) and (b) the full IL box.
   Pins hard-coded in the regression; notes in
   `tools/metalwalls_oracle/pin_notes_const-pot-il-edl-box.md`.

## Acceptance criteria

- [x] `cpu::SpeciesDensityProfile` compiles into `atomiverse`, Observer role
      (FrameImmutable|Deferred), throws on empty type list / bad bins /
      inverted window; density normalization: Σ_bins ρ_s(z)·A·dz = N_s when
      all atoms of species s are inside the window.
- [x] Unit tests pass: ctor validation, instantaneous binning on a known
      config, mean accumulation over ≥2 launches, normalization.
- [x] Regression `reg_const-pot-il-edl-box` (≤ 60 s, public API only; ~2 s
      measured):
  - initial BO solve matches the independent numpy oracle pin for both
    Q_left and Q_right (empty and with-IL solves) to ≤ 1e-6 e (CG tol
    1e-12 vs direct solve; agreement in practice to 12 digits);
  - screening: |Q| with IL > |Q| empty at the same ΔΨ (EDL capacitance
    enhancement over vacuum), from the pinned deterministic solves;
  - polarity: Q_left < 0 < Q_right at Ψ = ∓1 V; |ΣQ_elec| ≤ 1e-8;
  - short NVT (600 × 1 fs): finite E/K/q every chunk, wall atoms frozen,
    electrode neutrality ≤ 1e-6 throughout, polarity preserved at end;
  - Ψ(z) via `PotentialProfile` after the run: wall-to-wall potential
    difference has the applied sign (Ψ(z_R) − Ψ(z_L) > 0);
  - ρ_s(z) via `SpeciesDensityProfile` run-mean: normalization exact
    (12 + 12 ions inside the window), profiles finite and ≥ 0; layering
    diagnostics (first-layer cation/anion asymmetry per electrode) printed
    (measured: cat 5.0 vs an 2.4 at the −1 V wall and mirrored).
- [x] REFERENCE.md states `oracle: metalwalls` (MW-form algebra, offline
      numpy pin; no MW binary), oracle_version, cited model provenance
      (Scalfi/Reed BO form; MW JOSS 10.21105/joss.02373; CG IL family
      Roy–Maroncelli 2010 / Merlet 2012), units.
- [x] Registered: `atv_add_regression(md const-pot-il-edl-box)` and
      `atv_add_cpp_test(... test_species_density_profile ...)`.
- [x] MD-only build green: `atomiverse`, `test_cpu_observer_species_density_profile`,
      `reg_const-pot-il-edl-box` with HF/PW/GPW/CUDA/TORCH **OFF**;
      ctest filtered to the new unit test green; `run_reg_const-pot-il-edl-box`
      green.
- [x] perf: `perf_exempt=true` — MW-only physics (Gaussian const-pot
      electrode + IL EDL box), no MetalWalls binary installed, no LAMMPS
      ELECTRODE twin available.

## Non-goals

- No new Driver, no `*Config` bag, no factory functions.
- No 2D Ewald / slab-correction product work (Yeh–Berkowitz foothold exists).
- No polarizable IL (PIM) here; no nanoporous geometry.
- No claim of quantitative MW capacitance parity — pins are the BO algebra
  on the exact geometry plus physics invariants.

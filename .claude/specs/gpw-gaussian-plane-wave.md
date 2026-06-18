---
title: GPW (Gaussian and Plane Waves) DFT — CP2K-parity Quickstep path
slug: gpw-gaussian-plane-wave
status: draft
science_required: true
acceptance: gpw-gaussian-plane-wave.acceptance.md
---

## Summary

Implement the **GPW method** (Gaussian basis + density collocated onto a uniform
real-space FFT grid + FFT Poisson for the Hartree term + **analytic GTH
pseudopotentials**) as a new inner solver `cpu::gpw::GPWSolver` and an Evaluator
component `cpu::gpw::GPWEvaluator`. The goal is to reproduce **CP2K's Quickstep
GPW** results closely enough to make CP2K's getting-started tutorial(s) pass, then
migrate the same tutorials + a benchmark subset into Atomiverse under
`benchmarks/gpw/` following the existing `pw/si/validation` template.

GPW is the *third* electronic-structure path, sitting between the two we already
have: it keeps the Kohn–Sham problem in the **contracted Gaussian AO basis** (like
`cpu::scf`) but computes the Hartree + XC + long-range local-PP terms on a
**uniform plane-wave grid** (like `cpu::pw`). It is the method behind ~76 % of
CP2K's regression tests, so it is the only Atomiverse path that can target CP2K's
own reference numbers rather than PySCF (our Gaussian SCF target) or QE/ABINIT
(our plane-wave target).

## Domain basis (and the honest accuracy target)

CP2K Quickstep solves `H_KS C = S C ε` in a contracted GTO basis. The KS matrix is
split into an **analytic ("hard") part** and a **grid ("soft") part**:

```
H_KS = T            (kinetic, analytic GTO)
     + V_PP,nl      (GTH nonlocal separable projectors — analytic GTO; GTH
                      projectors are Gaussian×polynomial, so ⟨φ_μ|p_i⟩ is analytic)
     + V_PP,loc,sr  (GTH local short-range, analytic GTO Gaussian potential)
     + V_soft[ρ]    (Hartree + long-range local PP + XC, collocated/integrated
                      on a uniform FFT grid)
```

Density and potential move between representations exactly as CP2K's GPW flowchart:
`ρ(Gaussian) --collocate--> ρ(grid) --FFT--> ρ(G) --Poisson--> V_H(G) --FFT⁻¹-->
V_H(grid) --integrate--> V_μν(Gaussian)`. The total energy uses CP2K's GPW
electrostatic decomposition with a Gaussian core density `n_c` so the long-range
local PP and Hartree combine without a `G=0` divergence.

**"Completely identical" is bounded by physics, not effort.** Bit-for-bit 1e-13
agreement (CP2K's *self*-comparison regtest tolerance) is not reachable by an
independent code — different FFT, different collocation truncation, different
integration order all enter at ~1e-8…1e-6 Ha. The realistic, defensible target,
**provided we ingest CP2K's exact inputs** (same GTH parameters, same MOLOPT basis
exponents/coefficients, same `CUTOFF`/`REL_CUTOFF` grid, same libxc functional and
version), is:

| Quantity | Acceptance target | Rationale |
|---|---|---|
| Total energy, matched grid+basis+PP+XC | **≤ 1e-6 Ha** vs CP2K | collocation/integration discretization floor |
| Total energy, converged (tight CUTOFF) | ≤ 1e-5 Ha cross-code | what "agreement" means between real codes |
| Forces (per component) | ≤ 1e-4 Ha/bohr | finite-grid force noise |

To reach the tight end the XC functional must match CP2K's to machine precision.
Rather than link libxc at runtime, this is delivered by the companion spec
[[native-xc-libxc-parity]] — a native `cpu/xc/` library validated bit-for-bit
against the libxc oracle (PBE etc.) and against CP2K's own formula (PADE). See D7.

## Design

### Reuse map (what already exists — do NOT rebuild)

Directly reusable, verified present in the codebase:

- **Gaussian basis + one-electron integrals**: `cpu/scf/basis_set.h`
  (`build_shell_list`), `cpu/scf/one_electron.cpp` (S, T, nuclear-attraction
  primitives), `cpu/scf/basis_eval.h` (`poly_val/grad` — GTO values on grid
  points), Boys/Rys machinery (`cpu/scf/rys_roots.cpp`). The analytic GTH integrals
  reuse these primitives.
- **FFT + reciprocal space**: `cpu/pw/internal/fft.h` (`FFTPlan3D` r2c/c2r/c2c,
  `FFTWBuffer`), `cpu/pw/internal/g_vector.h` (`enumerate_g_vectors`),
  `cpu/pw/internal/poisson.h` (periodic `V_H(G)=4π ρ(G)/G²`).
- **XC**: `cpu/scf/xc_functional.h` (`evaluate_xc`, `evaluate_xc_batch`, `is_gga`,
  LSDA + PBE), reused pointwise on the uniform grid.
- **Linear algebra / SCF machinery**: `cpu/scf/internal/eigen_solver.cpp`
  (`eigh` generalized via `X=S^{-1/2}`), `cpu/scf/internal/diis_mixer.cpp`,
  density-matrix build `cpu/scf/internal/density_builder.cpp`, convergence
  `cpu/scf/internal/dual_convergence_check.cpp`.
- **Component/Driver framework**: `cpu/component/component.h` (`ComponentBase`),
  `cpu/driver/evaluate.h` + `cpu/driver/single_point.h`, `Frame`
  (`frame.scalars.energy` in eV, `frame.atoms.f{x,y,z}` in eV/Å). `GPWEvaluator`
  mirrors `cpu/pw/pw_evaluator.h` and `cpu/scf/scf_evaluator.h` 1:1.
- **CMake**: `fftw3_stub` target, `ATV_BUILD_RHF`/`ATV_BUILD_PW` gating pattern.

### New work (the actual deliverable)

- **D1 — CP2K data ingestion** (`cpu/gpw/internal/cp2k_io.*`):
  - GTH potential parser (one shared format across `data/POTENTIAL` and
    `data/GTH_POTENTIALS`) → analytic GTH params (`rloc`, `C_i` local;
    `r_l`, `h^l_ij` nonlocal per ℓ, `Zion`). M1 needs `GTH-PADE-q8` (Ar),
    `GTH-PADE-q1` (H) from `data/POTENTIAL`; M3 needs `GTH-PBE-q6/q1` from
    `data/GTH_POTENTIALS`.
  - GTH basis parser (one shared CP2K basis-file format) → contracted GTOs.
    M1 reads `DZVP-GTH-PADE` from `data/BASIS_SET`; M3 reads `DZVP-MOLOPT-GTH`
    (shared-exponent MOLOPT layout) from `data/BASIS_MOLOPT`. The two files use
    the same grammar, different contraction layouts — parse both.
  - Source files are read from CP2K's `data/` dir (path via config or
    `ATV_GPW_DATA_DIR`). Start by supporting only the elements the milestone
    inputs need (Ar, H, O), generalize later.
- **D2 — Uniform-grid density collocation + potential integration**
  (`cpu/gpw/internal/collocate.*`) — *the core new algorithm*. Map Gaussian
  products `φ_μ(r)φ_ν(r)` onto the uniform FFT grid and the adjoint integrate
  `V_μν = ∫ φ_μ(r) V_soft(r) φ_ν(r) dr`. Implement CP2K's `CUTOFF`/`REL_CUTOFF`
  **multigrid** assignment (each Gaussian to a grid level by exponent) and the
  per-Gaussian real-space cutoff radius. NOTE: this is *not* the existing
  atom-centered Becke/Lebedev grid in `cpu/scf/grid.h`; it is a new uniform grid
  built on the `cpu/pw` FFT infrastructure. This is what makes the numbers match
  CP2K.
- **D3 — Analytic GTH integrals in the Gaussian basis**
  (`cpu/gpw/internal/gth_integrals.*`): local short-range Gaussian-potential
  matrix elements and separable nonlocal projector matrix elements `Σ_ij h_ij
  ⟨φ_μ|p_i⟩⟨p_j|φ_ν⟩`. Built on the existing one-electron/Boys primitives.
- **D4 — Martyna–Tuckerman non-periodic Poisson**
  (`cpu/gpw/internal/mt_poisson.*`): reciprocal-space MT decoupling for isolated
  systems (`PSOLVER MT`, `PERIODIC NONE`). Required for the H₂O tutorial.
  Periodic Poisson is reused from `cpu/pw`.
- **D5 — GPW SCF orchestrator** (`cpu/gpw/gpw_solver.*` + `cpu/gpw/internal/*`):
  follow the codebase's actual design — each ES method is a **self-contained domain
  that mirrors its siblings** (`cpu::pw` has its own `internal/{convergence,
  density_builder,...}`; `cpu::scf` has its own; neither reaches into the other's
  `internal`). So GPW gets its OWN thin loop helpers under `cpu::gpw::internal`
  (mixing, convergence, density-from-MO, and a generalized-eigensolve wrapper
  calling xtensor-blas directly) — structurally parallel to `cpu::scf::SCFSolver`,
  NOT reaching into `scf::internal`. The only cross-domain reuse is **shared
  physics via public headers**: the XC functional — `cpu::xc::pade` (and
  `cpu::scf::evaluate_xc` for PBE), exactly as `cpu::pw` reuses
  `cpu::scf::evaluate_xc`. KS matrix is built from the GPW grid/analytic pieces
  already done: `one_electron` (T,S), `gth_integrals` (V_loc^sr, V_nl),
  `hartree_grid` (J), `xc_grid` (V_xc), plus the GTH long-range core-charge
  electrostatics. (The plane-wave `cpu::pw` solver is a different method — Davidson
  on plane-wave wavefunctions — and is not reused: GPW's eigenproblem is the small
  dense generalized one in the AO basis, and CP2K=GPW≠pure-PW.)
- **D6 — GPWEvaluator + GPWConfig** (`cpu/gpw/gpw_evaluator.*`,
  `cpu/gpw/gpw_config.h`): `role()=Evaluator`,
  `capabilities()=Energy|Force|DensityMatrix|Orbitals`; plugs into `Evaluate`,
  `md::MD`, `minimizer::Minimize` unchanged.
- **D7 — Native XC with libxc/CP2K bit-parity** — delegated to the companion spec
  [[native-xc-libxc-parity]]. GPW consumes the shared `cpu/xc/` library through the
  existing `evaluate_xc` interface; no runtime libxc dependency is shipped. PADE
  (CP2K-native, for M1) and PBE (libxc-routed, for M3) are the two functionals GPW
  needs first. This removes the XC-mismatch divergence between Atomiverse-GPW and
  CP2K to machine precision.
- **D8 — Analytic GPW forces** (`cpu/gpw/internal/gpw_forces.*`, milestone M4):
  Hellmann–Feynman (local PP) + Pulay (basis-derivative) + nonlocal-projector +
  grid (Hartree/XC) force terms. Enables geom-opt and MD tutorials.

### Milestones (each gated on a concrete CP2K reference for exact-match)

- **M1 — Energy, periodic, LDA.** Target `tests/QS/regtest-gpw-1/Ar.inp`
  (ref `E = -21.04944231395054` Ha) and `H2.inp` (ref `-1.06345574964568` Ha) from
  the CP2K repo. Exercises D1+D2+D3+D5 with periodic Poisson + LSDA (PADE), no MT,
  no GGA. **First numeric parity gate.**
- **M2 — GGA.** Add PBE gradient-on-grid (`σ=|∇ρ|²` via FFT of `ρ(G)`); match a
  periodic PBE GPW regtest.
- **M3 — Non-periodic + MT → the tutorial.** Add D4; **reproduce the
  getting-started `docs/getting-started/h2o.inp`** (isolated H₂O, PBE,
  `DZVP-MOLOPT-GTH`, `GTH-PBE-q6/q1`, `CUTOFF 400 REL_CUTOFF 50`, `PSOLVER MT`).
  This is the headline "入门教程跑通".
- **M4 — Forces.** D8; validate forces vs CP2K and vs Atomiverse's existing
  finite-difference oracle (`cpu/scf/gradient.h::scf_finite_difference_gradient`
  pattern). The geom-opt and MD "Next Steps" tutorials then run on `Minimize`/`MD`
  with `GPWEvaluator`.
- **M5 — Migration.** Add `examples/09_gpw.cpp` (H₂O single point) and migrate the
  tutorial + a benchmark subset into `benchmarks/gpw/` following
  `benchmarks/pw/si/validation` (CP2K deck template + reference CSV + numeric
  comparison harness + `SOURCES.md`).

## Files

New (under the new `cpu/gpw` module, mirroring `cpu/pw` layout):
- `include/atomiverse/cpu/gpw/{gpw_config.h,gpw_solver.h,gpw_evaluator.h}`
- `include/atomiverse/cpu/gpw/internal/{cp2k_io,collocate,gth_integrals,mt_poisson,gpw_forces}.h`
- `src/cpu/gpw/{gpw_solver.cpp,gpw_evaluator.cpp}`
- `src/cpu/gpw/internal/{cp2k_io,collocate,gth_integrals,mt_poisson,gpw_forces}.cpp`
- `examples/09_gpw.cpp`
- `benchmarks/gpw/h2o/validation/` (+ `reference/cp2k/`, `SOURCES.md`) — copy/adapt
  `benchmarks/pw/si/validation`
- `tests/cpu/gpw/test_cpu_gpw_*.cpp`

Edit:
- `CMakeLists.txt` — add `option(ATV_BUILD_GPW ...)` (requires `ATV_BUILD_RHF` +
  `ATV_BUILD_PW`), wire `src/cpu/gpw/*.cpp`, define `ATV_HAS_GPW`; optional
  `ATV_GPW_RUNTIME_LIBXC`.
- `examples/CMakeLists.txt`, `benchmarks/CMakeLists.txt`, `tests/CMakeLists.txt`
  — register new targets behind `ATV_BUILD_GPW`.
- `.claude/architecture.md` — add the `cpu/gpw` module + "third ES path" note.

## Tasks

> Runs **in parallel** with [[native-xc-libxc-parity]] (the two are decoupled).
> Phases 0–2 below need no XC and proceed independently; the XC dependency first
> bites at Phase 3 (M1 needs PADE) and Phase 5/6 (PBE), by which time XC tracks
> X2/X3 should be landing.

- [~] **Phase 0 (D1):** GTH potential + MOLOPT basis parsers; unit-test against
      hand-checked Ar/H/O entries from CP2K `data/`.
      - [x] GTH potential parser (`cpu/gpw/internal/cp2k_io`) — H GTH-PBE-q1,
            O GTH-PBE-q6, Ar GTH-PADE-q8 verified via Slurm (`test_cpu_gpw_cp2k_io`).
      - [x] GTH/MOLOPT Gaussian basis parser (`cpu/gpw/internal/cp2k_basis`) —
            DZVP-GTH-PADE + DZVP-MOLOPT-GTH verified (`test_cpu_gpw_basis`).
- [~] **Phase 1 (D2):** uniform-grid collocation + integration
      (`cpu/gpw/internal/collocate`) — single grid validated: `integrate(V≡1)`
      recovers analytic overlap S and `∫ρ = Tr(P S)` to 1e-4
      (`test_cpu_gpw_collocate`). TODO: CUTOFF/REL_CUTOFF multigrid, periodic image
      sums, FFT Poisson on this grid.
- [~] **Phase 2 (D3):** analytic GTH local-sr + nonlocal projector integrals
      (`cpu/gpw/internal/gth_integrals`) — primitive Cartesian Gaussians, s/p,
      closed-form validated (`test_cpu_gpw_gth_integrals`). TODO: contraction over
      shells, l>=2 projectors, and CP2K sign/normalization-convention parity
      (deferred to M1 validation).
- [~] **Phase 3 (D5 + periodic Poisson):** GPW SCF orchestrator
      (`cpu/gpw/gpw_solver`, self-contained in `cpu::gpw`, reuses public
      `cpu::scf::scf_utils` LA). First end-to-end self-consistency green: minimal
      H2 converges (`test_cpu_gpw_solver`). TODO for the **M1 gate** (`Ar.inp`,
      `H2.inp` to ≤1e-6 Ha): real DZVP basis (not single-s), exact GTH params via
      the parsers, erf-vs-point ion-ion / self-energy electrostatic corrections,
      and the absolute-energy bit-parity tuning vs CP2K.
- [ ] **Phase 4 (D6):** `GPWEvaluator`/`GPWConfig`; `cpu::single_point` runs GPW;
      example `09_gpw.cpp` (periodic first).
- [x] **Phase 5 (M2):** PBE/GGA on grid (`grid_fft` FFT grad/div + `xc_grid_pbe`,
      reusing `cpu::scf::evaluate_xc(PBE)`); FD-validated GGA potential.
- [x] **Phase 6 (D4, M3):** MT isolated Poisson (`mt_poisson`, truncated-Coulomb
      Rozzi-0D) + **`h2o.inp` tutorial reproduced to 0.43 mHa**
      (`test_cpu_gpw_h2o_tutorial`: -17.2127672 vs CP2K -17.21233714458516, every
      term sub-mHa). TODO: exact erf-screened MT for bit-parity.
- [ ] **Phase 7 (D7):** runtime libxc bridge; tighten H₂O parity to ≤ 1e-6 Ha.
- [ ] **Phase 8 (D8, M4):** analytic forces; geom-opt + MD tutorials run.
- [ ] **Phase 9 (M5):** migrate tutorial + benchmark into `benchmarks/gpw/`;
      `SOURCES.md`; comparison harness vs CP2K reference CSV.

## Testing

- TDD per Atomiverse convention (RED test referencing its `ac-NNN` before GREEN).
- Unit (`tests/cpu/gpw/`): parsers (D1), collocation round-trip & charge
  conservation (D2), GTH matrix elements (D3), MT Poisson vs analytic point-charge
  (D4).
- Integration / parity: per-milestone, compare Atomiverse total energy against the
  CP2K reference number from the corresponding `.inp` (M1 uses values embedded in
  `tests/QS/regtest-gpw-1/TEST_FILES.toml`; M3 uses CP2K run of `h2o.inp` via the
  locally-installed `cp2k` — see `~/.local/bin/cp2k`).
- Migration benchmark (`benchmarks/gpw/h2o/validation`): regenerate the CP2K
  reference with `run_*.sh`, assert Atomiverse within tolerance, emit CSV.

## Out of scope

- GAPW (all-electron Gaussian-augmented), hybrids/HFX in GPW, MP2/RPA, k-points
  for GPW (Γ-only first), CUDA GPW backend, smearing/metals. Each is a future spec.
- A general CP2K `.inp` parser — only the minimal sections the milestone inputs use
  are read; calculations are still specified by `GPWConfig` + ingested data files.

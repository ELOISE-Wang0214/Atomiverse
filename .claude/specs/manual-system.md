# Atomiverse Manual System — Planning Spec

status: draft
date: 2026-07-01
depends on: [manual-discovery](manual-discovery.md)
purpose: The plan for turning Atomiverse's docs/examples/benchmarks into one
coherent, industrial-grade, C++-centered computational-chemistry manual. Derived
from the real codebase (see manual-discovery), not an imposed template.

---

## 1. Manual vision

> Atomiverse = an industrial-grade molecular-simulation engine with a
> textbook-level C++ manual.

The manual teaches computational chemistry *through the Atomiverse architecture*.
Its through-line is the engine's own unifying idea:

> **One loop, many physics.** A `Driver` owns the loop; a `Component` with a
> `Force`/`Energy` capability supplies the physics; a `Frame` holds the state.
> You change method by swapping the evaluator, and change hardware by swapping
> the backend namespace — the loop, the tags, the result bundle stay the same.

Every chapter answers the six questions: (1) what is the CC concept, (2) how does
Atomiverse represent it, (3) which public C++ API exposes it, (4) what scientific
result to expect, (5) how it is validated, (6) how the same architectural pattern
recurs elsewhere. C++ is primary; Python/nanobind is documented as an adjunct only.

## 2. Terminology canon (single source of truth)

Fix these names everywhere; a new chapter inherits whatever is written first.

| concept | canonical name | NOT |
|---|---|---|
| loop owner | `Driver` / category base (`EvaluateDriver`, `DynamicsDriver`, `OptimizeDriver`, `ResponseDriver`) | "SinglePoint", "SCF driver" |
| result bundle | `SystemScalars` (`frame.scalars`) | "output struct" |
| radial distribution | `ATV::RDF` | `cpu::RDF`, `cpu::observer::RDF` |
| XYZ recorder | `cpu::XYZ` | `cpu::recorder::XYZ` |
| tag type | `Capability` (singular) | `Capabilities` |
| Gaussian SCF engine | `cpu::scf::SCFSolver` (inner) + `SCFEvaluator` (component) | `scf::SCF` driver, `scf::Phase` |
| precision | tiered `core/precision.h` (`StorageReal`/`WorkReal`/`AccReal`, `LocalInt`/`TagInt`/`ImageInt`/`BigInt`) | global `-DATV_REAL` |
| ML/torch layer | "differentiable-potential library" | "third backend" |

Header locations are canonical too: driver headers in `cpu/driver/`, Frame/Box in
`cpu/frame/`.

## 3. Tutorial roadmap

Each chapter maps to a *real* capability (status from manual-discovery §5).
"exists" = a current `docs/tutorial/*`; "new" = to be written; "needs code" =
blocked on a missing feature.

| # | chapter | concept | key public API | status of capability | doc state | validation / benchmark |
|---|---|---|---|---|---|---|
| **I. Foundations** |
| 0 | The Atomiverse architecture | Driver/Component/Frame, 2 axes | `cpu::Evaluate`, `ComponentBase` | supported | exists (00) — refresh terminology | — |
| 1 | Systems, units, precision, validation | Box/Atom/units, tiered precision, FD-gradient check | `cpu::Box`, `core/units.h`, `core/precision.h`, `scf_finite_difference_gradient` | supported | **new** | FD gradient benchmark |
| **II. Localized-basis electronic structure** |
| 2 | AOs → single-point RHF | Gaussians, integrals, MOs, SCF, RHF | `build_integrals`, `rhf_solve`, `SCFEvaluator`, `single_point` | supported | exists (01,02) | qc/rhf vs PySCF |
| 3 | Open-shell & spin | UHF, spin contamination, Mulliken | `uhf_solve`, `mulliken_charge/spin` | supported | exists (03) | vs PySCF |
| 4 | Basis sets & convergence | basis-set effect, STO-3G→def2-SVP | `BasisSet`, `SCFConfig` | supported (def2-TZVP partial) | **new** | basis-convergence benchmark |
| **III. DFT** |
| 5 | Kohn–Sham DFT | density rep, KS flow, XC, grids | `uks_solve`, `XCType`, `GridLevel` | supported | exists (04) — extend | vs PySCF; native-XC vs libxc 1e-12 |
| 6 | Hybrid functionals & gradients | exact-exchange admixture, analytic ∇ | `b3lyp_config`/`pbe0_config`, `scf_gradient` | supported | **new** | hybrid ∇ vs PySCF |
| **IV. Correlation methods (status chapter)** |
| 7 | Beyond mean field | MP2/CI/CC hierarchy, cost vs accuracy | — | **absent (planned)** | **new (honest status)** | — (documents intent only) |
| **V. Periodic / plane-wave** |
| 8 | Plane-wave DFT for solids | periodic cells, reciprocal space, cutoff, k-points, pseudopotentials | `cpu::pw::PWEvaluator`, `PWConfig` | experimental (LDA, Γ-validated) | exists (05) — add caveats + convergence | pw/si EOS vs CP2K/ABINIT |
| 9 | Gaussian-and-plane-waves (GPW) | Gaussian AO + grid Hartree/XC, GTH PP | `cpu::gpw::GPWEvaluator` | experimental (closed-shell, Γ) | **new** (examples 09/10 exist) | **new** benchmarks/gpw vs CP2K |
| **VI. Classical simulation** |
| 10 | Force fields & classical MD | pair/bonded, neighbor lists, PBC, NVE | `cpu::md::MD`, `set_pair`, `NeighborList`, `VelocityVerlet` | supported | exists (07) | md/lj vs LAMMPS + energy conservation |
| 11 | Ensembles: thermostats & barostats | NVT/NPT, Nosé–Hoover/Langevin/Berendsen | `set_thermostat`, `set_barostat`, `NoseHooverNPT` | supported | **new** | ensemble-means vs LAMMPS |
| 12 | Geometry optimization | LBFGS/CG on a PES | `cpu::minimizer::Minimize`, `LBFGS`/`CG` | supported | exists (08) | convergence + final-geometry |
| **VII. Analysis & observables** |
| 13 | Observables & trajectory analysis | energy conservation, RDF, **MSD/diffusion/VACF** | `ATV::RDF`, (new observers), `cpu::XYZ` | **partial (needs code)** | **new** | observable-convergence + reference |
| **VIII. Response & cross-method** |
| 14 | Vibrational analysis & thermochemistry | FD Hessian, frequencies, thermo | `cpu::vib::Hessian`, `thermochem.h` | supported | **new** (undocumented) | vs PySCF |
| 15 | Cross-method workflows | swap the evaluator; cost vs accuracy; backends | classical MD → AIMD → observable | supported | **new** | end-to-end golden run |

Appendices: **Backend guide** (cpu reference / cuda partial / torch adjunct),
**Validation policy** (§5), **Benchmark guide** (§6), **Method-status matrix**
(manual-discovery §5).

## 4. Example strategy

Examples are executable manual components. Each must have: scientific purpose,
C++-usage purpose, documented expected result, runtime category
(fast/slow/hardware), backend requirement, a manual-chapter link, and a
validation/regression link. Rules: public API only; no private internals (except
explicit developer chapters); no hidden/developer-specific paths (enforce the
§7 portability rule). Current examples 00–08 are compliant; 09/10 need index
entries + de-hardcoded data paths; new chapters (4,6,9,11,13,14,15) need new
examples.

## 5. Validation policy

Every result claim is classified, and the class dictates the check:

| class | meaning | mechanism | example tolerance |
|---|---|---|---|
| exact | must be bit-reproducible | golden file / `static_assert` | integer counts, phase order |
| numerical tolerance | float round-off only | abs/rel compare | native XC vs libxc **1e-12**; RHF energy vs PySCF **1e-7**; analytic ∇ vs PySCF **1e-5 Ha/Bohr** |
| physical approximation | method/basis error | compare to converged/reference | pw Si EOS a0 ≤0.08 Å, B0 ≤25 GPa |
| stochastic variation | RNG-seeded dynamics | seed control + ensemble means | Langevin ⟨T⟩ vs LAMMPS ≤2% |
| backend-dependent | cpu vs cuda mixed precision | per-quantity tolerance, documented | see §6 |
| performance noise | timing variance | median-of-N + regression band | perf CSV thresholds |

Requirements: deterministic seeds for all stochastic examples; reference data
committed with a generator + provenance (`SOURCES.md`); reference regeneration
must target the committed path (fix the qc-generator drift); precision defaults
documented (accumulators always double). No result is reported without its class.

## 6. Backend-consistency policy

- **CPU is the reference.** Every physical quantity has a CPU golden value.
- **CUDA is compared, not assumed identical.** Because CUDA runs mixed precision
  (`real_t=float` storage), cpu↔cuda is **not** bitwise; the policy defines a
  per-quantity tolerance (energies/forces to the storage-precision floor) and a
  committed `benchmarks/*/backend/` parity case. **A cpu↔cuda CI job is
  currently missing and is a required deliverable.**
- **Torch is an adjunct**, validated by gradient checks (autodiff vs FD), not by
  driver-level parity — it has no driver.
- Mixed-precision on CPU (new tiered `StorageReal=float`) gets the same treatment:
  a documented tolerance vs the all-double reference + an energy-drift check.

## 7. Industrial-grade rules (enforced in cleanup + CI)

- No developer-machine-only assumptions: no hardcoded absolute paths (2 fixed;
  audit the rest), data-dependent tests self-skip when data is absent.
- Deterministic examples where physics allows; documented expected output.
- Stable public API in tutorials; private internals only in developer chapters.
- Errors surface as catchable C++ exceptions, never process aborts (see the
  molrs `xyz_read_first_frame` abort — fix upstream).
- Fast checks in CI; expensive/GPU/nightly checks clearly marked.

## 8. C++ tutorial style guide

- Construct systems in C++ (Frame from molrs I/O or builders); show units.
- Use categorical setters (`set_pair`/`set_integrator`/…) as the primary API;
  show `add()`/`on(phase).add()` only for advanced hooks.
- One concept per snippet; every snippet compiles against the real public header.
- End each with the expected numerical result + its validation class (§5).
- Cross-reference the recurring architectural pattern (question 6).

## 9. Migration plan (what moves / stays / deprecates)

- **Stays:** tutorials 00–08 (accurate); examples 00–08; benchmark layout
  `<domain>/<case>/{validation,performance}`; Doxygen setup.
- **Fix in place (cleanup, Task #5):** the manual-discovery §11 fix list
  (namespaces, developer-doc paths, SCF-section rewrite, getting-started row,
  GPW indexing, libtorch quarantine, small slips).
- **Add:** new chapters 1,4,6,7,9,11,13,14,15 + their examples; GPW/kspace/vib/
  minimizer benchmark cases; cpu↔cuda backend-consistency case + CI job;
  perf-regression gate; register `test_xc_vs_libxc`; commit or delete the
  `md/pe` MSU reference data.
- **Deprecate/reframe:** `libtorch-backend.md` → "differentiable potentials"
  adjunct, not a backend; the toy pw `V_loc` path documented as non-physical.

## 10. Missing-feature roadmap (gates specific tutorials)

Implement the *smallest correct* version, only when a chapter needs it. **Policy
([[feedback-infrastructure-in-molrs]]): infrastructure — I/O, data models,
trajectory analysis — is added in molrs (Rust) + exposed via FFI, NEVER
reimplemented in C++. Only physics/compute kernels (evaluators, integrators,
minimizers, XC) are added in the C++/CUDA backends.**

Discovery correction: the analysis primitives below are **already implemented in
molrs** (`molrs/src/compute/`, a full freud-class library behind the `Compute`/
`Fit` traits); the "feature" is FFI exposure + a thin C++ observer wrapper, not a
new implementation. See [molrs-analysis-ffi](molrs-analysis-ffi.md).

| feature | unlocks | where it lives / minimal work | validation | priority |
|---|---|---|---|---|
| **MSD + Einstein self-diffusion** | Ch 13 | molrs `compute::msd` + `fit::EinsteinDiffusion` → FFI + `cpu::observer::MSD` wrapper | vs known-D reference + molrs Rust test | **high** |
| **VACF + Green-Kubo D** | Ch 13 | molrs `fit::VACF`/`GreenKuboDiffusion` → FFI wrapper | consistency with Einstein D | medium |
| **RDF via molrs** (replace hand-rolled C++ RDF) | Ch 13 | molrs `compute::rdf` → FFI; parity vs C++ RDF | parity, then deprecate C++ RDF | medium |
| **Static structure factor S(k)** | Ch 13 | molrs `compute::diffraction` → FFI wrapper | vs RDF transform | low |
| DCD/TRR/PDB trajectory I/O | Ch 10/13 | already in molrs → FFI readers | round-trip vs XYZ | low |
| FIRE minimizer | Ch 12 | **C++ Optimizer component** (algorithm/physics stays C++) | vs LBFGS on same PES | low |
| cpu↔cuda parity harness | §6, backend appendix | benchmark case + CI job | per-quantity tol | high (industrial) |
| tiered-precision migration (Frame→StorageReal) | Ch 1 | Task #9 (C++ storage types) | mixed-vs-double drift | medium |

Correlation methods (MP2/CC) are **out of scope**; Ch 7 documents them as planned.

---

## Execution order

1. Cleanup (Task #5): terminology + developer-doc drift + GPW indexing + libtorch
   reframe + reference-data fixes. (Unblocks accurate new chapters.)
2. Analysis observers (Task #8: MSD/diffusion first) — unlocks Ch 13.
3. New chapters on existing capabilities (1,4,6,9,11,14,15) + examples.
4. Benchmark build-out (Task #7): backend-consistency CI, convergence, regression,
   missing-domain cases.
5. Precision migration (Task #9) + Ch 1 precision section.

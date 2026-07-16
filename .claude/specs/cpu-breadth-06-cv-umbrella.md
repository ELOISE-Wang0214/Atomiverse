---
title: CPU breadth 06 — collective variables, harmonic bias, and umbrella sampling
slug: cpu-breadth-06-cv-umbrella
status: approved
execution: queued
last_progress_update: 2026-07-16
science_required: true
created: 2026-07-15
acceptance: cpu-breadth-06-cv-umbrella.acceptance.md
depends_on: cpu-breadth-05-remd-component-scope
---

# CPU breadth 06 — collective variables, harmonic bias, and umbrella sampling

## Summary

Add the smallest useful enhanced-sampling stack on CPU:

1. a public value/gradient CollectiveVariable protocol;
2. one PBC-aware, stable-tag `DistanceCV`;
3. one `HarmonicBias` evaluator usable by Evaluate, MD, Minimize, and MC;
4. multiple fixed-center umbrella windows using the replica executor from link
   05; and
5. one-dimensional WHAM analysis in molrs.

The probe tests cross-workflow evaluator reuse, sparse global state access,
requested Energy versus Force writes, replica parameterization without
exchange, and the engine/analysis ownership boundary. Atomiverse measures and
propagates; molrs bins, normalizes, solves WHAM, and reports the PMF.

## Progress (2026-07-16)

**Execution: queued after link 05 — 0/9 tasks and 0/10 active acceptance
criteria verified. GPU remains deferred.**

No CV, bias, umbrella-window, sample-streaming, or WHAM implementation is
claimed. This link waits for the per-replica program ownership established by
05; every active task remains open.

## Domain basis

For collective variable `xi(R)`, center `xi0`, and force constant `k`,

`U_bias = 1/2 k (xi-xi0)^2`,

`F_bias = -k (xi-xi0) grad_R xi`.

The first CV is the minimum-image distance between two stable atom tags. Its
gradients on the two atoms are equal and opposite. Units are `xi` [Å], `k`
[eV/Å²], bias energy [eV], force [eV/Å].

Umbrella sampling follows Torrie & Valleau, J. Comput. Phys. 23, 187 (1977).
WHAM follows Kumar et al., J. Comput. Chem. 13, 1011 (1992). The public
regression reconstructs the analytic radial PMF of a two-particle harmonic
model, including the radial Jacobian term, over only the connected/overlapped
sampled interval.

## Recon / current evidence

- The roadmap has no CV layer, bias component, umbrella method, or WHAM.
- `IdentityReplication` already associates one scalar parameter with each
  replica; after link 05 it has real per-replica program instances.
- Atomiverse's component capabilities say what a launch may write, but there is
  no general per-launch request distinguishing MC energy-only evaluation from
  force-requiring MD/Minimize/Hessian.
- MC rollback restores the proposal changeset, not arbitrary force fields. A
  bias that writes trial forces during an energy-only request would leave force
  inconsistent after rejection.
- Frame must not become a cache for CV values/gradients or WHAM histograms.
- Project ownership rules place trajectory analysis and statistical reduction
  in molrs, not duplicated C++ math.

## Architecture questions

1. Is a CV a low-level value/gradient protocol rather than a Driver,
   Component taxonomy, or Frame field?
2. Can one bias evaluator obey an explicit requested-capability channel instead
   of inspecting which workflow/phase invoked it?
3. Can umbrella windows reuse replica expansion/program scope without copying
   REMD scheduling or asking users for N Frames/drivers?
4. Can streaming measurements cross to molrs while all histogram/WHAM math
   remains outside the simulation engine?

## Design

### CollectiveVariable is a narrow protocol

Create a cohesive CPU `bias/` domain. `CollectiveVariable` has constructor/
build-time resolution and one allocation-free evaluation contract: return a
scalar value and, when requested, write gradients into caller-owned/pre-sized
storage. It is not a Driver and is not itself a scheduled Component.

`DistanceCV(tag_i, tag_j)` resolves stable tags, applies the canonical Box
minimum image, and refreshes owner/local indices after migration generation
changes. On MPI ranks it communicates only the two tagged positions and scatters
gradient contributions to the owners; no global Frame copy is introduced.
Missing/duplicate tags and zero separation fail with explicit diagnostics.

Users can implement the public protocol outside the library. `HarmonicBias`
owns an explicit CV object through the component-lifetime rules from link 01.

### Requested capabilities, not workflow tests

Add a narrow per-launch requested-capabilities service in the execution
context. The active Driver sets it before evaluator stages:

- MD, Minimize, and Hessian request Energy/Force as required;
- MC trial evaluation requests Energy only;
- Evaluate exposes an explicit request with a documented default.

`HarmonicBias` advertises possible `Energy|Force` output. It always computes
the CV value and bias energy; it asks for the gradient and writes force only
when Force is requested. It never switches on `PhaseKind` or Driver type. In an
energy-only MC trial, rejection restores the proposal and leaves force arrays
untouched.

This link supplies one real caller of the request channel. Link 07 decides
whether existing SCF/PW/other evaluators should adopt it generally.

### Explicit binding across workflows

HarmonicBias declares no workflow-specific default binding. The same component
type is explicitly bound at Evaluate.evaluate, MD.force, Minimize.force,
Hessian.force, or MC.energy. Unit/public conformance proves identical energy and
force math; MC additionally proves exact rollback and energy-only behavior.

### Umbrella windows reuse replica infrastructure

`IdentityReplication(centers)` generates windows. A fixed common temperature,
one `k`, and one exchange-free replica program run independent MD in each
window using link 05's per-replica program scope and executor. It is not REMD
with disabled acceptance hidden in a special branch: the coupling stage is an
ordinary no-exchange window observer/scheduler. Each window RNG derives from a
master seed and window index.

Atomiverse streams records `{window, center, k, temperature, cv_value}` through
the existing analysis boundary. Frame remains the sole simulation state
currency; PMF and diagnostics are analysis results.

### WHAM belongs to molrs

molrs adds a one-dimensional WHAM implementation/streaming accumulator that:

- consumes window metadata and CV samples;
- iterates to a declared free-energy tolerance with a maximum iteration cap;
- returns convergence, iterations, PMF/bin probabilities, effective counts,
  and adjacent-window overlap;
- detects disconnected overlap graphs and empty/effectively empty bins, failing
  loud instead of returning a plausible curve;
- exposes the result through molrs Python and the existing CXX analysis boundary
  only as needed for streaming/finalization.

No histogram, normalization, or WHAM iteration is implemented in Atomiverse
C++.

## Files to create or modify

- `include/atomiverse/cpu/bias/` and `src/cpu/bias/` (new cohesive domain)
- requested-capability service in CPU execution/launch context and Driver stage
  setup
- explicit binding/tests in Evaluate, MD, Minimize, Hessian, and MC
- umbrella window preset/program using `cpu/replica/` infrastructure
- molrs compute/free-energy WHAM implementation, bindings, tests, and docs
- Atomiverse streaming analysis bridge/observer for CV records
- `regressions/sampling/umbrella_harmonic/` (new, with `REFERENCE.md`)
- lifecycle/performance benchmark and `.claude/notes/notes.md`

## Tasks

- [ ] Write failing DistanceCV analytic/finite-difference, PBC, stable-tag,
      migration, missing-tag, zero-distance, and 1/2-rank tests.
- [ ] Write HarmonicBias energy/force finite differences and the four-workflow
      plus MC rollback/requested-write conformance tests.
- [ ] Write synthetic WHAM oracle/disconnected-overlap tests in molrs and the
      public umbrella regression before production implementation.
- [ ] Implement the CV protocol and DistanceCV with caller-owned gradient
      storage and zero launch allocation.
- [ ] Implement the requested-capability channel and HarmonicBias without any
      Driver/PhaseKind branch.
- [ ] Build exchange-free umbrella windows on IdentityReplication and link 05's
      per-replica programs; derive independent deterministic RNG streams.
- [ ] Stream CV samples to molrs and implement all histogram/WHAM math there.
- [ ] Run CPU unit, 1/2-rank CV tests, public umbrella regression, molrs tests,
      allocation, and performance gates.
- [ ] Record the request-channel, CV ownership, and analysis boundary evidence
      in `.claude/notes/notes.md` for link 07.

## Testing strategy

- **Unit**: DistanceCV values/gradients, bias analytic math, requested writes,
  explicit binding across workflows, MC rollback, window/RNG independence, and
  streaming schema.
- **molrs unit/scientific**: synthetic histograms with known PMF, convergence,
  effective counts, overlap graph, disconnected/empty failure paths.
- **Regression**: 4–8 overlapping windows on an analytic two-particle radial
  potential, public API, connected-range PMF RMSE < 0.3 kBT, ≤60 seconds.
- **Performance**: each window tree builds once; CV/bias/window propagation has
  zero steady allocation; sample streaming is batched/off the hot path.

## Deferred — GPU

- Device CV value/gradient protocols and CUDA bias kernels.
- GPU umbrella window programs, streaming, and CPU/GPU PMF parity.
- multi-GPU tags/owner lookup for sparse CV endpoints.

No CUDA/Torch source is modified. CPU replica-parallel MPI umbrella execution
and multidimensional WHAM are ordinary future CPU scope.

## Out of scope

- PLUMED, metadynamics, ABF, steered MD, adaptive windows, replica exchange
  umbrella sampling, or multidimensional CVs.
- Storing CV values/gradients/histograms/PMF in Frame.
- C++ histogram, normalization, WHAM, or plotting code.
- More built-in CVs than stable-tag distance.

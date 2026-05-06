---
title: LAMMPS-style categorical driver API — set_* setters, capability-driven auto-injection
status: in-progress
created: 2026-05-06
---

# LAMMPS-style categorical driver API

## Summary

Reorganize the user-facing surface of the three CPU drivers
(`cpu::Evaluate`, `cpu::md::MD`, `cpu::minimizer::Minimize`) into a
**LAMMPS-shaped, C++-typed, project-namespaced categorical API**:

- One `set_*(value)` method per LAMMPS command kind.
- LAMMPS structure (one pair, one neighbor, one integrator;
  many fixes, observers, records).
- Project-native names — no LAMMPS jargon (`dump` / `compute` /
  `fix`) leaks into Atomiverse code.
- Typed values, no string parsing.

The previous draft's two failures (fluent `with_*` toggles + auto-
deriving NeighborList from `pair_style` cutoff) are corrected:

- `set_*` is the C++ idiom for "configure category X with this typed
  value." Same shape across all driver categories.
- NeighborList stays **explicit and user-set**, exactly like LAMMPS
  (`neighbor 0.3 bin` is its own command). The driver does not
  invent a NeighborList from the pair cutoff.

`add()` and `on(phase).add(c, order)` from `Driver<>` remain
**public** as the documented escape hatch for user-defined custom
hooks (LAMMPS analogue: writing a custom `compute` / `fix` and
issuing `compute` / `fix` against it). The categorical setters are
the canonical surface; raw `add()` / `on()` are for code the
project does not yet have a category for.

Auto-injection is restricted to two narrow cases:

1. Capability-dependency resolution (e.g. NPT declares
   `reads() == Stress` → driver auto-installs `VirialStress` if no
   stress provider is registered).
2. MPI halo exchange when `ATV_BUILD_MPI` is defined and a pair
   evaluator is set (compile-time gated, not runtime).

Everything else the user names explicitly.

## Domain basis

Not applicable — API ergonomics over an existing physics
implementation. Numerics are byte-identical to a correctly-written
manual graph today.

## Design

### Categorical setters (project-native names)

Same shape on every driver, with category-appropriate setters
present or absent. Naming uses `set_*` consistently — singular
categories configure-or-throw on conflict; plural categories append.

| Setter | Cardinality | Accepts | LAMMPS analogue |
|---|---|---|---|
| `set_pair(P)`              | singular | one pair evaluator (LJ, future Coulomb-real-space) | `pair_style …` |
| `set_bond(B)`              | singular | one bond evaluator | `bond_style …` |
| `set_angle(A)`             | singular | one angle evaluator | `angle_style …` |
| `set_dihedral(D)`          | singular | one dihedral evaluator | `dihedral_style …` |
| `set_neighbor(N)`          | singular | one neighbor list (`NeighborList<>{cutoff, skin, every, …}`) | `neighbor` + `neigh_modify` |
| `set_integrator(I)`        | singular | one integrator (`VelocityVerlet`, `NoseHooverNPT`) | `fix … nve|nvt|npt` |
| `set_thermostat(T)`        | append   | thermostat (`NoseHoover`, `FixLangevin`) | `fix … langevin` |
| `set_barostat(B)`          | append   | barostat (`Barostat` Berendsen) | `fix … berendsen` |
| `set_observer(O)`          | append   | derived-quantity observer (`RDF`, `MullikenCharges`, `VirialStress` as output) | `compute …` |
| `set_record(R)`            | append   | recorder (`XYZ`, `Zarr`, thermo file) | `dump …` |
| `set_thermo(int every)`    | singular | thermo print cadence | `thermo N` |
| `set_optimizer(M)`         | singular | optimizer (`LBFGS`, `CG`, `FIRE`); **`Minimize` only** | `min_style …` |

`set_*` semantics:
- **Singular**: `set_pair(LJ{ε,σ,6.8})` twice with the same value is
  idempotent. With a different value it throws `std::invalid_argument`,
  message names both candidates.
- **Append (plural)**: each call adds to the category's ordered list.
  Order preserved (LAMMPS fix-id ordering).

`add()` and `on(phase).add(c, order)` remain **public**. They are the
escape hatch for user-defined hooks the project does not have a
category for. Two examples where a user reaches for them:

```cpp
// Custom compute that doesn't fit set_observer (e.g. writes a
// non-frame field on the side):
sim.on(ATV::md::Phase::end_of_step).add(MyCustomHook{}, +50);

// Pre-register an experimental component that hasn't picked a role yet:
sim.add(ExperimentalThing{});
```

The categorical setters are implemented in terms of `add()` /
`on().add()` internally, so any combination of categorical + raw
calls produces a single coherent Plan.

### Driver-by-driver setter map

| Driver | Setters present |
|---|---|
| `cpu::md::MD` | `set_pair`, `set_bond`, `set_angle`, `set_dihedral`, `set_neighbor`, `set_integrator`, `set_thermostat`, `set_barostat`, `set_observer`, `set_record`, `set_thermo` |
| `cpu::minimizer::Minimize` | `set_pair`, `set_bond`, `set_angle`, `set_dihedral`, `set_neighbor`, `set_optimizer`, `set_observer`, `set_record`, `set_thermo` |
| `cpu::Evaluate` | `set_pair`, `set_bond`, `set_angle`, `set_dihedral`, `set_neighbor`, `set_observer`, `set_record` |

### Concrete example — what user code looks like

Before (current `tests/mpi/test_pe_msu_pipeline.cpp`):

```cpp
sim.on(ATV::md::Phase::pre_force).add(cpu::md::mpi::HaloExchangeForward{lj_cut}, -20);
sim.on(ATV::md::Phase::post_force).add(cpu::md::mpi::HaloExchangeReverse{}, -1000);
sim.on(ATV::md::Phase::force).add(cpu::md::NeighborList<>{lj_cut}, -10);
sim.add(cpu::md::VirialStress{});
sim.add(cpu::md::VelocityVerlet{dt});
sim.add(cpu::md::LJ{eps, sigma, lj_cut});
sim.add(cpu::md::NoseHooverNPT{T, P, tau_T, tau_P});
sim.build();
sim.run(frame, 1000);
```

After:

```cpp
cpu::md::MD sim;
sim.set_pair(cpu::md::LJ{eps, sigma, lj_cut});
sim.set_neighbor(cpu::md::NeighborList<>{lj_cut, /*skin=*/0.3});
sim.set_integrator(cpu::md::NoseHooverNPT{T, P, tau_T, tau_P});  // reads()==Stress → VirialStress auto
sim.set_thermo(100);
sim.build();
sim.run(frame, 1000);
```

`lj_cut` appears twice — once on the pair force field (its own
parameter) and once on the neighbor list (its own parameter), same as
LAMMPS. NPT does not name `VirialStress` — capability resolution
provides it. HaloExchange does not appear at all — compile-time
gated.

### Auto-injection — exactly two cases, narrowly scoped

**Case A — Capability dependency resolution.** Add to `ComponentBase`:

```cpp
struct ComponentBase {
    /// Capabilities this component READS from the frame each launch.
    /// Default: None. NoseHooverNPT / Berendsen / Deform return Stress.
    /// (LAMMPS analogue: a barostat needing `compute pressure`.)
    virtual Capability reads() const { return Capability::None; }

    // existing: virtual Capability capabilities() const = 0;  // what it WRITES
};
```

In `install_infra()`, the driver computes:

```cpp
provided = ⋃ component.capabilities()
required = ⋃ component.reads()
for each bit in (required \ provided):
    install canonical_provider(bit)   // looked up in a static map, before adding
```

Canonical-provider table (compile-time):

| Required bit | Canonical provider |
|---|---|
| `Stress` | `cpu::md::VirialStress` |
| (others) | none today — adding a bit requires adding a provider entry |

Dedup is at the **driver** level, not the Plan level: before
injecting the canonical provider, the driver checks whether the
pending component set already contains a component **with that
provider's `name()`** (e.g. for Stress, `name() == "virial_stress"`).
If yes, skip injection. This makes the injection rule observable
and testable.

Note on Stress dedup specifically: LJ and bonded force evaluators
also advertise `Capability::Stress` because they accumulate the
pair / bond *virial* into `frame.stress` during force evaluation.
But they do not contribute the *kinetic* term `Σ m·v_α·v_β` that
NPT's pressure feedback needs — that is `VirialStress`'s sole job.
So the dedup is by canonical-provider name, not by capability bit:
LJ being present does not satisfy the NPT pressure dependency, but
an explicit `set_observer(VirialStress{})` does.

**Case B — MPI halo exchange.** Compile-time gated by
`ATV_BUILD_MPI`. When defined AND `set_pair(P)` was called, the
driver installs:

- `cpu::md::mpi::HaloExchangeForward{P.pair_cutoff()}` at `pre_force @ -20`
- `cpu::md::mpi::HaloExchangeReverse{}` at `post_force @ -1000`

When `ATV_BUILD_MPI` is **not** defined, the halo classes don't exist
in the build and no injection runs. When defined and `world_size() == 1`,
the halo components run as no-ops (existing Phase 1 behavior — the
single-rank no-op path is already exercised by every non-MPI test
in the MPI-enabled CI build). No runtime branch in the driver.

This corrects the previous draft's runtime `world_size() > 1`
predicate, which was the wrong question — single- vs multi-rank is
not knowable at `build()` (an `mpiexec -n 4` of a binary built
without MPI is impossible; with MPI it's at user discretion). The
right question is "is MPI in this build?", which is compile-time.

NeighborList is **not** in either auto-injection case. The user must
call `set_neighbor(...)` if they registered a pair / bond / angle /
dihedral force that needs one. The driver throws at `build()` time
when a force evaluator that requires a neighbor list is set without
a corresponding `set_neighbor`. (How "requires NL" is signaled: a
new `bool needs_neighbor_list() const` virtual on `ComponentBase`,
default `false`; pair / bonded evaluators override to `true` when
applicable. Lighter than a full Capability bit and only the driver
reads it.)

### Where everything wires up

- `CPUDriver<P>::build()` calls `install_infra()` before
  `validate_components()`. Default `install_infra()` is no-op.
- Each concrete driver overrides `install_infra()` to apply Cases
  A and B and to throw on missing required NeighborList.
- The `set_*` setters are thin: each stores the typed value(s) in a
  driver-private slot (singular: `std::optional<T>`; append:
  `std::vector<T>`) and calls `add()` / `on().add(c, order)` exactly
  once with the project's documented bindings.

### What does **not** change

- `Driver<>` template, `Plan`, `Stage`, `Backend`, `LaunchContext`,
  `BuildContext` — untouched.
- `ComponentRole`, `Capability`, `has_capability`, `is_none` — untouched.
- `validate_components()` rules — same role + capability requirements.
- Component lifecycle (`build`/`launch`/`destroy`) — untouched.
- `add()` and `on()` remain **public** for user hooks.
- Numerical results — bit-identical to a correctly-written manual
  graph today.
- CUDA backend — out of scope.

## Frame changes

None.

## Files to create / modify

**Create:**
- `include/atomiverse/cpu/pair_evaluator.h` — `PairEvaluator` base at
  backend root (not `cpu/md/`) since cutoff-pair forces are shared
  across MD / minimizer / Evaluate / future MC drivers. Concrete
  subclasses (today: `cpu::md::LJ`) live in their domain directory.
- `include/atomiverse/cpu/minimizer/minimize.h` — extract `Minimize`
  from `lbfgs.h` (it co-habits with the LBFGS algorithm component
  today, leftover from the algorithm-as-driver era).

**Modify:**
- `include/atomiverse/facade/component_role.h` (or wherever
  `ComponentBase` lives) — add `reads()` and `needs_neighbor_list()`
  virtuals, default `Capability::None` and `false`.
- `include/atomiverse/cpu/cpu_driver.h` — `install_infra()` virtual
  hook, called from `build()` before `validate_components()`. Default
  no-op. `add()` / `on()` stay public.
- `include/atomiverse/cpu/md/md.h` — full set of `set_*` setters from
  the table above. `install_infra()` override applying Case A,
  Case B (compile-time MPI gate), and the missing-NL throw.
- `include/atomiverse/cpu/minimizer/minimize.h` (new) — same setters
  appropriate to optimization (no `set_integrator` / `set_thermostat` /
  `set_barostat` — has `set_optimizer` instead).
- `include/atomiverse/cpu/evaluate.h` — same setters appropriate
  to one-pass evaluation (no `set_integrator` / `set_thermostat` /
  `set_barostat` / `set_optimizer` / `set_thermo`).
- `include/atomiverse/cpu/md/lj.h`, `bonded/*.h` — override
  `needs_neighbor_list()` to return `true` for pair-cutoff
  evaluators that need a neighbor walker.
- `include/atomiverse/cpu/md/nose_hoover_npt.h`,
  `barostat.h`, `deform.h` — override `reads()` to return
  `Capability::Stress`.

**Tests / examples / docs:**
- Migrate every test under `tests/cpu/md/`, `tests/mpi/`, and the
  AIMD tests under `tests/cpu/scf/` to the categorical setters.
- Migrate `examples/` and `docs/cpu/` tutorials.
- Update CLAUDE.md "Example" block.

## Dependency graph

- Independent of `parallel-md-engine.md` (in progress) — only touches
  user surface.
- `cuda-driver-taxonomy.md` (planned) adopts the same shape from day
  one. CUDA halo (when implemented) goes under the same `ATV_BUILD_MPI`
  + `set_pair` rule.
- `python-bindings.md` — `mp.MD` exposes the same setters as Python
  methods: `sim.set_pair(mp.LJ(...))`, `sim.set_integrator(...)`. No
  translation layer.

## Tasks

- [x] T1 — `Capability ComponentBase::reads()` + `bool ComponentBase::needs_neighbor_list()`
      virtuals (defaults `None` / `false`). Existing components compile unchanged.
      *Done 2026-05-06; defaults verified by `tests/cpu/test_install_infra_hook.cpp`.*
- [x] T2 — `CPUDriver<P>::install_infra()` virtual hook called from
      `build()` before `validate_components()`. `add()` / `on()` stay public.
      *Done 2026-05-06.*
- [x] T3 — Categorical setters on `cpu::md::MD`: `set_pair`, `set_bond`,
      `set_angle`, `set_dihedral`, `set_neighbor`, `set_integrator`,
      `set_thermostat`, `set_barostat`, `set_observer`, `set_record`,
      `set_thermo`. Singular setters store `std::optional<T>` and
      throw on second-call-with-different-value; append setters keep
      insertion order. *Done 2026-05-06.*
- [x] T4 — `MD::install_infra()`: (Case A) Stress dep → VirialStress
      with driver-level dedup by canonical-provider name (`virial_stress`);
      (Case B) `#ifdef ATV_BUILD_MPI` halo exchange auto-bind when
      `set_pair` is set; throw on missing `set_neighbor` when a force
      evaluator with `needs_neighbor_list()==true` is registered.
      *Done 2026-05-06.*
- [x] T5 — Override `reads()` on NoseHooverNPT / Barostat
      to return `Capability::Stress`. Override `needs_neighbor_list()`
      on LJ via the new `cpu::PairEvaluator` base.
      Deform deferred (does not actually read stress today; spec wording
      was overly broad). *Done 2026-05-06.*
- [x] T6 — Mirror setters on `cpu::minimizer::Minimize`. Move
      `Minimize` into `minimizer/minimize.h`. Add `set_optimizer`.
      *Done 2026-05-06.*
- [x] T7 — Mirror setters on `cpu::Evaluate`. No `set_integrator` /
      `set_thermostat` / `set_barostat` / `set_optimizer` / `set_thermo`.
      *Done 2026-05-06.*
- [ ] T8 — **Deferred.** Migrate every test, example, tutorial to the
      categorical surface. Manual `set_neighbor` everywhere LAMMPS
      would say `neighbor 0.3 bin`. No manual `add(VirialStress{})`
      in NPT tests. The new API is additive, so existing tests still
      compile and pass; migration improves readability but is
      mechanical and can land in follow-up PRs.
- [x] T9 — CLAUDE.md "Example" block updated to advertise the new
      API and the auto-injection rules. `docs/cpu/*` tutorial
      migration deferred alongside T8. *Done 2026-05-06.*
- [x] T10 — **Plan-equivalence test.** `tests/cpu/md/test_md_plan_equivalence.cpp`
      proves that NVE-LJ and NPT-LJ simulations built via the
      categorical setter API produce a byte-identical compiled
      component graph compared to the raw `add()` / `on().add()`
      path with explicit infrastructure. *Done 2026-05-06.*
- [x] T11 — Singular-setter discipline. Verified for `set_integrator`
      (different type → throws), `set_pair` / `set_neighbor` /
      `set_optimizer` (same name idempotent). Different-class re-call
      for `set_pair` deferred until a second pair evaluator (e.g.
      Coulomb-real-space) exists in the codebase.
      *Covered by `tests/cpu/md/test_md_set_api.cpp` and
      `tests/cpu/test_minimize_evaluate_set_api.cpp`. Done 2026-05-06.*
- [x] T12 — Capability-dep resolution + driver-level dedup verified
      via `IntrospectMD` snapshot of pending-component name counts
      at validate-time. *Done 2026-05-06.*
- [ ] T13 — **Deferred.** MPI compile-time gating CI matrix
      (`ATV_BUILD_MPI=ON`/`OFF`). The compile-time `#ifdef` lives in
      `MD::install_infra()`; the runtime no-op single-rank path is
      already exercised by every non-MPI test in the MPI-enabled
      build. A dedicated CI job toggling `ATV_BUILD_MPI=OFF` is the
      remaining piece — adds ~5 minutes to CI matrix time, lands in
      a follow-up PR.
- [x] T14 — Missing-NL throw verified across `MD`, `Minimize`,
      `Evaluate`. Message names the offending evaluator and the
      `set_neighbor` instruction. *Done 2026-05-06.*
- [x] T15 — Escape-hatch coexistence: a user-defined hook registered
      via `sim.on(end_of_step).add(CustomEndOfStepHook{}, +50)`
      appears in the compiled Plan alongside the `set_*` registrations.
      *Done 2026-05-06.*

## Testing strategy

**Plan equivalence (load-bearing).** T10 proves the new API is sugar:
same components, same phases, same orders as a hand-written graph.

**Numerical regression.** `tests/mpi/test_pe_msu_pipeline.cpp` and
`tests/mpi/test_npt_multirank.cpp` migrated to the new surface
produce identical thermodynamic traces — drift on total energy ≤ 1e-12.

**Singular-setter discipline.** T11.

**Capability-dep resolution + dedup.** T12 — both that injection
fires and that explicit registration coalesces.

**MPI gating semantics.** T13 — compile-time, not runtime. CI matrix
must include both `-DATV_BUILD_MPI=ON` and `=OFF` configurations and
both must pass.

**Missing-NL throw.** T14 — protects against the user removing
`set_neighbor` and silently ending up with O(N²) or no-force results.

**Escape hatch.** T15 — `add()` / `on()` still work. Custom hooks
remain a first-class citizen.

**Tutorial / example.** Tutorials in `docs/cpu/` compile and run
end-to-end with the new API.

## Out of scope

- CUDA mirror — separate spec, adopts same shape from day one.
- Multi-rank `cpu::Evaluate` — when there's a real workflow, it
  reuses Case B by adding `set_pair` halo binding to `Evaluate::install_infra()`.
- Variable-cell minimization (`Minimize` does not need stress today).
- Group selection (LAMMPS `fix N group_id …`). All current
  Atomiverse components apply to all atoms; introducing groups is
  a separate spec.
- `pair_modify` / `kspace_style` / `comm_modify` — full LAMMPS
  modifier surface. Each future modifier is its own incremental
  spec when there is a real use case.
- Auto-deriving `set_neighbor` skin / cadence from `set_pair`
  cutoff. The user supplies them, like LAMMPS. Sensible defaults
  (`skin=0.3`, `every=10`, `bin`) live on `NeighborList<>`'s
  constructor defaults.
- Removing or hiding `add()` / `on()`. Kept public for user hooks.

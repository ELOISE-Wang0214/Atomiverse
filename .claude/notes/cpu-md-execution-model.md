# CPU MD Execution Model & Architecture Normalization

Status: **design doc** (write-first; the refactor it specifies is staged separately).
Reference: LAMMPS *Developer flow* — https://docs.lammps.org/Developer_flow.html

This is the developer-flow reference for the CPU MD engine: the explicit
execution skeleton, what every phase is responsible for (**even when empty**),
how components bind into it, and the naming/structure conventions the code
should converge on. The goal is that a reader can see the entire timestep at a
glance — like reading `Verlet::run()` in LAMMPS — without chasing call sites.

---

## 1. Top-level flow (driver lifetime)

Every driver (`MD`, `Minimize`, `Evaluate`, …) runs the same three-stage shape,
inherited from `CPUDriver` (`cpu/cpu_driver.h`):

```
run(frame, n_steps):
    setup_pass(frame)              # once, on the first run() of a trajectory
        backend.init()             #   allocate backend, timer
        backend.load(plan)         #   the compiled Plan (phase-ordered slots)
        backend.set_frame(frame)   #   bind frame; partition_atoms (multi-rank)
        backend.build_components() #   each component's build(): one-time alloc
    for s in [0, n_steps):         # the run loop
        step_once(s)               #   backend.step(s) — one full timestep
    backend.drain_analysis()       # flush async observers
    # teardown_pass() only on driver destruction / frame switch
```

LAMMPS map: `setup_pass`↔`init()`+`setup()`, the loop↔`Verlet::run()`,
`drain_analysis`/`teardown_pass`↔`post_run()`. **Continuable:** consecutive
`run()` calls extend one trajectory; `build()` happens once so stateful
integrators (NHC chains, RNG) keep their state.

---

## 2. The timestep — phases in order

One `backend.step(idx)` walks the compiled **Plan**: a list of *slots*, one per
`(phase, order)`, sorted by phase then ascending order. Each slot holds the
components bound there; each is `launch(ctx)`-ed in order. A component may also
be frequency-gated (`run every N steps`). This is AV's analogue of LAMMPS
iterating its fix list at each `*_integrate`/`*_force` callback.

`md::Phase` (`facade/stage.h`) — the 6 phases AV defines today:

| # | phase | responsibility (one line) | LAMMPS equivalent |
|---|---|---|---|
| 0 | `initial_integrate` | first half-kick `v += ½dt·f/m`, drift `x += dt·v`; barostat box push | `initial_integrate` |
| 1 | `pre_force` | everything that must precede force eval: migrate, build ghosts, build neighbor list, zero accumulators | `post_integrate`+`pre_exchange`+`pre_neighbor`+`pre_force` (merged) |
| 2 | `force` | pair + bonded force/energy/virial accumulation | `pair/bond/...->compute()` |
| 3 | `post_force` | reverse ghost-force comm; virial finalize; external/constraint forces | `pre_reverse`+reverse comm+`post_force` (merged) |
| 4 | `final_integrate` | second half-kick `v += ½dt·f/m`; barostat velocity scale | `final_integrate` |
| 5 | `end_of_step` | global reductions (energy, KE, virial), thermostat/barostat bookkeeping, observers, recorders, load balance | `end_of_step`+output |

### 2.1 What actually binds where (current)

Ordering uses small integers; **negative = infrastructure that must run first**,
0 = the user-facing op, **positive = consumers that run after**.

```
initial_integrate:
   0     integrator            velocity_verlet | nve_limit | nose_hoover_npt
  +200   deform                (box deformation, stage-5 tension)

pre_force:
  -1000  zero_forces           f[*] = 0 over owned rows
   -999  zero_energy           scalar energies = 0
   -998  zero_stress           virial = 0
   -30   atom_migration        re-home atoms that crossed a rank boundary
   -20   halo_forward          send owned → ghost positions (+ zero ghost forces)
   -10   neighbor_list         (re)build the cell list / pair list

force:
   0     pair                  lj (pair_lj_cut) | morse
   0     bonded                harmonic_bond | harmonic_angle | multi_harmonic_dihedral

post_force:
   -50   virial_stress         finalize per-atom → global virial
  -1000  halo_reverse          sum ghost forces back to owners

final_integrate:
   0     integrator            (same component, second half)

end_of_step:
  -1000  zero_kinetic          KE = 0
   -100  kinetic_energy        accumulate ½Σm v² (owned)
   -50   reduce_energy         one batched MPI_Allreduce of {E, per-term, KE, virial}
   -45   barostat              (berendsen) box rescale
   -40   nose_hoover_npt       NHC + box bookkeeping
   +10   snapshot_publisher    publish frame snapshot to async workers
   +20   async observers       rdf, thermo, recorders (default slot)
  +1000  rcb_balancer          recursive-coordinate-bisection re-partition
```

### 2.2 The empty-phase rule (the LAMMPS lesson)

LAMMPS lists *every* callback in `Verlet::run()` even when no fix uses it, so the
skeleton is legible. AV must do the same: **all 6 phases always exist and are
always walked**, even with zero components bound (the slot loop simply skips
empty slots). A reader sees the full timestep from the `Phase` enum + the table
above; "nothing happens in `post_force` for a bare NVE run" is then an explicit,
documented fact, not an absence to be inferred.

---

## 3. Component contract

A component (`cpu::ComponentBase`) is the unit of work, ~LAMMPS `Fix`/`Pair`:

```cpp
std::string name() const;                       // timing/identity key, snake_case
ComponentRole role() const;                     // Integrator | Evaluator | None ...
Capabilities capabilities() const;              // Force | Energy | Stress ...
std::vector<DefaultBinding> default_bindings();  // {(phase, order), ...}
void build(const BuildContext&);                // one-time alloc (NOT in the loop)
void launch(const LaunchContext&);              // per-step work (NO allocation)
void destroy(const DestroyContext&);            // teardown
```

Invariants (already enforced; keep): **no allocation in `launch()`** (warm-up in
`build()`; growth only on rare capacity events), forces accumulate into the
shared `frame.atoms.f*`, owned rows are `[0, n_local)` and ghosts
`[n_local, n_local+n_ghost)` (`owned_range.h`, `ResolveView`).

---

## 4. Gaps vs LAMMPS (what the normalization fixes)

1. **Merged phases hide structure.** AV folds `post_integrate`, `pre_exchange`,
   `pre_neighbor`, `post_neighbor` into `pre_force`, and `pre_reverse`+reverse
   into `post_force`, ordering them with magic negative integers (`-30`, `-20`,
   `-10`). A reader cannot tell "migrate before neighbor build" from the enum;
   it is encoded only in the integer constants scattered across `default_bindings`.
2. **Magic order integers.** `-1000`, `-999`, `-50`, `+1000` are unnamed. Their
   relative order is load-bearing (zero before accumulate; reduce before
   publish) but undocumented at the binding site.
3. **Naming drift.** Mixed lengths/styles: `HaloExchangeForward`,
   `DihedralMultiHarmonic`, `DomainDecomposition`, `NoseHooverNPT` vs the clean
   `NeighborList`, `VelocityVerlet`. Component `name()` strings
   (`multi_harmonic_dihedral`, `vib_hessian_accumulator`) are long.

---

## 5. Normalized design

### 5.1 Explicit phase set (self-documenting skeleton)

Adopt the LAMMPS-granular phases as **named ordering anchors** so binding sites
read in English, not integers. Two options (pick in review):

- **(A) Add real phases** to `md::Phase`: `initial_integrate, post_integrate,
  pre_exchange, pre_neighbor, neighbor, pre_force, force, post_force,
  final_integrate, end_of_step`. Components bind to the *named* phase; intra-
  phase order shrinks to {-1,0,+1}. Most new phases are usually empty — that is
  the point (legible skeleton). Cost: touch every `default_bindings`.
- **(B) Named order bands** (cheaper): keep 6 phases, replace magic integers
  with a documented enum `enum Order { Migrate=-30, Ghost=-20, Neigh=-10,
  Zero=-1000, ... }` in one header. Binding sites read `{pre_force, Order::Neigh}`.

Recommendation: **(A)** — it is the LAMMPS model the request asks for; the empty
phases are the explicit "does nothing here" the user wants.

### 5.2 Naming convention

Rule: **one or two words, lower-domain-jargon, no redundant suffixes.** Class
names `UpperCamel`, component `name()` strings `snake_case` matching the class.

| current | proposed | why |
|---|---|---|
| `HaloExchangeForward` / `…Reverse` | `HaloSend` / `HaloRecv` | 2 words, directional, no "Exchange" filler |
| `DihedralMultiHarmonic` | `Dihedral` (style = ctor arg) or `MultiDihedral` | the "Harmonic"/"Multi" is a coeff detail |
| `BondHarmonic` / `AngleHarmonic` | `Bond` / `Angle` | only one style ships; style is data |
| `NoseHooverNPT` | `BaroNHC` (or keep — it is a proper noun) | optional |
| `DomainDecomposition` | `Domain` | 1 word; it *is* the domain |
| `AtomMigration` | `Migrate` | verb, 1 word |
| `NeighborList` | `Neighbors` | keep / shorten |
| `VirialStress` | `Virial` | 1 word |
| `ReduceEnergy` | `Reduce` | scope is energies |
| `SnapshotPublisher` | `Publish` | verb |
| `RCBBalancer` | `Rebalance` | verb |
| `polyethylene_tension` (bench) | `pe_md` | bench, not core |

(Renames are mechanical but wide; do them per-component with the OpenAPI/test
grep in one commit each, behavior-preserving.)

### 5.3 Module layout (already mostly right)

Keep the pillar split (`md/` integrators+infra, `potentials/` force kernels,
`comm/` halo+migration, shared infra at `cpu/` root). The execution-model types
(`Phase`, `Plan`, `ComponentBase`, `Backend`) are the public "how a step runs"
surface and should be the first thing `architecture.md` links to.

---

## 6. Refactor sequence (when approved)

1. Land §5.1(A): expand `md::Phase`, migrate `default_bindings` to named phases,
   delete magic integers. Gate: full MPI + non-MPI suite green, timing
   breakdown unchanged.
2. Land §5.2 renames one component at a time (class + `name()` + call sites +
   tests), each its own behavior-preserving commit.
3. Update `architecture.md` to point here as the execution-model reference.

Each step is test-gated and behavior-preserving; no physics change.

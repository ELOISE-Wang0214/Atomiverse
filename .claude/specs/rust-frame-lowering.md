---
title: Rust-First Frame Lowering and Typed Requirement Refactor
status: draft
created: 2026-04-10
depends_on:
  - python-bindings.md
# frame-block-refactor: satisfied (spec retired 2026-05-19; cpu::Frame block
# model + atoms/bonds/etc. naming already ships).
---

## Summary

Atomiverse will move to a Rust-first control plane:

```text
py.Frame -> rs.Frame -> cpu::Frame / cuda::Frame / torch::Frame
```

The canonical, cross-language data model is `rs::Frame` (`molrs::Frame` or a compatible Rust wrapper around it). Backend-native frames are execution artifacts, not public interchange types.

The public execution API becomes:

```cpp
Spec -> compile() -> Plan -> Driver(plan)
```

Lowering from `rs::Frame` into backend-native frames is driven by a **typed requirement system**, not by hardcoded field lists or string slot keys.

This spec replaces the earlier string-based semantic slot direction. Public requirement APIs must be type-driven. Any string names or DLPack descriptors are internal lowering details only.

## Goals

1. Make `rs::Frame` the sole canonical host-side frame for Python and Rust.
2. Treat `cpu::Frame`, `cuda::Frame`, and `torch::Frame` as peer backend artifacts.
3. Split compile-time composition from runtime execution:
   - `Spec` collects components.
   - `Plan` is the compiled schedule plus lowering requirements.
   - `Driver` only executes a compiled plan.
4. Eliminate public string-based requirement APIs.
5. Let components declare dependencies in a typed, open-ended way.
6. Keep Rust lowering generic over arbitrary `rs::Frame` fields and types.
7. Align built-in semantic naming with MolPy naming conventions.
8. Preserve a transitional compatibility path for existing `MD`, `SCF`, `LBFGS`, `SD`, `CG` wrappers.

## Non-goals

1. Do not make DLPack the public semantic API.
2. Do not require a single central enum or central string registry for every possible semantic slot.
3. Do not force all backends to share the same physical storage layout.
4. Do not make `cpu::Frame` the mandatory ingress path for CUDA or Torch.
5. Do not expose backend-internal tensors or device buffers as the canonical cross-language frame model.

## Problem

Atomiverse currently mixes three separate concerns:

1. Public composition of simulation pipelines.
2. Runtime execution on CPU or CUDA.
3. Field-level assumptions used for frame materialization.

The old design has several issues:

1. `Driver` historically owned both the pre-runtime builder state and runtime execution state.
2. Python conversion logic hardcodes concrete fields such as `atoms.x`, `atoms.y`, `atoms.z`, `atoms.type`.
3. Existing requirement experiments used string slot identity such as `"atoms.position"`, which is still hardcoded and not type-safe.
4. `cpu::Frame` still occupies a privileged architectural position that the user explicitly wants removed.

The new design must fix all four.

## Architectural direction

### Canonical layer vs backend layer

The architecture is intentionally asymmetric:

```text
Python / Rust API
    |
    v
rs::Frame
    |
    +--> lower(plan.requirements(), backend=cpu)   -> cpu::Frame
    +--> lower(plan.requirements(), backend=cuda)  -> cuda::Frame
    +--> lower(plan.requirements(), backend=torch) -> torch::Frame
```

`rs::Frame` is the truth layer.

`cpu::Frame`, `cuda::Frame`, and `torch::Frame` are backend-native execution caches:

1. They are allowed to use fixed layouts for performance.
2. They are allowed to omit fields not required by the compiled plan.
3. They are allowed to have backend-specific workspaces and runtime state.

### Public execution API

New canonical API:

```cpp
cpu::md::Spec spec;
spec.add(VelocityVerlet{0.001});
spec.on(ATV::md::Phase::force).add(NeighborList{6.8}, -10);
spec.on(ATV::md::Phase::force).add(LJ{0.0104, 3.4, 6.8});

auto plan = spec.compile();

// Rust side, not C++:
// cpu::Frame frame = rs::materialize<cpu::Frame>(rs_frame, plan.requirements());

cpu::md::Driver driver(std::move(plan));
driver.run(frame, 1000);
```

Legacy wrappers such as `cpu::md::MD` remain temporarily available as compatibility facades over `Spec + Driver`.

## Requirement model

### Public requirement API must be typed

Public requirement identity is:

```cpp
Slot<Domain, Tag, ValueSpec>
```

The first two axes capture semantics:

1. `Domain`: which namespace or entity the field belongs to
2. `Tag`: which semantic field is being requested

The third axis captures representation requirements:

3. `ValueSpec`: dtype + shape

Example:

```cpp
Slot<Atoms, X, Tensor<F64, Scalar>>
Slot<Atoms, Velocity, Tensor<F64, Fixed<3>>>
Slot<Atoms, Force, Tensor<F64, Fixed<3>>>
Slot<System, Box, Tensor<F64, Fixed<3, 3>>>
```

Component code declares dependencies through the type system:

```cpp
RequirementSet req;
req.read<Slot<Atoms, X, Tensor<F64, Scalar>>>();
req.read<Slot<Atoms, Y, Tensor<F64, Scalar>>>();
req.read<Slot<Atoms, Z, Tensor<F64, Scalar>>>();
req.read<Slot<Atoms, Mass, Tensor<F64, Scalar>>>();
req.accumulate<Slot<Atoms, Force, Tensor<F64, Fixed<3>>>>();
```

No public API should accept:

1. `const char*`
2. `std::string`
3. path strings such as `"atoms.position"`
4. a central string table as the canonical requirement identity

### Built-in semantic names follow MolPy naming conventions

Built-in domains and tags should track MolPy naming convention concepts one-to-one.

The recommended built-in semantic types are:

```cpp
namespace ATV::molpy {

struct Atoms {};
struct Bonds {};
struct Angles {};
struct Dihedrals {};
struct Impropers {};
struct System {};
struct Neighbors {};

struct X {};
struct Y {};
struct Z {};
struct Velocity {};
struct Force {};
struct Mass {};
struct Type {};
struct Number {};
struct Charge {};

struct Atomi {};
struct Atomj {};
struct Atomk {};
struct Atoml {};

struct Box {};
struct Energy {};
struct Kinetic {};
struct PairList {};

}
```

Notes:

1. The C++ types are PascalCase because they are types.
2. Their semantic meaning maps directly onto MolPy keys such as `atoms`, `x`, `number`, `atomi`, `box`.
3. Built-in types are convenience only. New components are free to define their own domains and tags locally.

### ValueSpec must include dtype and shape

`ValueSpec` should be split into dtype and shape at the type level:

```cpp
namespace req::dtype {
struct F32 {};
struct F64 {};
struct I32 {};
struct I64 {};
struct Bool {};
}

namespace req::shape {
struct Scalar {};
template<std::size_t... Dims> struct Fixed {};
template<std::size_t Rank> struct Dynamic {};
}

template<class DType, class Shape>
struct Tensor {};

template<class T>
struct Opaque {};
```

Examples:

```cpp
using AtomX = Slot<molpy::Atoms, molpy::X,
                   Tensor<req::dtype::F64, req::shape::Scalar>>;

using AtomVelocity = Slot<molpy::Atoms, molpy::Velocity,
                          Tensor<req::dtype::F64, req::shape::Fixed<3>>>;

using AtomForce = Slot<molpy::Atoms, molpy::Force,
                       Tensor<req::dtype::F64, req::shape::Fixed<3>>>;

using SystemBox = Slot<molpy::System, molpy::Box,
                       Tensor<req::dtype::F64, req::shape::Fixed<3, 3>>>;

using NeighborPairs = Slot<molpy::Neighbors, molpy::PairList,
                           Opaque<PairListHandle>>;
```

This is preferred over using `std::tuple` as the public `ValueSpec`.

Reason:

1. `tuple` is a product type, not a shape descriptor.
2. Shape types must support fixed rank, dynamic rank, and opaque runtime artifacts.
3. Lowering needs explicit access to shape metadata, not just an arbitrary nested type.

### Type erasure is internal only

Runtime aggregation and lowering can erase the type identity into internal descriptors.

Internal shape:

```cpp
struct TypeToken;

struct SlotDescriptor {
    TypeToken domain;
    TypeToken tag;
    TypeToken dtype;
    TypeToken shape;
};

struct RequirementDescriptor {
    SlotDescriptor slot;
    bool read;
    bool write;
    bool accumulate;
};

class RequirementSet;          // phase-local typed declarations
class CompiledRequirementSet;  // phase-bound internal descriptor table
```

The component author should only see the typed API.

## Why DLPack is not the public requirement API

DLPack is useful, but at the wrong layer.

DLPack describes runtime tensors:

1. dtype
2. ndim
3. shape
4. strides
5. device

It does not encode semantic identity such as:

1. `atoms.x`
2. `atoms.v`
3. `bonds.atomi`
4. `system.box`

Therefore:

1. Public semantic requirement API uses `Slot<Domain, Tag, ValueSpec>`.
2. Internal lowering may translate `ValueSpec` into a DLPack-compatible runtime descriptor.
3. Torch / CUDA interop may use DLPack at the handoff boundary.

In short: DLPack is a transport/runtime descriptor, not a semantic slot system.

## Lowering model

### Plan requirements drive materialization

Each component contributes phase-local typed requirements.

`compile()` aggregates them into:

1. execution schedule
2. compiled lowering requirements

The Rust lowering layer consumes:

```cpp
plan.requirements()
```

and materializes only the required subset of `rs::Frame` into the chosen backend frame.

### Lowering must stay open-ended

The Rust lowering system must not require a central schema table for every possible slot.

Instead, the model is:

1. built-in Atomiverse semantics ship with built-in lowering traits
2. custom slots can provide their own lowering traits or adapters
3. the lowering engine consumes internal descriptors, not public strings

Conceptually:

```cpp
template<class SlotT>
struct LoweringTraits;
```

Possible responsibilities of `LoweringTraits<SlotT>`:

1. map the slot to a MolPy/MolRS namespace and field when such a mapping exists
2. describe the target backend storage class
3. describe copy, pack, or scatter behavior
4. optionally declare whether the field is optional or synthesized

Built-in slots such as `Slot<Atoms, X, Tensor<F64, Scalar>>` can lower to:

```text
rs_frame["atoms"]["x"]
```

but that mapping is an internal trait, not the public requirement identity.

### Backend-specific lowering

Lowering is backend-specific and plan-driven:

1. CPU lowering builds `cpu::Frame`
2. CUDA lowering builds `cuda::Frame`
3. Torch lowering builds `torch::Frame`

The same slot may lower differently by backend:

1. `Velocity` may be three SoA arrays in CPU and CUDA
2. `Velocity` may be one tensor view in Torch
3. `PairList` may be skipped during initial materialization and produced lazily by `NeighborList`

## Driver and facade refactor

### Spec / Plan / Driver split

Public responsibilities:

#### `Spec`

Owns:

1. component instances
2. default bindings
3. explicit phase bindings

Does not own:

1. runtime buffers
2. frames
3. backend state

#### `Plan`

Owns:

1. sorted stage schedule
2. compiled lowering requirements
3. component ownership transferred from `Spec`

Does not own:

1. active runtime context
2. backend frame memory
3. execution statistics

#### `Driver`

Owns:

1. compiled `Plan`
2. runtime backend object
3. optional timing state

Does not own:

1. composition-time builder state
2. unresolved bindings
3. materialization policy

### Legacy compatibility layer

Existing domain wrappers remain as transitional facades:

1. `cpu::md::MD`
2. `ATV::md::MD`
3. `ATV::scf::SCF`
4. `cpu::minimizer::LBFGS`
5. `ATV::minimizer::SD`
6. `ATV::minimizer::CG`

These wrappers should compose:

```text
Spec + lazy compile() + Driver
```

and expose the old:

```cpp
add()
on()
build()
run()
```

API only as a temporary compatibility shim.

## Rust and Python integration

### Python path

Target Python path:

```text
molpy.Frame -> Rust -> rs::Frame -> lower(plan.requirements()) -> backend frame
```

Python should stop manually building `cpu::Frame` field-by-field.

The current hardcoded Python converter becomes transitional and will be removed after the Rust materializer is stable.

### Rust responsibilities

Rust owns:

1. canonical frame validation
2. lowering based on compiled requirements
3. backend selection
4. caching and synchronization policy
5. future multi-language public API

Hot numerical kernels remain in C++ / CUDA / Torch backends.

### Synchronization policy

Backends do not implicitly round-trip everything on every run.

Synchronization is explicit and mask-driven:

1. positions
2. velocities
3. forces
4. energy
5. kinetic
6. topology if needed

This matters most for CUDA and Torch repeated-run workflows.

## Transitional rules

### Allowed transitional pieces

The following are allowed temporarily:

1. existing legacy `MD/SCF/LBFGS/...` wrappers
2. string-based or field-hardcoded Python adapter code that is clearly marked transitional
3. temporary dual-path materialization during migration

### Forbidden future-facing interfaces

The following must not be extended:

1. public string slot APIs
2. central string slot registries
3. new component APIs that declare requirements via strings
4. new Rust materialization entry points that hardcode concrete field lists in the public API

## Acceptance criteria

The refactor is complete when all of the following are true:

1. Public requirement APIs are typed.
2. `Plan::requirements()` returns compiled typed descriptors, not string-keyed requirements.
3. Built-in component requirements use typed slots.
4. New built-in semantics align with MolPy naming conventions.
5. Rust lowering consumes compiled descriptors and no longer depends on a public string slot API.
6. `Spec -> compile() -> Plan -> Driver(plan)` is the canonical public execution path.
7. Legacy facades remain functional during migration.
8. `rs::Frame -> cpu::Frame` and `rs::Frame -> cuda::Frame` are independent lowering paths.

## Roadmap

### Phase 0: Reset the requirement direction

Goal:

1. Remove the string-slot direction from the design.
2. Freeze any public API that still treats requirement identity as strings.

Deliverables:

1. This spec.
2. `requirements.h` rewritten around typed `Slot<Domain, Tag, ValueSpec>`.

### Phase 1: Stabilize typed requirement substrate

Goal:

1. Implement type-level `Domain`, `Tag`, `DType`, and `Shape`.
2. Implement internal type erasure for compiled requirements.

Deliverables:

1. `RequirementSet`
2. `CompiledRequirementSet`
3. `SlotDescriptor`
4. tests for merging and lookup of typed requirements

### Phase 2: Land canonical built-in semantics

Goal:

1. Add built-in semantic types aligned with MolPy naming conventions.
2. Keep the system open for custom domains and tags.

Deliverables:

1. `facade/molpy_semantics.h`
2. built-in dtype/shape vocabulary
3. conventions doc for custom slot authors

### Phase 3: Complete the Spec / Plan / Driver split

Goal:

1. Make `Driver` a pure execution object.
2. Move compile-time responsibilities into `Spec`.
3. Keep legacy wrappers working as transitional facades.

Deliverables:

1. updated generic facade
2. `cpu::md::Spec`, `cpu::md::Driver`
3. `ATV::md::Spec`, `ATV::md::Driver`
4. `ATV::scf::Spec`, `ATV::scf::Driver`
5. `cpu::minimizer::Spec`, `cpu::minimizer::Driver`
6. `ATV::minimizer::Spec`, `ATV::minimizer::Driver`

### Phase 4: Convert built-in components to typed requirements

Goal:

1. Stop returning ad hoc empty or string-driven requirements.
2. Teach core components to describe the fields they need.

Initial component set:

1. `ZeroForces`
2. `ZeroEnergy`
3. `ZeroKinetic`
4. `KineticEnergy`
5. `VelocityVerlet`
6. `NeighborList`
7. `LJ`

Follow-up:

1. minimizer components
2. SCF components
3. observer components

### Phase 5: Add lowering traits and Rust materializer

Goal:

1. Map built-in typed slots onto `rs::Frame` field access.
2. Materialize backend-native frames directly from compiled requirements.

Deliverables:

1. built-in lowering traits
2. `rs::materialize<cpu::Frame>(rs_frame, plan.requirements())`
3. `rs::materialize<cuda::Frame>(rs_frame, plan.requirements())`
4. sync-back policy

### Phase 6: Python migration

Goal:

1. Move Python off hardcoded `molpy.Frame -> CpuFrame` conversion.
2. Route Python through Rust-first materialization.

Deliverables:

1. Python binding path updated to `molpy.Frame -> Rust -> rs::Frame`
2. existing converter marked deprecated
3. integration tests for MD and SCF from `molpy.Frame`

### Phase 7: Backend-native repeated-run workflows

Goal:

1. Reuse backend frames across repeated runs.
2. Avoid rebuilding or retransferring host state unnecessarily.

Deliverables:

1. explicit materialize API
2. explicit sync API
3. repeated-run benchmarks

### Phase 8: Remove transitional public paths

Goal:

1. delete dead compatibility layers once Rust-first path is stable.

Deliverables:

1. remove public string slot remnants
2. retire old Python converter path
3. reduce or remove legacy driver facade where no longer needed

## Test plan

### Requirement substrate tests

1. typed requirement merging
2. compiled phase-bound aggregation
3. no string-based public APIs in new component declarations

### Facade tests

1. `Spec -> compile() -> Driver(plan)` on CPU
2. legacy `MD/SCF/LBFGS/...` wrappers still function
3. `plan.requirements()` is stable and queryable

### Lowering tests

1. `rs::Frame -> cpu::Frame` based on compiled requirements
2. `rs::Frame -> cuda::Frame` based on compiled requirements
3. no forced `cpu::Frame -> cuda::Frame` ingress in the canonical Rust path

### Integration tests

1. `molpy.Frame -> Rust -> MD.run()`
2. `molpy.Frame -> Rust -> SCF.run()`
3. repeated-run backend reuse

## Files expected to change

Core facade:

1. `include/atomiverse/facade/requirements.h`
2. `include/atomiverse/facade/molpy_semantics.h`
3. `include/atomiverse/facade/plan.h`
4. `include/atomiverse/facade/driver.h`
5. `include/atomiverse/facade/driver_impl.h`

Domain wrappers:

1. `include/atomiverse/cpu/md/md.h`
2. `include/atomiverse/cuda/md/md.cuh`
3. `include/atomiverse/cuda/scf/scf.cuh`
4. `include/atomiverse/cpu/minimizer/lbfgs.h`
5. `include/atomiverse/cuda/minimizer/sd.cuh`
6. `include/atomiverse/cuda/minimizer/cg.cuh`

Built-in components:

1. MD init / integrator / neighbor / LJ headers
2. minimizer components
3. SCF components

Rust / Python path:

1. Rust materializer crate
2. Python binding entry points
3. old `_convert.py` path


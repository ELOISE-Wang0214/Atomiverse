---
title: Python Binding Layer (nanobind + DLPack)
status: in-progress
created: 2026-04-04
---

## Summary

Python bindings for Atomiverse via nanobind. `mp.Frame` (molpy) is the sole user-facing Frame type. Atomiverse is a compute backend that borrows data from `mp.Frame`, runs on CPU or CUDA, and writes results back. DLPack is the zero-copy exchange protocol: numpy for CPU, `torch.from_dlpack()` for GPU tensors. No Atomiverse-specific FFI to molrs::Frame — numpy is the lingua franca.

## mp.Frame data model

molpy's Frame is a **dict of Blocks**, not a flat struct:

```python
frame = mp.Frame()
frame["atoms"] = mp.Block({
    "x": np.array([0.0, 1.0, 2.0]),   # (N,) SoA — same layout as cpu::Frame
    "y": np.array([0.0, 0.0, 0.0]),
    "z": np.array([0.0, 0.0, 0.0]),
    "mass": np.array([12.0, 1.0, 1.0]),
    "type": np.array([0, 1, 1]),
})
frame["bonds"] = mp.Block({
    "atomi": np.array([0, 0]),
    "atomj": np.array([1, 2]),
    "type": np.array([0, 0]),
})
frame.box = mp.Box(np.eye(3) * 30.0)

n_atoms = frame["atoms"].nrows        # derived, never a constructor arg
```

**Key insight:** molpy positions are already SoA (`x`, `y`, `z` as separate 1D arrays). No AoS↔SoA conversion needed — the transfer to `cpu::Frame` is direct memcpy per field.

### Naming: molcrafts conventions

Atomiverse is part of molcrafts. Field names follow the molpy naming conventions (`molpy/docs/getting-started/naming-conventions.md`). No mapping table — same names on both sides.

**Canonical field names** (used in molpy blocks AND Atomiverse Frames):

| block | fields | dtype |
|-------|--------|-------|
| `atoms` | `x`, `y`, `z`, `vx`, `vy`, `vz`, `fx`, `fy`, `fz`, `mass`, `type` | Real / Int |
| `bonds` | `atomi`, `atomj`, `type` | Int |
| `angles` | `atomi`, `atomj`, `atomk`, `type` | Int |
| `dihedrals` | `atomi`, `atomj`, `atomk`, `atoml`, `type` | Int |
| `box` | 3×3 matrix | Real |
| (scalars) | `energy`, `kinetic`, `stress[6]` | Real |

**Current divergences in Atomiverse** (to be aligned):

| current cpu::Frame name | molcrafts name | change needed |
|---|---|---|
| `type_id` | `type` | rename |
| `bond_i`, `bond_j`, `bond_type` | `atomi`, `atomj`, `type` | rename, scope under bonds |
| `angle_i`, `angle_j`, `angle_k`, `angle_type` | `atomi`, `atomj`, `atomk`, `type` | rename, scope under angles |
| `dihedral_i/j/k/l`, `dihedral_type` | `atomi/j/k/l`, `type` | rename, scope under dihedrals |

The flat `cpu::Frame` struct has naming collisions (`type` in atoms vs bonds). Resolution: introduce sub-structs mirroring the block model:

```cpp
struct Frame {
    // atoms block
    struct { Real *x, *y, *z, *vx, *vy, *vz, *fx, *fy, *fz, *mass; Int *type; } atoms;
    // bonds block
    struct { Int *atomi, *atomj, *type; } bonds;
    // angles block
    struct { Int *atomi, *atomj, *atomk, *type; } angles;
    // dihedrals block
    struct { Int *atomi, *atomj, *atomk, *atoml, *type; } dihedrals;
    // system
    Real energy, kinetic, stress[6];
    Mat3 box;
};
```

This is a **prerequisite refactor** — not part of this spec, but required before the Python binding can use uniform naming.

### Storage layout principle

Layout is driven by computation, not convention. Both sides use the same field **names**; the memory **layout** may differ:

| field | molpy layout | Atomiverse layout | rationale |
|-------|-------------|-------------------|-----------|
| `x`, `y`, `z` | separate 1D arrays (SoA) | separate 1D arrays (SoA) | CUDA coalescing; per-dimension PBC ops |
| `vx`, `vy`, `vz` | separate 1D arrays | separate 1D arrays | same pattern as positions |
| `mass`, `type` | 1D arrays | 1D arrays | identical |
| topology | 1D arrays per field | 1D arrays per field | identical |
| `box` | (3,3) matrix | Mat3 / 9-element flat | same data, trivial conversion |

Currently both sides are SoA, so transfer is memcpy per field. If Atomiverse later merges xyz into `(N,3)` or `(3,N)` for batch matrix ops, the transfer layer absorbs the layout conversion — field names stay the same.

## Device model

PyTorch-style: `mp.Frame` carries a `.device` attribute. `.to(device)` returns a new Frame on the target device. Data access on non-CPU devices raises an error — no implicit transfer.

```python
frame = mp.Frame()
frame["atoms"] = mp.Block({"x": x, "y": y, "z": z, "mass": m, "type": t})
frame.box = mp.Box(np.eye(3) * 30.0)
frame.device                         # "cpu"

cuda_frame = frame.to("cuda")       # H2D, returns new mp.Frame(device="cuda")
cuda_frame["atoms"]["x"]            # nb::ndarray<Real, device::cuda> (DLPack)
cuda_frame["atoms"]["x"].numpy()    # RuntimeError: data on cuda

# DLPack to PyTorch (zero-copy, no transfer)
x = torch.from_dlpack(cuda_frame["atoms"]["x"])  # same CUDA memory
```

Devices:

| device | backing storage | field access | DLPack device type |
|--------|----------------|-------------|-------------------|
| `"cpu"` | molrs::Frame (Rust) | numpy arrays | `kDLCPU` |
| `"cuda"` | `cuda::Frame` (device ptrs) | DLPack ndarrays | `kDLCUDA` |

`"torch"` is not a separate device. CUDA-resident DLPack arrays are directly consumable by PyTorch/JAX/CuPy via `from_dlpack()`. This is strictly more general than a dedicated torch device.

## Data flow

### Construction (always CPU, dict-of-blocks)

```python
frame = mp.Frame()
frame["atoms"] = mp.Block({
    "x": x_coords, "y": y_coords, "z": z_coords,
    "mass": masses, "type": types
})
frame["bonds"] = mp.Block({"atomi": bi, "atomj": bj, "type": bt})
frame.box = mp.Box(np.eye(3) * 30.0)
```

### .to("cuda")

```
mp.Frame (cpu, dict-of-blocks)
  │  frame["atoms"]["x"] → numpy (N,)   ← already SoA!
  │  frame["atoms"]["y"] → numpy (N,)
  │  frame["atoms"]["z"] → numpy (N,)
  │  ...
  │
  │  C++: memcpy per field into cpu::Frame (staging buffer, never exposed)
  ▼
cpu::Frame (xtensor, SoA)   [internal, transient]
  │
  │  C++: cpu::to_cuda()
  ▼
cuda::Frame (device ptrs, SoA)
  │
  │  nanobind: wrap device ptrs as nb::ndarray<device::cuda>
  ▼
mp.Frame (cuda)
  │  ._handle = CudaFrameWrapper (Python-ref-counted)
  │  ["atoms"]["x"] → nb::ndarray<Real, device::cuda>  (DLPack-ready)
```

No AoS↔SoA conversion — both sides are SoA. Transfer is memcpy per field.

### .to("cpu")

```
mp.Frame (cuda)
  │  ._handle.cuda_frame
  ▼
cuda::Frame
  │  C++: cuda::to_cpu()
  ▼
cpu::Frame (xtensor, SoA)   [internal, transient]
  │  memcpy per field back into numpy arrays
  ▼
mp.Frame (cpu, dict-of-blocks)
  │  ["atoms"]["x"] → numpy (N,)
```

### Driver.run()

```python
sim = av.MD()
sim.run(frame, 1000)
```

Two paths:

| frame.device | internal behavior |
|---|---|
| `"cpu"` | C++ `Driver::run(cpu::Frame&, n)` — full lifecycle (materialize → run → retrieve) |
| `"cuda"` | C++ `Backend::run(n)` on existing cuda::Frame — no transfer overhead |

When frame is CPU, `run()` writes results back into the same `mp.Frame` in-place (atoms block updated with new positions/velocities/forces, metadata updated with energy/stress).

When frame is CUDA, data stays on device. User calls `.to("cpu")` explicitly to read results.

## Transfer mechanics

Both molpy and Atomiverse use SoA for positions (`x`, `y`, `z` as separate 1D arrays). No layout conversion needed — transfer is memcpy per field.

### CPU transfer (numpy → cpu::Frame → numpy)

```cpp
// Bulk transfer: accept dict of numpy arrays, populate cpu::Frame
// Called from Python with frame["atoms"] fields
void populate_cpu_frame(cpu::Frame& dst,
                        nb::ndarray<Real, nb::ndim<1>, nb::device::cpu, nb::c_contig> x,
                        nb::ndarray<Real, nb::ndim<1>, nb::device::cpu, nb::c_contig> y,
                        nb::ndarray<Real, nb::ndim<1>, nb::device::cpu, nb::c_contig> z,
                        nb::ndarray<Real, nb::ndim<1>, nb::device::cpu, nb::c_contig> mass,
                        nb::ndarray<Int, nb::ndim<1>, nb::device::cpu, nb::c_contig> type_id) {
    Int N = dst.n_atoms();
    std::memcpy(dst.x.data(), x.data(), N * sizeof(Real));
    std::memcpy(dst.y.data(), y.data(), N * sizeof(Real));
    std::memcpy(dst.z.data(), z.data(), N * sizeof(Real));
    std::memcpy(dst.mass.data(), mass.data(), N * sizeof(Real));
    std::memcpy(dst.type_id.data(), type_id.data(), N * sizeof(Int));
}
```

### CUDA DLPack exposure

On CUDA, individual SoA arrays are exposed as DLPack ndarrays. These are **non-owning views** of `cuda::Frame` device memory:

```python
# Access via block-style interface (mirrors molpy)
cuda_frame["atoms"]["x"]   # nb::ndarray<Real, device::cuda>, shape (N,)
cuda_frame["atoms"]["y"]   # nb::ndarray<Real, device::cuda>, shape (N,)

# Direct DLPack → PyTorch (zero-copy)
x = torch.from_dlpack(cuda_frame["atoms"]["x"])
y = torch.from_dlpack(cuda_frame["atoms"]["y"])
z = torch.from_dlpack(cuda_frame["atoms"]["z"])

# Compose (N,3) if model expects it:
pos = torch.stack([x, y, z], dim=1)
```

## nanobind C++ module (`_atomiverse`)

### CpuFrame (internal staging buffer, not user-facing)

Users never construct or see `CpuFrame`. It exists only as a transient shuttle between numpy arrays and `cuda::Frame`.

```cpp
// Not exported to Python as a class — used internally by transfer functions
// CpuFrame is allocated inside to_cuda() and to_cpu(), then discarded.
```

### CudaFrameWrapper (Python-ref-counted CUDA handle)

Wraps `cuda::Frame`. Exposes fields as DLPack ndarrays. Python ref-counting triggers `cuda::destroy()`.

```cpp
struct CudaFrameWrapper {
    cuda::Frame frame;
    ~CudaFrameWrapper() { cuda::destroy(frame); }

    // Block-style field access for Python __getitem__ dispatch
    nb::ndarray<Real, nb::device::cuda> get_real_field(const std::string& field) {
        Real* ptr = resolve_real_ptr(frame, field);  // "x" → frame.x, etc.
        size_t shape[1] = {(size_t)frame.n_atoms()};
        return nb::ndarray<Real, nb::device::cuda>(ptr, 1, shape, nb::handle());
    }

    nb::ndarray<Int, nb::device::cuda> get_int_field(const std::string& field) {
        Int* ptr = resolve_int_ptr(frame, field);
        size_t n = resolve_count(frame, field);  // n_atoms or n_bonds etc.
        size_t shape[1] = {n};
        return nb::ndarray<Int, nb::device::cuda>(ptr, 1, shape, nb::handle());
    }
};
```

### Transfer functions

```cpp
// mp.Frame(cpu) → mp.Frame(cuda)
// Python passes dict-of-dicts: {"atoms": {"x": np, "y": np, ...}, "bonds": {...}}
m.def("frame_to_cuda", [](nb::dict blocks, nb::object box_matrix) -> CudaFrameWrapper* {
    // Extract atoms block (required)
    nb::dict atoms = nb::cast<nb::dict>(blocks["atoms"]);
    auto x = nb::cast<nb::ndarray<Real, nb::ndim<1>, nb::device::cpu, nb::c_contig>>(atoms["x"]);
    Int N = x.shape(0);

    cpu::Frame staging(N);
    // Same names: memcpy per field directly
    std::memcpy(staging.atoms.x, x.data(), N * sizeof(Real));
    std::memcpy(staging.atoms.y, get_field<Real>(atoms, "y"), N * sizeof(Real));
    std::memcpy(staging.atoms.z, get_field<Real>(atoms, "z"), N * sizeof(Real));
    std::memcpy(staging.atoms.mass, get_field<Real>(atoms, "mass"), N * sizeof(Real));
    std::memcpy(staging.atoms.type, get_field<Int>(atoms, "type"), N * sizeof(Int));
    // optional: vx, vy, vz if present
    // topology blocks: bonds, angles, dihedrals — same pattern
    // box ...
    return new CudaFrameWrapper{cpu::to_cuda(staging)};
});

// mp.Frame(cuda) → (blocks_dict, box, scalars)
m.def("frame_to_cpu", [](CudaFrameWrapper& src) -> nb::tuple {
    Int N = src.frame.n_atoms();
    cpu::Frame staging(N);
    cuda::to_cpu(src.frame, staging);

    // Return blocks as dict-of-dicts with same field names
    nb::dict atoms;
    atoms["x"] = make_numpy_copy(staging.atoms.x, N);
    atoms["y"] = make_numpy_copy(staging.atoms.y, N);
    // ... all atom fields ...
    nb::dict blocks;
    blocks["atoms"] = atoms;
    // bonds, angles, dihedrals if present ...

    nb::dict scalars;
    scalars["energy"] = staging.energy;
    scalars["kinetic"] = staging.kinetic;
    // ...
    return nb::make_tuple(blocks, make_box_numpy(staging.box), scalars);
});
```

### Driver + components

```cpp
nb::class_<md::MD>(m, "MD")
    .def(nb::init<>())
    .def("add", &md::MD::add<md::VelocityVerlet>)  // one overload per component
    .def("on", ...)   // returns phase binder
    .def("build", &md::MD::build)
    .def("run_cpu", [](md::MD& sim, cpu::Frame& frame, int n) {
        sim.run(frame, n);
    })
    .def("run_cuda", [](md::MD& sim, CudaFrameWrapper& w, int n) {
        // direct GPU path — no transfer
        sim.run_on_device(w.frame, n);
    })
    ;

// Components — simple value classes, direct binding
nb::class_<md::VelocityVerlet>(m, "VelocityVerlet")
    .def(nb::init<Real>(), nb::arg("dt"));

nb::class_<LJ>(m, "LJ")
    .def(nb::init<Real, Real, Real>(),
         nb::arg("epsilon"), nb::arg("sigma"), nb::arg("cutoff"));

// ... NeighborList, NoseHoover, Thermo, Recorder, SCF components ...
```

## Python package (`atomiverse`)

Two-layer design: thin C++ nanobind layer (`_atomiverse`) + Pythonic OOP wrapper.

**Core principle:** `mp.Frame` is the sole user-facing type. Atomiverse never exposes its own Frame. `mp.Frame ↔ CpuFrame` conversion is internal.

```
python/src/
  atomiverse/
    __init__.py           # public API: MD, RHF, UHF, UKS, components
    _binding.py           # C++ import shim (_atomiverse → _cpp)
    _convert.py           # mp.Frame ↔ CpuFrame conversion (internal)
    md.py                 # MD driver: generic add(), string phases, auto-build
    scf.py                # RHF/UHF/UKS solvers, Integrals, result classes
    _device.py            # device handlers for molpy registry (opt-in)
  bindings/
    bind.cpp              # NB_MODULE(_atomiverse, m) — top-level
    bind_frame.cpp        # CpuFrame (internal), numpy transfer
    bind_md.cpp           # MD driver + MD components (per-type add methods)
    bind_scf.cpp          # SCF integrals + solvers (free functions)
    cpu_frame_wrapper.h   # CpuFrame C++ wrapper struct
```

### Python API

```python
import molpy as mp
import atomiverse as av
import numpy as np

# Data layer: molpy — mp.Frame is the ONLY user-facing type
frame = mp.Frame({"atoms": mp.Block({
    "x": [0.0, 0.74], "y": [0, 0], "z": [0, 0],
    "type": [1, 1], "mass": [1.008, 1.008],
})})
frame.box = mp.Box(np.eye(3) * 100.0)

# Compute layer: atomiverse — accepts mp.Frame, converts internally
sim = av.MD(dt=0.001)
sim.add(av.InitVelocity(temperature=10.0, seed=42))
sim.on("force").add(av.NeighborList(cutoff=8.5), order=-10)
sim.on("force").add(av.LJ(epsilon=0.01, sigma=3.4, cutoff=8.5))
sim.run(frame, 1000)   # mp.Frame in → CpuFrame → compute → write back
# frame["atoms"]["x"] now has updated positions
# frame.metadata["energy"] has the total energy

# SCF — OOP solvers, also accept mp.Frame
rhf = av.RHF(frame, basis="sto-3g")
result = rhf.solve(tol=1e-10)
result.energy           # property
result.eigenvalues      # property (numpy)
result.density_matrix   # property (numpy)
rhf.integrals.S         # overlap matrix

# UHF — charge/multiplicity
uhf = av.UHF(frame, multiplicity=2)
result = uhf.solve()

# DFT
dft = av.UKS(frame, xc="LSDA", multiplicity=1)
result = dft.solve()
result.xc_energy

# Optional: register device handlers for frame.to("cuda")
av.register_molpy()   # raises if molpy protocol not ready
```

### `__init__.py`

```python
from atomiverse._atomiverse import (
    CudaFrameWrapper, frame_to_cuda, frame_to_cpu,
    MD, VelocityVerlet, InitVelocity, NeighborList, LJ, NoseHoover,
    Thermo, Recorder, RDF,
    SCF, Hamiltonian, EigenSolver, DensityBuilder, ConvergenceCheck,
)
from atomiverse._device import register_devices

register_devices()   # plug "cuda" into molpy's device registry
```

### `_device.py` — molpy integration

```python
import numpy as np
import molpy as mp
from atomiverse._atomiverse import frame_to_cuda, frame_to_cpu

class CudaDeviceHandler:
    """Transfers mp.Frame data to/from CUDA via Atomiverse."""

    def from_cpu(self, frame: mp.Frame) -> mp.Frame:
        """mp.Frame(cpu) → mp.Frame(cuda)"""
        # Same field names on both sides — just pass blocks through
        blocks = {}
        for name in ("atoms", "bonds", "angles", "dihedrals"):
            if name in frame:
                blocks[name] = {k: v for k, v in frame[name].items()}

        box = np.asarray(frame.box.matrix) if frame.box else None

        # C++ side: iterate blocks, memcpy per field → cpu::Frame → cuda::Frame
        cuda_handle = frame_to_cuda(blocks, box)

        return mp.Frame(_device="cuda", _handle=cuda_handle)

    def to_cpu(self, frame: mp.Frame) -> mp.Frame:
        """mp.Frame(cuda) → mp.Frame(cpu)"""
        # C++ does: cuda::Frame → cpu::Frame → dict-of-dicts of numpy arrays
        # Returns {"atoms": {"x": np, "y": np, ...}, "bonds": {...}, ...}
        blocks, box, scalars = frame_to_cpu(frame._handle)

        result = mp.Frame()
        for name, fields in blocks.items():
            result[name] = mp.Block(fields)   # same field names, no mapping
        if box is not None:
            result.box = mp.Box(box)
        result.metadata.update(scalars)       # energy, kinetic, stress

        return result


def register_devices():
    mp.register_device("cuda", CudaDeviceHandler())
```

## mp.Frame integration contract

molpy must implement:

```python
class Frame(MutableMapping[str, Block]):
    _blocks: dict[str, Block]
    _device: str = "cpu"
    _handle: Any = None          # None for cpu (data in molrs), CudaFrameWrapper for cuda
    box: Box | None = None
    metadata: dict[str, Any]

    _device_registry: ClassVar[dict] = {}

    @classmethod
    def register_device(cls, name: str, handler) -> None:
        cls._device_registry[name] = handler

    @property
    def device(self) -> str:
        return self._device

    def to(self, device: str) -> "Frame":
        if device == self._device:
            return self
        if device == "cuda":
            return self._device_registry["cuda"].from_cpu(self)
        elif device == "cpu":
            return self._device_registry[self._device].to_cpu(self)
        raise ValueError(f"unknown device: {device}")

    def __getitem__(self, key: str) -> Block:
        if self._device != "cpu":
            # CUDA frame: return DLPack-backed block from _handle
            return self._handle.get_block(key)
        return self._blocks[key]
```

This is the **minimal contract**. molpy defines the protocol and the dict-of-blocks structure; atomiverse implements `CudaDeviceHandler`. molpy has no CUDA dependency.

When `device == "cuda"`, block field access returns `nb::ndarray<device::cuda>` instead of numpy. This preserves the `frame["atoms"]["x"]` API on both devices — only the backing storage differs.

## DLPack interop with PyTorch

Zero-copy path for ML potentials:

```python
import torch
import molpy as mp
import atomiverse as av

frame = mp.Frame()
frame["atoms"] = mp.Block({"x": x, "y": y, "z": z, "mass": m, "type": t})
frame.box = mp.Box(np.eye(3) * 10.0)

cuda_frame = frame.to("cuda")

# Zero-copy: cuda::Frame device ptr → PyTorch tensor
atoms = cuda_frame["atoms"]
x = torch.from_dlpack(atoms["x"])    # (N,) CUDA tensor, same memory
y = torch.from_dlpack(atoms["y"])
z = torch.from_dlpack(atoms["z"])
pos = torch.stack([x, y, z], dim=1)  # (N, 3)

# Forward + backward through ML model
pos.requires_grad_(True)
energy = mlip_model(pos)
forces = -torch.autograd.grad(energy, pos)[0]

# Write forces back (same device memory if non-leaf)
# Or: copy into cuda_frame forces via DLPack
```

No `torch::Frame` needed in the Python path. DLPack subsumes it.

## CMake integration

```cmake
option(ATV_BUILD_PYTHON "Build Python bindings (nanobind)" OFF)

if(ATV_BUILD_PYTHON)
    find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
    # nanobind via FetchContent (recommended) or find_package
    FetchContent_Declare(
        nanobind
        GIT_REPOSITORY https://github.com/wjakob/nanobind.git
        GIT_TAG v2.4.0
    )
    FetchContent_MakeAvailable(nanobind)
    nanobind_add_module(_atomiverse
        src/python/bind.cpp
        src/python/bind_frame.cpp
        src/python/bind_md.cpp
        src/python/bind_scf.cpp
        src/python/bind_potentials.cpp
        src/python/bind_observers.cpp
    )
    target_link_libraries(_atomiverse PRIVATE atomiverse)
    install(TARGETS _atomiverse LIBRARY DESTINATION atomiverse)
endif()
```

Build via scikit-build-core for pip-installable wheels:

```toml
# pyproject.toml
[build-system]
requires = ["scikit-build-core>=0.10", "nanobind>=2.0"]
build-backend = "scikit_build_core.build"

[project]
name = "atomiverse"
requires-python = ">=3.9"
dependencies = ["numpy", "molpy"]
```

## Files to create

| Path | Purpose |
|------|---------|
| `src/python/bind.cpp` | `NB_MODULE(_atomiverse, m)` top-level module |
| `src/python/bind_frame.cpp` | CpuFrame, CudaFrame, SoA↔AoS, transfer |
| `src/python/bind_md.cpp` | MD driver + component bindings |
| `src/python/bind_scf.cpp` | SCF driver + component bindings |
| `src/python/bind_potentials.cpp` | LJ, HarmonicBond, HarmonicAngle, ... |
| `src/python/bind_observers.cpp` | Thermo, Recorder, RDF, ... |
| `python/atomiverse/__init__.py` | Public API + device registration |
| `python/atomiverse/_device.py` | CudaDeviceHandler (mp.Frame ↔ cuda) |
| `python/atomiverse/md.py` | MD convenience wrappers (optional) |
| `python/atomiverse/scf.py` | SCF convenience wrappers (optional) |
| `pyproject.toml` | scikit-build-core build config |
| `CMakeLists.txt` | Add `ATV_BUILD_PYTHON` option |

## Dependency graph

```
molpy  (runtime, pip)
  └── molrs (Rust, PyO3)

atomiverse (pip wheel, scikit-build-core)
  ├── _atomiverse.so (nanobind C++ extension)
  │   └── libatomiverse (C++/CUDA static lib)
  │       └── molrs_cxx (CXX FFI, trajectory I/O — NOT for Frame transfer)
  ├── numpy (runtime)
  └── molpy (runtime, for device registration)
```

Key: the Frame data bridge is **numpy (Python level)**, not molrs CXX FFI (C++ level). Two independent FFI paths coexist without conflict:
- molrs CXX FFI: C++ trajectory I/O (MolRec, XYZ, Zarr) — existing, unchanged
- numpy/DLPack: Python Frame data exchange — new, no Rust involvement

## Testing strategy

| Level | Test | Verifies |
|-------|------|----------|
| Unit | `import atomiverse` succeeds | nanobind module loads |
| Unit | `CpuFrame(100).n_atoms == 100` | basic construction |
| Round-trip | numpy → CpuFrame → CudaFrame → CpuFrame → numpy | data integrity within tolerance |
| DLPack | `torch.from_dlpack(cuda_frame.x).data_ptr() == cuda_frame.x.data_ptr()` | zero-copy (same pointer) |
| Layout | scatter (N,3) → SoA → gather → (N,3) matches original | AoS↔SoA correctness |
| Device | `frame.to("cuda").device == "cuda"` | device tracking |
| Driver | MD 100-step via Python matches C++ reference energy | end-to-end correctness |
| SCF | RHF H2 via Python matches -1.1168 Ha | SCF correctness |
| Lifecycle | CudaFrame ref-count → 0 triggers cuda::destroy() | no GPU memory leak |
| Topology | bonds/angles/dihedrals survive round-trip | topology transfer |

## Open questions

1. ~~**mp.Frame field names**~~ **Resolved.** molpy uses `x`, `y`, `z` (SoA), `mass`, `type` in atoms block. Topology: `atomi`, `atomj`, `atomk`, `atoml`, `type` in bonds/angles/dihedrals blocks. Field mapping table above aligns with naming conventions doc.

2. ~~**mp.Frame topology API**~~ **Resolved.** Dict-of-blocks: `frame["bonds"]`, `frame["angles"]`, `frame["dihedrals"]` each a Block with `atomi`, `atomj`, etc.

3. **Autograd write-back**: After `torch.autograd.grad()` computes forces, what's the cleanest path to write them back into `cuda::Frame.fx/fy/fz`? DLPack views from `nb::ndarray` should be writable (non-const pointer). If `torch.from_dlpack()` produces a view of the original device memory, writes to that tensor propagate directly. Need to verify: does `torch.from_dlpack()` guarantee a writable view, or can it copy?

    > As zero-copy as possible is ideal, but we may need a small CUDA kernel if PyTorch tensors aren't guaranteed to be writable views of the original memory.

4. **CPU-only mode**: When CUDA is not available (`ATV_BUILD_CUDA=OFF`), `.to("cuda")` should raise a clear error. The nanobind module should still load for CPU-only workflows.

    > `raise RuntimeError("Atomiverse was built without CUDA support")` in `to_cuda()` if CUDA is unavailable.

5. **mp.Frame device protocol**: molpy currently has no `_device`, `_handle`, `to()`, or `register_device()`. This protocol needs to be added to molpy. Spec defines the contract; molpy implementation is a separate task.

6. **CUDA Block access**: When `frame.device == "cuda"`, `frame["atoms"]["x"]` should return a DLPack ndarray, not a numpy array. This requires either (a) a CudaBlock wrapper class that mimics Block but returns device ndarrays, or (b) mp.Frame dispatches `__getitem__` to the _handle. Design choice to be made in the molpy protocol spec.

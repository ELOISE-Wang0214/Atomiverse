# Tutorial 4: The LibTorch Backend -- Differentiable Simulation

This tutorial is for readers who know basic C++ and have a passing familiarity
with GPU programming, but have **no background in automatic differentiation or
PyTorch**. It assumes you have read the earlier tutorials or at least understand
that Atomiverse runs molecular simulations using `cpu::Frame` (host arrays) and
`cuda::Frame` (GPU device pointers).

The goal is to explain what the third backend -- `ATV::torch` -- adds, why you
would use it, and how every piece of its API works.

---

## 1. Introduction

### The problem: simulation parameters are not free

When you run a molecular dynamics simulation, the physics comes from the force
field: a set of equations and parameters that tell you how atoms push and pull on
each other. The Lennard-Jones potential, for example, uses two parameters:
$\varepsilon$ (the depth of the energy well, in eV) and $\sigma$ (the
collision diameter, in Angstrom). If these parameters are wrong, the simulation
produces wrong energies, wrong pressures, wrong melting points.

Where do the right parameters come from? Historically, a domain expert tunes
them by hand or fits them to a table of experimental measurements. That works,
but it is slow. Modern approaches use **gradient-based optimization**: define a
loss function that measures how far the simulation output is from a reference
(quantum-mechanical calculation, experimental data), then use calculus to
compute the derivative of that loss with respect to every parameter, and follow
the gradient downhill.

The problem is that a simulation is a long chain of floating-point operations.
Computing $\partial E / \partial \varepsilon$ by finite differences (run twice,
perturb $\varepsilon$, compare) requires two full simulations per parameter.
For a force field with hundreds of parameters, this is prohibitively expensive.

### The solution: automatic differentiation

Automatic differentiation (AD) is a technique from numerical computing that
records every arithmetic operation during a forward computation and then
replays them in reverse to produce exact derivatives. It is the engine behind
all modern deep-learning frameworks, including PyTorch. If you can express a
computation as a sequence of differentiable operations, AD gives you the
gradient of any scalar output with respect to any input -- in a single backward
pass, regardless of how many parameters there are.

### What the torch backend provides

Atomiverse's `ATV::torch` backend wraps the existing CUDA simulation kernels
in LibTorch's autograd graph. The forward pass calls the same GPU kernels that
the `cuda::` backend uses, so the numerical output is identical to machine
precision. The backward pass provides analytically-derived gradient kernels that
LibTorch's AD engine calls automatically when you invoke `.backward()`.

This enables three workflows that were previously impossible inside Atomiverse:

1. **Parameter fitting.** Given a reference energy $E_\text{ref}$, define
   $\mathcal{L} = (E - E_\text{ref})^2$, call `L.backward()`, and read
   $\partial \mathcal{L} / \partial \varepsilon$ directly.

2. **Force matching.** Given reference forces $F_\text{ref}$ from a quantum
   calculation, compute forces from a parameterized potential, and minimize
   $\sum_i \|F_i - F_{i,\text{ref}}\|^2$ with gradient descent.

3. **ML interop.** A neural-network potential produces energy and forces as
   `at::Tensor` objects. The torch backend lets those tensors flow seamlessly
   into the same pair lists and topology structures that Atomiverse already uses.

### What the torch backend does not do

The torch backend is **not** a simulation driver. It does not run MD time steps
or SCF iterations. It provides differentiable *operators* -- individual energy
and force evaluations -- that you call from your own training loop or fitting
script. Simulation drivers (`md::MD`, `scf::SCF`, `minimizer::SD`) remain in
the `cuda::` backend.

---

## 2. Build Setup

### CMake flags

The torch backend is gated behind a single CMake option. When it is `OFF`
(the default), no LibTorch headers are included, no LibTorch libraries are
linked, and the core engine compiles identically to before.

```bash
# Standard build (no torch)
cmake -B build -DATV_BUILD_TESTS=ON

# Build with LibTorch backend
cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_TORCH=ON \
      -DCMAKE_PREFIX_PATH=/path/to/libtorch
```

The `CMAKE_PREFIX_PATH` must point to the root of your LibTorch installation
(the directory containing `share/cmake/Torch/TorchConfig.cmake`). You can
download a CUDA-enabled LibTorch from https://pytorch.org/get-started/locally/.

### What happens in CMakeLists.txt

When `ATV_BUILD_TORCH=ON`, the build system does four things:

1. Calls `find_package(Torch REQUIRED)` to locate LibTorch.
2. Adds the torch backend source files to the `atomiverse` library:
   `src/torch/transfer.cu`, `src/torch/stream_guard.cu`,
   `src/torch/md/lj.cu`, `src/torch/md/neighbor_list.cu`,
   `src/torch/bonded/harmonic_bond.cu`.
3. Links `${TORCH_LIBRARIES}` and adds `${TORCH_INCLUDE_DIRS}`.
4. Defines the preprocessor macro `ATV_HAS_TORCH` on the `atomiverse` target.

All torch code -- headers and sources -- is wrapped in `#ifdef ATV_HAS_TORCH`.
If you build without the flag, these files are syntactically invisible.

### Build matrix

| `ATV_BUILD_TORCH` | LibTorch required? | `ATV_HAS_TORCH` defined? | Torch tests built? |
|---|---|---|---|
| `OFF` (default) | No | No | No |
| `ON` | Yes | Yes | Yes |

Both `ATV_REAL=float` and `ATV_REAL=double` are supported. The compile-time
constant `ATV::torch::kRealDtype` maps `float` to `ScalarType::Float` and
`double` to `ScalarType::Double`.

### Version compatibility

LibTorch must be compiled against the same CUDA toolkit version as Atomiverse.
A mismatch (e.g., LibTorch built with CUDA 12.1, Atomiverse built with CUDA
11.8) will produce linker errors or runtime crashes. Check with:

```bash
python -c "import torch; print(torch.version.cuda)"   # LibTorch's CUDA
nvcc --version                                         # Your toolkit
```

---

## 3. Architecture

### Three backends

Atomiverse now has three backends, each with the same structural pattern -- a
Frame type, a Context type, and transfer functions:

```
Facade (Driver, Plan, Stage)
  |
  +-- cpu::              host device (xtensor)
  |     cpu::Frame       xt::xtensor arrays
  |     cpu::Context     (empty for now)
  |
  +-- cuda::             GPU device (raw pointers)
  |     cuda::Frame      Real*, Int* device pointers
  |     cuda::Context    CUDA stream handles
  |
  +-- torch::            differentiable device (LibTorch)
        torch::Frame     at::Tensor (CUDA, with autograd)
        torch::Context   autograd configuration
```

The three backends are **fully decoupled**. No `cpu::` header includes a
`cuda::` or `torch::` header. No `cuda::` header includes a `torch::` header.
Transfer functions are the only bridge, and they live in their own headers.

### Frame comparison

The table below shows how the same physical quantity is represented in each
backend. All three store the same data; only the container type differs.

| Quantity | `cpu::Frame` | `cuda::Frame` | `torch::Frame` |
|---|---|---|---|
| Position x | `xt::xtensor<Real,1>` | `Real*` (device) | `at::Tensor` (CUDA) |
| Force fx | `xt::xtensor<Real,1>` | `Real*` (device) | `at::Tensor` (CUDA) |
| Energy | `Real` (scalar) | `Real*` (device, 1 elem) | `at::Tensor` (CUDA, shape [1]) |
| Pair indices | not stored | `PairList` (Int* device) | `at::Tensor` (CUDA, Int) |
| Box | `Mat3` (host) | `Real*` (device, 9 elem) | `at::Tensor` (CUDA, shape [9]) |
| Memory management | xtensor RAII | manual `cudaMalloc`/`cudaFree` | Torch reference counting |
| Autograd | N/A | N/A | `requires_grad` per tensor |
| Atom count | `n_atoms()` (host) | `n_atoms()` (host-cached) | `n_atoms()` (host-cached) |

The distinguishing property of `torch::Frame` is that position and parameter
tensors can carry `requires_grad=true`. When a differentiable operator reads
these tensors, LibTorch records the operation in a computational graph. Calling
`.backward()` on the output then walks that graph in reverse to produce
gradients.

### Namespace aliasing

The namespace `ATV::torch` collides with the global `::torch` namespace that
LibTorch itself defines. Inside the `ATV::torch` namespace, unqualified names
like `Tensor` resolve to `ATV::torch::Tensor` (which is aliased to
`::torch::Tensor`). Outside `ATV::torch`, you must write `ATV::torch::Tensor`
or use the full `::torch::Tensor`.

The standard pattern in tests and user code is a namespace alias:

```cpp
namespace atorch = ATV::torch;
```

This avoids the collision entirely. You write `atorch::Frame`,
`atorch::lj()`, `atorch::Tensor` without ambiguity. The project's test files
all follow this convention.

### Type aliases inside `ATV::torch`

The `torch/frame.cuh` header defines four aliases so that code inside
`ATV::torch` never needs to write `::torch::`:

```cpp
namespace ATV::torch {
    using Tensor     = ::torch::Tensor;
    using ScalarType = ::at::ScalarType;
    using Device     = ::torch::Device;
    using Stream     = ::c10::cuda::CUDAStream;
}
```

It also defines two compile-time constants that map Atomiverse's `Real` and
`Int` to LibTorch scalar types:

```cpp
inline constexpr ScalarType kRealDtype =
    std::is_same_v<Real, float> ? ScalarType::Float : ScalarType::Double;

inline constexpr ScalarType kIntDtype =
    std::is_same_v<Int, int> ? ScalarType::Int : ScalarType::Long;
```

These constants ensure that every tensor created by the torch backend has the
correct element type, matching the compile-time precision of the rest of
Atomiverse.

---

## 4. Quick Start: Differentiable LJ

This section is a complete, standalone example. It creates a two-atom system,
builds a pair list, computes the Lennard-Jones energy through the torch backend,
calls `.backward()`, and reads the gradient of the energy with respect to
$\varepsilon$.

### The physics in one paragraph

The Lennard-Jones (LJ) 12-6 potential describes the interaction between two
neutral atoms separated by a distance $r$:

$$
V(r) = 4\varepsilon \left[\left(\frac{\sigma}{r}\right)^{12}
      - \left(\frac{\sigma}{r}\right)^{6}\right] - V_\text{shift}
$$

where $\varepsilon$ is the well depth (in eV), $\sigma$ is the collision
diameter (in Angstrom), and $V_\text{shift} = V(r_\text{cut})$ is subtracted so
that the potential is exactly zero at the cutoff distance $r_\text{cut}$. For
two argon-like atoms at $r = 1.5$ Angstrom with $\varepsilon = 1.0$ eV,
$\sigma = 1.0$ Angstrom, and $r_\text{cut} = 3.0$ Angstrom, the energy is
approximately $-0.3881$ eV.

### Step-by-step code

```cpp
// Build: cmake -B build -DATV_BUILD_TORCH=ON -DCMAKE_PREFIX_PATH=/path/to/libtorch
// Run:   ./build/examples/torch_lj

#include "atomiverse/torch/md/lj.cuh"
#include "atomiverse/torch/frame.cuh"

#include <cuda_runtime.h>
#include <iostream>

namespace atorch = ATV::torch;
using ATV::Real;
using ATV::Int;

int main() {
    // ---- 1. Tensor options ----
    // Every tensor needs a dtype and a device. We create two option bundles:
    // one for floating-point data (positions, forces, parameters) and one for
    // integer data (atom indices, pair indices).
    auto real_opts = ::torch::TensorOptions()
        .dtype(atorch::kRealDtype)
        .device(::torch::kCUDA);
    auto int_opts = ::torch::TensorOptions()
        .dtype(atorch::kIntDtype)
        .device(::torch::kCUDA);

    // ---- 2. Build the Frame ----
    // Two atoms on the x-axis: atom 0 at x=0, atom 1 at x=1.5 Angstrom.
    Real r = 1.5;

    atorch::Frame frame;
    frame.n_atoms_ = 2;
    frame.x  = ::torch::tensor({Real(0), r}, real_opts);
    frame.y  = ::torch::zeros({2}, real_opts);
    frame.z  = ::torch::zeros({2}, real_opts);
    frame.fx = ::torch::zeros({2}, real_opts);
    frame.fy = ::torch::zeros({2}, real_opts);
    frame.fz = ::torch::zeros({2}, real_opts);
    frame.energy = ::torch::zeros({1}, real_opts);

    // ---- 3. Populate the pair list ----
    // In a real workflow you would call atorch::neighbor_list(). Here we
    // construct the single pair (0,1) manually.
    frame.n_pairs_ = 1;
    frame.pair_i    = ::torch::tensor({Int(0)}, int_opts);
    frame.pair_j    = ::torch::tensor({Int(1)}, int_opts);
    frame.pair_dx   = ::torch::tensor({r}, real_opts);
    frame.pair_dy   = ::torch::zeros({1}, real_opts);
    frame.pair_dz   = ::torch::zeros({1}, real_opts);
    frame.pair_dist = ::torch::tensor({r}, real_opts);

    // ---- 4. Create differentiable parameters ----
    // epsilon and sigma are scalar CUDA tensors with requires_grad=true.
    // This tells LibTorch to track all operations involving them.
    auto epsilon = ::torch::tensor({Real(1.0)}, real_opts).requires_grad_(true);
    auto sigma   = ::torch::tensor({Real(1.0)}, real_opts).requires_grad_(true);
    Real cutoff  = 3.0;

    // ---- 5. Forward pass ----
    // atorch::lj() calls the same CUDA kernel as cuda::LJ, but wraps it in
    // an autograd::Function so that backward() works.
    atorch::Tensor energy = atorch::lj(frame, epsilon, sigma, cutoff);
    cudaDeviceSynchronize();

    Real E = energy.to(::torch::kCPU).item<Real>();
    std::cout << "Energy: " << E << " eV\n";

    // ---- 6. Backward pass ----
    // One call to .backward() propagates gradients through the entire graph.
    energy.backward();
    cudaDeviceSynchronize();

    Real dE_deps = epsilon.grad().to(::torch::kCPU).item<Real>();
    Real dE_dsig = sigma.grad().to(::torch::kCPU).item<Real>();

    std::cout << "dE/d(epsilon): " << dE_deps << "\n";
    std::cout << "dE/d(sigma):   " << dE_dsig << "\n";

    return 0;
}
```

### What each section does

**Section 1** creates `TensorOptions` objects. These are the LibTorch equivalent
of specifying `float` vs `double` and `CPU` vs `CUDA`. The constant
`atorch::kRealDtype` ensures the tensor dtype matches the compile-time `Real`.

**Section 2** constructs a `torch::Frame` by assigning tensors to each field.
Unlike `cpu::Frame`, which takes an atom count in its constructor and
preallocates arrays, `torch::Frame` is a plain struct -- you fill in the tensors
yourself. The host-cached metadata `n_atoms_` must be set manually.

**Section 3** populates the pair list. In this minimal example there is one
pair: atoms 0 and 1. The displacement vector is $(1.5, 0, 0)$ and the distance
is $1.5$ Angstrom. In a real workflow, `atorch::neighbor_list()` does this
automatically.

**Section 4** creates the LJ parameters as tensors with `requires_grad_(true)`.
This is the key step that enables differentiation. Any tensor without
`requires_grad` is treated as a constant by autograd.

**Section 5** calls `atorch::lj()`. Internally, this function:
1. Extracts raw `Real*` pointers from each tensor via `.data_ptr<Real>()`.
2. Launches `cuda::kernels::lj_force` -- the exact same `__global__` kernel
   that the non-differentiable CUDA backend uses.
3. Registers the operation in LibTorch's autograd graph so that `.backward()`
   knows how to compute gradients.

**Section 6** calls `energy.backward()`, which triggers the reverse pass. The
autograd engine calls `LJFunction::backward()`, which computes
$\partial E / \partial \varepsilon$ and $\partial E / \partial \sigma$ using
analytically-derived gradient kernels, and stores the results in
`epsilon.grad()` and `sigma.grad()`.

### Expected output

```
Energy: -0.388051 eV
dE/d(epsilon): -0.388051
dE/d(sigma):   3.06558
```

The derivative $\partial E / \partial \varepsilon$ equals
$E / \varepsilon = -0.3881 / 1.0 = -0.3881$ because the LJ potential is
linear in $\varepsilon$. The derivative with respect to $\sigma$ involves
the full analytical expression for the shifted LJ pair gradient.

---

## 5. Transfer Functions

The torch backend connects to the other two backends through five transfer
functions, all in the `ATV::torch` namespace.

### `to_torch(cuda::Frame)` -- zero-copy wrap

```cpp
ATV::torch::Frame to_torch(const ATV::cuda::Frame& src);
```

This is the fast path. It wraps each `Real*` device pointer in the
`cuda::Frame` as an `at::Tensor` using `torch::from_blob`. No data is copied.
The returned tensors are non-owning views: they borrow the memory from the
`cuda::Frame`. The caller must ensure the `cuda::Frame` outlives the
`torch::Frame`.

This function is designed for cross-device components that need to call a torch
operator on data that is already resident on the GPU. The overhead is just a
handful of `from_blob` calls -- no GPU memory allocation, no DMA.

### `to_cuda(torch::Frame, cuda::Frame)` -- pointer extraction

```cpp
void to_cuda(const ATV::torch::Frame& src, ATV::cuda::Frame& dst);
```

The reverse of `to_torch`. It extracts the raw `Real*` from each tensor via
`.data_ptr<Real>()` and stores it in the corresponding field of the
`cuda::Frame`. Before extracting, it validates that every tensor is CUDA-
resident, contiguous, and has the correct dtype. If any check fails, it throws
`std::invalid_argument`.

This function is non-owning in the other direction: the `cuda::Frame` borrows
pointers from the `torch::Frame`. The `torch::Frame`'s tensors must remain
alive.

### `from_cpu(cpu::Frame)` -- host-to-device copy

```cpp
ATV::torch::Frame from_cpu(const ATV::cpu::Frame& src);
```

Allocates new CUDA tensors and copies data from the `cpu::Frame`'s xtensor
arrays. This is a full H2D transfer. The returned `torch::Frame` owns its own
memory (managed by LibTorch's CUDA caching allocator).

Use this to bootstrap a standalone torch workflow from user-facing `cpu::Frame`
data.

### `to_cpu(torch::Frame, cpu::Frame)` -- device-to-host copy

```cpp
void to_cpu(const ATV::torch::Frame& src, ATV::cpu::Frame& dst);
```

Copies tensor data back to host xtensor arrays. This is a full D2H transfer.
The destination `cpu::Frame` must already be allocated with the correct number
of atoms.

### `destroy(torch::Frame)` -- release all tensors

```cpp
void destroy(ATV::torch::Frame& frame);
```

Resets every tensor field to an empty `Tensor()` and zeros the metadata counts.
Because LibTorch uses reference counting, this releases the GPU memory as soon
as no other tensor shares the allocation.

### Ownership summary

| Function | Direction | Copy? | Who owns memory after? |
|---|---|---|---|
| `to_torch` | cuda -> torch | No (zero-copy) | `cuda::Frame` owns; `torch::Frame` borrows |
| `to_cuda` | torch -> cuda | No (pointer extraction) | `torch::Frame` owns; `cuda::Frame` borrows |
| `from_cpu` | cpu -> torch | Yes (H2D) | `torch::Frame` owns (Torch allocator) |
| `to_cpu` | torch -> cpu | Yes (D2H) | `cpu::Frame` owns (xtensor) |
| `destroy` | -- | -- | Memory released |

### Round-trip example

```cpp
namespace atorch = ATV::torch;

// Start with user data on CPU
ATV::cpu::Frame cpu_frame(4);
cpu_frame.x = xt::xtensor<Real, 1>({1.0, 2.0, 3.0, 4.0});
cpu_frame.box = Mat3::cubic(10.0);

// Transfer to torch backend (allocates CUDA tensors)
atorch::Frame tframe = atorch::from_cpu(cpu_frame);

// ... perform differentiable operations on tframe ...

// Transfer results back to CPU
ATV::cpu::Frame result(4);
atorch::to_cpu(tframe, result);

// result.x[0] == 1.0, result.x[3] == 4.0
```

### Zero-copy round-trip through `cuda::Frame`

```cpp
// Start from an existing cuda::Frame (e.g., inside a cross-device component)
ATV::cuda::Frame cframe = ATV::cpu::to_cuda(cpu_frame);

// Zero-copy wrap
atorch::Frame tframe = atorch::to_torch(cframe);
// tframe.x.data_ptr<Real>() == cframe.x  (same GPU address)

// Extract back
ATV::cuda::Frame cframe2;
atorch::to_cuda(tframe, cframe2);
// cframe2.x == cframe.x  (same GPU address again)

// Cleanup: destroy the cuda::Frame's allocation (not the torch::Frame's)
ATV::cuda::destroy(cframe);
```

---

## 6. Differentiable Operators

### What "differentiable" means in practice

An operator is **differentiable** if it participates in LibTorch's autograd
graph. Concretely, this means two things:

1. The forward pass is wrapped in a `torch::autograd::Function` subclass. When
   called, LibTorch records the inputs, outputs, and a pointer to the backward
   function.

2. A backward function exists that, given the gradient of the output (a scalar
   $\partial \mathcal{L} / \partial E$), produces the gradients of all inputs
   that have `requires_grad=true`.

### LJ 12-6 (`ATV::torch::lj`)

```cpp
Tensor lj(Frame& frame, const Tensor& epsilon, const Tensor& sigma, Real cutoff);
```

The forward pass zeros `frame.energy`, then launches `cuda::kernels::lj_force`
-- the same `__global__` kernel used by the CUDA backend's `LJ` component.
This kernel loops over all pairs in the pair list, computes the LJ energy and
force for each pair, and accumulates the results via `atomicAdd`. The function
returns `frame.energy` as a tensor suitable for `.backward()`.

The backward pass computes three gradients:

**Gradient with respect to positions** ($\partial E / \partial x_i$). Since
the force is the negative gradient of the energy,
$\partial E / \partial x_i = -f_i$, this gradient is already available from
the forward pass. The backward function simply scales the forces by the
upstream gradient and flips the sign.

**Gradient with respect to** $\varepsilon$. The LJ potential is linear in
$\varepsilon$, so $\partial E / \partial \varepsilon$ is computed by
summing over all pairs:

$$
\frac{\partial E}{\partial \varepsilon}
= \sum_{\text{pairs}} \left[4\left(\frac{\sigma}{r}\right)^{12}
  - 4\left(\frac{\sigma}{r}\right)^6\right]
  - \frac{V_\text{shift}}{\varepsilon}
$$

**Gradient with respect to** $\sigma$. This requires differentiating each
power of $\sigma / r$:

$$
\frac{\partial E}{\partial \sigma}
= \sum_{\text{pairs}} 4\varepsilon
  \left[\frac{12\,\sigma^{11}}{r^{12}}
  - \frac{6\,\sigma^5}{r^6}\right]
$$

Both parameter gradients are computed by a dedicated CUDA kernel
(`lj_param_grad`) that runs a parallel reduction over pairs, using shared memory
for the block-level sum and `atomicAdd` for the final accumulation.

### Harmonic Bond (`ATV::torch::harmonic_bond`)

```cpp
Tensor harmonic_bond(Frame& frame, const Tensor& K, const Tensor& r0);
```

The forward pass launches `cuda::kernels::harmonic_bond_k`, which computes:

$$
E = \sum_{\text{bonds}} K_t \,(r - r_{0,t})^2
$$

where $t$ is the bond type index, $K_t$ is the spring constant (eV /
Angstrom$^2$), and $r_{0,t}$ is the equilibrium distance (Angstrom). The
displacement uses minimum image convention (orthorhombic box).

The backward pass computes gradients with respect to three sets of inputs:

$$
\frac{\partial E}{\partial x_i} = 2 K_t \frac{r - r_{0,t}}{r}\,\Delta x
\qquad
\frac{\partial E}{\partial K_t} = (r - r_{0,t})^2
\qquad
\frac{\partial E}{\partial r_{0,t}} = -2 K_t (r - r_{0,t})
$$

The backward kernel runs one thread per bond, using `atomicAdd` to accumulate
per-atom position gradients and per-type parameter gradients.

### Neighbor List (`ATV::torch::neighbor_list`) -- non-differentiable

```cpp
void neighbor_list(Frame& frame, Real cutoff, Real skin = Real(-1));
```

This function populates the pair list fields (`pair_i`, `pair_j`, `pair_dx`,
`pair_dy`, `pair_dz`, `pair_dist`) in the frame. It delegates to the existing
`cuda::` backend's `md::NeighborList` component internally: it creates a
temporary `cuda::Frame`, runs the full CUDA neighbor-list pipeline (Morton
sort, AABB, grid CSR, fused pair search), then clones the results into
self-contained tensors on the `torch::Frame`.

The neighbor list is **not differentiable**. All returned tensors have
`requires_grad=false`. Gradients through discrete pair membership (which atoms
are neighbors) have no physical meaning -- this is a combinatorial operation,
not a smooth function.

The `skin` parameter adds extra distance beyond the cutoff to avoid rebuilding
the pair list every step. A negative value (the default) uses `0.3 * cutoff`
automatically.

### What is differentiable and what is not

| Operator | Differentiable? | Differentiable w.r.t. |
|---|---|---|
| `lj` | Yes | positions, $\varepsilon$, $\sigma$ |
| `harmonic_bond` | Yes | positions, $K$, $r_0$ |
| `neighbor_list` | No | -- |

This distinction is fundamental: the pair list is *infrastructure* that feeds
into differentiable operators, but is not itself differentiable. You build the
pair list once, then call `lj()` (which is differentiable) on the result.

---

## 7. Testing Patterns

### Namespace aliasing

Every test file begins with:

```cpp
namespace atorch = ATV::torch;
```

This eliminates the `ATV::torch` vs `::torch` collision in test code.
References to raw LibTorch always use the `::torch::` prefix:

```cpp
auto real_opts = ::torch::TensorOptions().dtype(atorch::kRealDtype).device(::torch::kCUDA);
```

### The gradcheck pattern

The standard way to validate a backward pass is to compare the autograd
gradient against a finite-difference approximation. The pattern appears
throughout the test suite:

```cpp
TEST(TorchLJ, GradcheckEpsilon) {
    Real eps_val = 0.5, sig_val = 1.0, cutoff = 3.0, r = 1.5;
    auto real_opts = ::torch::TensorOptions()
        .dtype(atorch::kRealDtype).device(::torch::kCUDA);

    // 1. Finite-difference reference
    Real h = 1e-4;
    Real numerical = (analytical_lj_energy(r, eps_val + h, sig_val, cutoff)
                    - analytical_lj_energy(r, eps_val - h, sig_val, cutoff))
                   / (2.0 * h);

    // 2. Autograd gradient
    auto frame = make_lj_frame(r);
    auto t_eps = ::torch::tensor({eps_val}, real_opts).requires_grad_(true);
    auto t_sig = ::torch::tensor({sig_val}, real_opts);

    auto energy = atorch::lj(frame, t_eps, t_sig, cutoff);
    energy.backward();
    cudaDeviceSynchronize();

    Real autograd = t_eps.grad().to(::torch::kCPU).item<Real>();

    // 3. Compare
    EXPECT_NEAR(autograd, numerical, 1e-2);
}
```

The three steps are always the same:

1. Compute the expected gradient using a known analytical formula and central
   finite differences with step size $h$.
2. Run the forward pass with `requires_grad_(true)` on the parameter of
   interest, then call `.backward()`.
3. Compare `param.grad()` against the finite-difference result.

The tolerance depends on the step size $h$ and the magnitude of higher-order
derivatives. Typical values are `1e-2` for LJ parameters (where the potential
has large second derivatives near the minimum) and `1e-4` for smoother
potentials.

### Analytical validation

For the forward pass, tests compare against hand-computed values. For example,
the shifted LJ energy for two atoms at distance $r = 1.5$ Angstrom with
$\varepsilon = 1.0$ eV, $\sigma = 1.0$ Angstrom, $r_\text{cut} = 3.0$
Angstrom is:

$$
V = 4 \cdot 1.0 \left[\left(\frac{1.0}{1.5}\right)^{12}
  - \left(\frac{1.0}{1.5}\right)^{6}\right]
  - 4 \cdot 1.0 \left[\left(\frac{1.0}{3.0}\right)^{12}
  - \left(\frac{1.0}{3.0}\right)^{6}\right]
\approx -0.3881 \text{ eV}
$$

This analytical reference function is defined as a helper:

```cpp
Real analytical_lj_energy(Real r, Real eps, Real sig, Real cutoff) {
    Real sr_cut = sig / cutoff;
    Real sr_cut6 = sr_cut * sr_cut * sr_cut * sr_cut * sr_cut * sr_cut;
    Real sr_cut12 = sr_cut6 * sr_cut6;
    Real v_shift = 4.0 * eps * (sr_cut12 - sr_cut6);

    Real sr = sig / r;
    Real sr6 = sr * sr * sr * sr * sr * sr;
    Real sr12 = sr6 * sr6;
    return 4.0 * eps * (sr12 - sr6) - v_shift;
}
```

### Edge case: zero pairs

Every operator must handle the degenerate case of zero interactions. The test
creates a frame with `n_pairs_ = 0` (or `n_bonds_ = 0`) and verifies that
the energy is exactly zero and no GPU kernel is launched (or, if launched, it
returns immediately without writing).

### The `#ifdef` guard

Every test file wraps its body in `#ifdef ATV_HAS_TORCH ... #else ... #endif`.
When building without LibTorch, a single `GTEST_SKIP()` test replaces all the
real tests, ensuring the test binary still compiles and reports a clear skip
message.

---

## 8. What's Next

The current implementation (Phase 1) delivers the foundational infrastructure:
`Frame`, `Context`, transfer functions, `StreamGuard`, and three operators (LJ,
HarmonicBond, NeighborList). The following phases are planned:

### Phase 2: ML Force Field Integration

**MLIPForce component.** A cross-device component that loads a TorchScript
model file (`.pt`), wraps `cuda::Frame` data into `torch::Frame` tensors at
each force evaluation step, runs `model.forward()` to get energy and forces,
and writes the results back. This is the main path for running neural-network
potentials inside an Atomiverse MD simulation:

```cpp
md::MD sim;
sim.add(md::VelocityVerlet{0.001});
sim.on(md::Phase::force)
   .add(NeighborList{5.0}, -10)
   .add(ATV::torch::MLIPForce{"model.pt"});
sim.run(frame, 1000);
```

**HarmonicAngle and MultiHarmonicDihedral operators.** The same pattern as
HarmonicBond: forward reuses the existing CUDA kernel, backward provides
analytical gradients.

### Phase 3: Python Bindings and Training Workflows

**Python bindings.** Expose `torch::Frame`, transfer functions, and
differentiable operators to Python via `pybind11`. This connects Atomiverse's
simulation engine to PyTorch's training infrastructure (optimizers, data
loaders, learning-rate schedulers).

**Differentiable simulation experiments.** Longer-range goal: differentiate
through multiple MD steps for end-to-end training of force fields from
trajectory-level loss functions.

### Kernel reuse guarantee

All future differentiable operators follow the same pattern established in
Phase 1: the forward pass calls the existing `cuda::kernels::` function
(pointer-level, backend-agnostic), and the backward pass adds new code in
`src/torch/`. No numerical computation is duplicated between backends.

---

## Appendix A: `torch::Frame` Complete Reference

```cpp
namespace ATV::torch {
struct Frame {
    // Per-atom [N], kRealDtype, CUDA
    Tensor x, y, z;              // Positions [Angstrom]
    Tensor vx, vy, vz;           // Velocities [Angstrom/fs]
    Tensor fx, fy, fz;           // Forces [eV/Angstrom]
    Tensor mass;                 // Masses [amu]
    Tensor type_id;              // Atom types (kIntDtype, no grad)

    // Topology (kIntDtype, no grad)
    Tensor bond_i, bond_j, bond_type;
    Tensor angle_i, angle_j, angle_k, angle_type;
    Tensor dihedral_i, dihedral_j, dihedral_k, dihedral_l, dihedral_type;

    // Pair list (kIntDtype indices, kRealDtype displacements)
    Tensor pair_i, pair_j;
    Tensor pair_dx, pair_dy, pair_dz;
    Tensor pair_dist;

    // System scalars
    Tensor energy;               // [1], potential energy [eV]
    Tensor kinetic;              // [1], kinetic energy [eV]
    Tensor stress;               // [6], Voigt [eV/Angstrom^3]
    Tensor box;                  // [9], row-major Mat3 [Angstrom]

    // Host-cached metadata
    Int n_atoms() const;
    Int n_bonds() const;
    Int n_pairs() const;
};
} // namespace ATV::torch
```

## Appendix B: `torch::Context` Complete Reference

```cpp
namespace ATV::torch {
struct Context {
    bool grad_enabled = true;  // Controls autograd graph construction
};
} // namespace ATV::torch
```

## Appendix C: StreamGuard

When a torch operator runs inside a CUDA backend simulation, the simulation
has its own CUDA stream (stored as `void*` in `cuda::Context`). LibTorch uses
a separate default stream. The `StreamGuard` class bridges the two:

```cpp
namespace ATV::torch {
class StreamGuard {
public:
    explicit StreamGuard(void* atv_stream);
    ~StreamGuard();
};
} // namespace ATV::torch
```

Construction captures the current Torch stream, then sets Torch's current
stream to the one passed in. Destruction restores the original. This is RAII:
even if an exception occurs, the original stream is restored.

```cpp
// Inside a cross-device component's launch():
void launch(const LaunchContext& ctx) override {
    ATV::torch::StreamGuard guard(ctx.runtime.stream);
    auto tframe = ATV::torch::to_torch(ctx.frame());
    auto energy = ATV::torch::lj(tframe, eps, sig, cutoff);
    // All Torch ops run on the simulation's stream
}
// guard destroyed here -> previous stream restored
```

## Appendix D: Directory Layout

```
include/atomiverse/torch/
    frame.cuh                  torch::Frame definition
    context.cuh                torch::Context definition
    transfer.cuh               to_torch, to_cuda, to_cpu, from_cpu, destroy
    stream_guard.cuh           StreamGuard RAII
    md/
        lj.cuh                 Differentiable LJ operator
        neighbor_list.cuh      NeighborList (non-differentiable)
    bonded/
        harmonic_bond.cuh      Differentiable HarmonicBond operator

src/torch/
    transfer.cu                Transfer function implementations
    stream_guard.cu            StreamGuard implementation
    md/
        lj.cu                  LJ forward (kernel reuse) + backward + autograd
        neighbor_list.cu       NeighborList (delegates to cuda component)
    bonded/
        harmonic_bond.cu       HarmonicBond forward + backward + autograd
```

This mirrors the conventions of the other backends (`cuda/md/`, `cuda/bonded/`,
etc.) with `torch/` as the top-level backend directory and functionality
organized underneath.

# User Guide

Read this manual if your goal is to **run** Atomiverse — build it, run a
calculation, modify an example — without first learning the internal
architecture.

It assumes no computational-chemistry background and only basic modern C++:
structs, loops, function calls. Every example in the repository assumes
**C++17 or later**.

## The one mental model

Everything in Atomiverse is the same four steps:

1. Fill a [`cpu::Frame`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/frame/frame.h)
   with geometry, masses, atom types, and optional topology.
2. Pick a **driver** — the kind of calculation you want (single point, dynamics,
   optimization, vibrations, sampling).
3. Add **components** — the physics, in execution order.
4. Run it, then read `frame`, the output files, or the result struct.

The driver owns the loop; the components supply the physics. You change *method*
by swapping a component, never by switching to a different program.

## Read in this order

1. **[Getting Started](getting-started.md)** — from a fresh checkout to a
   converged calculation, including the HPC-module gotchas that trip people up.
2. **[Methods](methods.md)** — what the engine can compute today, and which
   driver plus which component runs each one.
3. **[SCF Learning Path](scf-roadmap.md)** — how the electronic-structure
   chapters fit together, if quantum chemistry is what brought you here.
4. **[Tutorial](../tutorial/index.md)** — the long-form book; each chapter pairs
   with a runnable program under `examples/`.
5. **[Switching Backends](backend-switching.md)** — moving the same problem
   between CPU and CUDA.
6. **[LibTorch Backend](libtorch-backend.md)** — only if you need differentiable
   operators or PyTorch interop.

## When to switch manuals

Use the **[Developer Manual](../developer/index.md)** when you want to know how
`Driver`, `Plan`, `Component`, and `Frame` fit together, or when you need to add
a feature. Use the **[API Reference](../api/index.md)** when you already know the
class you need and want the exact signature.

# Developer Manual

The developer manual explains how Atomiverse is put together internally. Read
this part of the documentation when you want to change the codebase rather than
just use it.

This manual still assumes the reader may be new to computational chemistry.
Where chemistry terms appear, the goal is to explain why the software is
organized the way it is, not to assume prior domain knowledge. Code examples
and contributor guidance use **C++17 or later** throughout.

## What To Read

Start with [Architecture](architecture.md). That chapter explains the main
design choice of the repository: a unified facade API drives two different
execution backends, while the user continues to work with a host-side
`cpu::Frame`.

Then read [Extending Components](extending-components.md). That page explains
how to add a new component, how the lifecycle rules work, and where tests and
examples should go.

For CPU-specific work, two further pages drill into the details:
[CPU Backend Internals](cpu-internals.md) explains *why* the CPU side
is split into three driver categories, why role and capability are
orthogonal axes, why SCF is an inner solver rather than a driver
subclass, and how the iteration-status channel signals convergence.
[Extending CPU Subsystems](extending-cpu-subsystems.md) is the
companion handbook with three case studies — adding an MD integrator,
adding a minimizer, and adding an SCF inner-step helper — at the level
of "what file do I create, what do I inherit, what does the test look
like?"

For molecular dynamics specifically, [MD Execution Model](md-execution-model.md)
is the Atomiverse counterpart of the LAMMPS *Developer flow* page: it walks one
timestep as a fixed sequence of named **phases**, documents what each phase is
responsible for (even when empty), and explains the `(phase, order)` binding
model a component uses to attach itself to the step. Read it before adding or
moving any MD component.

## Codebase Map

The repository is easier to navigate if you group it by responsibility instead
of by file extension.

| Directory | Role |
|---|---|
| `include/atomiverse/core` | Small shared vocabulary types such as `Real`, `Int`, `Mat3`, `Scheduling`, and device tags |
| `include/atomiverse/facade` | Backend-agnostic orchestration types such as `Driver`, `Plan`, and phase enums |
| `include/atomiverse/cpu`, `src/cpu` | CPU data structures, CPU backend, CPU SCF solvers, CPU-side analysis and output |
| `include/atomiverse/cuda`, `src/cuda` | CUDA backend, staged components, transfers, and async observer bridges |
| `examples` | End-to-end programs that show supported workflows |
| `tests` | Regression tests for planning, backends, MD, SCF, RHF, UHF, UKS, and integrals |
| `docs` | User guide, developer manual, API reference entry page, and tutorials |

## Development Mindset

If you are deciding where a change belongs, ask these questions in order.

1. Is this a shared concept that both backends need to talk about? It probably belongs in `core` or `facade`.
2. Is this algorithm CPU-only, especially xtensor-based linear algebra or post-processing? It probably belongs under `cpu`.
3. Is this staged GPU execution, device memory management, or async host transfer? It probably belongs under `cuda`.
4. Does this change need a runnable example or a regression test? In Atomiverse, the answer is usually yes.

## Maintaining The Docs

The documentation now has three source-level entry points:

| Path | Purpose |
|---|---|
| `docs/user-guide` | task-oriented pages for users |
| `docs/developer` | architecture and contributor guidance |
| `docs/api` | landing page into generated API reference |

The generated reference site is assembled by Doxygen from `include/` and
`docs/`. The main landing page lives in `docs/mainpage.dox`, and the navigation
layout is configured in `docs/DoxygenLayout.xml`.

If you add or rename major documentation sections, update all three of these:

1. the markdown page under `docs/...`
2. the top-level landing page in `docs/mainpage.dox`
3. the source documentation pointers in `README.md`

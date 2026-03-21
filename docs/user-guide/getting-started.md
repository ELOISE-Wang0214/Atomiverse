# Getting Started

Atomiverse is a C++/CUDA molecular simulation codebase. The same repository
contains three kinds of workflows: molecular dynamics, minimization, and
electronic-structure calculations. You do not need to learn all three at once.
The easiest way to start is to build one example, run it, and then read the
matching tutorial.

All build recipes and examples in this repository assume **C++17 or later**.
The CUDA side is also compiled in C++17 mode.

## What You Need

At minimum, you need CMake and a C++17-capable compiler. Some workflows also
need BLAS and LAPACK because the CPU reference solvers use dense linear algebra.
CUDA workflows additionally need the CUDA toolkit.

Important CMake options:

- `ATV_BUILD_CUDA`: build the CUDA backend
- `ATV_BUILD_RHF`: enable RHF, UHF, and UKS electronic-structure code
- `ATV_BUILD_EXAMPLES`: build the example programs
- `ATV_BUILD_TESTS`: build the test suite
- `ATV_USE_MATPLOT`: enable plot generation in selected examples

## Recommended Build Profiles

If you want the smallest CPU-only starting point, build the toy SCF example:

```bash
cmake -B build -DATV_BUILD_CUDA=OFF -DATV_BUILD_EXAMPLES=ON
cmake --build build --target toy-scf
```

If you want the CPU Hartree-Fock examples as well:

```bash
cmake -B build -DATV_BUILD_CUDA=OFF -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
cmake --build build --target h2-rhf
```

If you want the CUDA workflows:

```bash
cmake -B build -DATV_BUILD_CUDA=ON -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
cmake --build build
```

If you only want the documentation:

```bash
cmake -B build/docs -DATV_DOCS_ONLY=ON
cmake --build build/docs --target docs
```

## Which Example Should You Run First?

Choose the example that matches the question you are trying to answer.

| Program | When to use it | What you learn |
|---|---|---|
| `toy-scf.cpp` | You are completely new to SCF | The feedback loop without chemistry details |
| `scf.cu` | You want to see the staged CUDA driver model | How SCF is expressed as `Driver + Component` |
| `h2-rhf.cpp` | You want the first real chemistry result | How integrals and the CPU RHF solver fit together |
| `h2_rhf.cu` | You want CPU and CUDA RHF in one program | How the same `cpu::Frame` feeds both backends |
| `lj.cu` | You want a full molecular simulation workflow | Minimization, MD, thermo output, trajectory, and RDF |

## The One Data Structure To Remember

Most beginner confusion comes from trying to track too many types at once. Start
with one rule instead: the user-facing system description lives in
[`cpu::Frame`](../../include/atomiverse/cpu/frame.h).

You fill that frame with atom positions, masses, atom types, and optionally
topology such as bonds or angles. CPU workflows operate on it directly. CUDA
workflows copy it to a device-side frame internally and copy results back when
the run finishes.

That means you do not have to manually manage GPU memory just to run a standard
simulation.

## A Good Reading Order For New Users

Read the documentation in this order if you are starting from zero:

1. This page
2. [SCF Learning Path](scf-roadmap.md)
3. [`docs/tutorial/01-scf.md`](../tutorial/01-scf.md)
4. [`docs/tutorial/02-rhf.md`](../tutorial/02-rhf.md)
5. [`docs/tutorial/03-uhf.md`](../tutorial/03-uhf.md)
6. [`docs/tutorial/04-uks.md`](../tutorial/04-uks.md)

The four tutorial chapters are intentionally progressive. The first chapter
teaches the SCF loop as a matrix feedback problem. The later chapters keep the
same loop and add physically meaningful ingredients one at a time.

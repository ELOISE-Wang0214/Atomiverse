# User Guide

The user guide is the part of the documentation you should read if your goal is
to run Atomiverse, modify an example, or understand what the software is doing
without first learning the full internal architecture.

This guide assumes two things. First, you may have no computational chemistry
background. Second, you may know only basic modern C++, at roughly the level of
reading simple structs, loops, and function calls. All examples in this
repository assume **C++17 or later**.

Atomiverse has one central idea that makes the rest of the code easier to read:
you always describe the system in a host-side [`cpu::Frame`](../../include/atomiverse/cpu/frame.h),
and then you choose a driver that runs a staged workflow on either the CPU or
the CUDA backend.

## Start Here

If you are new to the project, read these pages in order:

1. [Getting Started](getting-started.md)
2. [SCF Learning Path](scf-roadmap.md)
3. [SCF Deep-Dive Tutorials](../tutorial/index.md)

The getting-started page explains what to build and which example to run first.
The roadmap explains how the SCF chapters fit together. The tutorial track then
teaches the self-consistent field workflow from the toy model upward to RHF,
UHF, and UKS.

## What Atomiverse Can Do Today

Atomiverse currently exposes three major workflow families.

- Molecular dynamics through staged CUDA drivers such as `md::MD`
- Energy minimization through staged CUDA drivers such as `minimizer::SD`
- Electronic-structure workflows through CPU solvers and CUDA SCF components

If you only need one mental model, use this one:

1. Fill a `cpu::Frame` with geometry, masses, atom types, and optional topology.
2. Choose a driver or a CPU solver.
3. Add components in execution order.
4. Run the workflow and inspect `frame`, output files, or result structs.

## When To Switch To Other Manuals

Use the [Developer Manual](../developer/index.md) when you want to understand how
`Driver`, `Plan`, `Backend`, and `Component` fit together, or when you need to
add a new feature.

Use the [API Reference](../api/index.md) when you already know which class or
function you need and want the exact header, signature, and member list.

# API Reference

This page is the front door to the generated API reference. It is not meant to
duplicate the class and file pages that Doxygen already generates. Its job is to
help you find the right part of the reference quickly.

All public code examples in Atomiverse assume **C++17 or later**.

## How To Use The Generated Reference

Use the **Classes**, **Namespaces**, and **Files** tabs when you already know the
name of the type you want. Those pages are the most authoritative source for
signatures, inheritance, public members, and documented headers.

Use this page when you know the task you want to perform, but not yet the exact
class or function name.

## Entry Points By Task

### Describe a system

Start with these types:

- [`cpu::Frame`](../../include/atomiverse/cpu/frame/frame.h)
- [`Mat3`](../../include/atomiverse/core/mat3.h)
- [`Real`](../../include/atomiverse/core/real.h)
- [`Int`](../../include/atomiverse/core/int.h)

### Build staged workflows

Start with these facade headers:

- [`Driver`](../../include/atomiverse/facade/driver.h)
- [`Plan`](../../include/atomiverse/facade/plan.h)
- [`Phase` enums](../../include/atomiverse/facade/stage.h)

### Run molecular dynamics and minimization

Start with these CPU driver and component headers (the ones the tutorials
use; each has a CUDA counterpart under `cuda/` with the same name):

- [`cpu::md::MD`](../../include/atomiverse/cpu/md/md.h)
- [`cpu::md::VelocityVerlet`](../../include/atomiverse/cpu/md/velocity_verlet.h)
- [`cpu::md::NeighborList`](../../include/atomiverse/cpu/neighbor/neighbor_list.h)
- [`cpu::PairLJCut`](../../include/atomiverse/cpu/potentials/pair_lj_cut.h)
- [`cpu::minimizer::Minimize`](../../include/atomiverse/cpu/minimizer/minimize.h) + `LBFGS` / `CG`

The CUDA equivalents are `cuda::md::MD` (`cuda/md/md.cuh`),
`cuda::PairLJCut`, `cuda::minimizer::Minimize`, etc. — same API, GPU-resident.

### Run electronic-structure calculations on the CPU

Start with these CPU headers:

- [`IntegralPack` and `build_integrals`](../../include/atomiverse/cpu/scf/integrals.h)
- [`rhf_solve`](../../include/atomiverse/cpu/scf/rhf_solve.h)
- [`uhf_solve`](../../include/atomiverse/cpu/scf/uhf_solve.h)
- [`uks_solve`](../../include/atomiverse/cpu/scf/uks_solve.h)

### Run SCF workflows on the CUDA backend

CUDA SCF is an inner solver (`cuda::scf::SCFSolver`) wrapped by an evaluator
component, not a staged driver. Start with these headers:

- [`cuda::scf::SCFEvaluator`](../../include/atomiverse/cuda/scf/scf_evaluator.cuh)
- [`cuda::scf::SCFSolver`](../../include/atomiverse/cuda/scf/scf_solver.cuh)
- [`cuda::Evaluate`](../../include/atomiverse/cuda/driver/evaluate.cuh)
- [`cuda::single_point`](../../include/atomiverse/cuda/driver/single_point.cuh)

### Record output and analyze results

Start with these types:

- [`cpu::Thermo`](../../include/atomiverse/cpu/recorder/thermo.h)
- [`trajectory_append` (plain-XYZ helper)](../../include/atomiverse/cpu/recorder/trajectory.h)
- [`ATV::RDF`](../../include/atomiverse/cpu/observer/rdf.h)
- [`MolRec`](../../include/atomiverse/cpu/io/molrec.h)
- [`cpu::XYZ`](../../include/atomiverse/cpu/recorder/xyz.h)
- [`mulliken_charge` and related helpers](../../include/atomiverse/cpu/observer/mulliken.h)
- [`DensityGrid`](../../include/atomiverse/cpu/observer/density.h)

## When To Leave The API Reference

If you need a conceptual explanation of what the software does, go back to the
[User Guide](../user-guide/index.md).

If you need to understand why the code is structured this way or where to place
a new feature, switch to the [Developer Manual](../developer/index.md).

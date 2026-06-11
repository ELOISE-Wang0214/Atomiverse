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

- [`cpu::Frame`](../../include/atomiverse/cpu/frame.h)
- [`Mat3`](../../include/atomiverse/core/mat3.h)
- [`Real`](../../include/atomiverse/core/real.h)
- [`Int`](../../include/atomiverse/core/int.h)

### Build staged workflows

Start with these facade headers:

- [`Driver`](../../include/atomiverse/facade/driver.h)
- [`Plan`](../../include/atomiverse/facade/plan.h)
- [`Phase` enums](../../include/atomiverse/facade/stage.h)

### Run molecular dynamics and minimization

Start with these CUDA-side driver and component headers:

- [`md::MD`](../../include/atomiverse/cuda/md/md.cuh)
- [`VelocityVerlet`](../../include/atomiverse/cuda/md/velocity_verlet.cuh)
- [`NeighborList`](../../include/atomiverse/cuda/neighbor/neighbor_list.cuh)
- [`PairLJCut`](../../include/atomiverse/cuda/potentials/pair_lj_cut.cuh)
- [`minimizer::SD`](../../include/atomiverse/cuda/minimizer/sd.cuh)

### Run electronic-structure calculations on the CPU

Start with these CPU headers:

- [`IntegralPack` and `build_integrals_sto3g`](../../include/atomiverse/cpu/scf/integrals.h)
- [`rhf_solve`](../../include/atomiverse/cpu/scf/rhf_solve.h)
- [`uhf_solve`](../../include/atomiverse/cpu/scf/uhf_solve.h)
- [`uks_solve`](../../include/atomiverse/cpu/scf/uks_solve.h)

### Run staged SCF workflows on the CUDA backend

Start with these headers:

- [`scf::SCF`](../../include/atomiverse/cuda/scf/scf.cuh)
- [`Hamiltonian`](../../include/atomiverse/cuda/scf/hamiltonian.cuh)
- [`EigenSolver`](../../include/atomiverse/cuda/scf/eigen_solver.cuh)
- [`DensityBuilder`](../../include/atomiverse/cuda/scf/density_builder.cuh)
- [`RHFFock`](../../include/atomiverse/cuda/scf/rhf_fock.cuh)
- [`RHFDensity`](../../include/atomiverse/cuda/scf/rhf_density.cuh)
- [`ConvergenceCheck`](../../include/atomiverse/cuda/scf/convergence_check.cuh)

### Record output and analyze results

Start with these types:

- [`Thermo`](../../include/atomiverse/cuda/recorder/thermo.cuh)
- [`Recorder`](../../include/atomiverse/cuda/recorder/recorder.cuh)
- [`RDF`](../../include/atomiverse/cuda/observer/rdf.cuh)
- [`MolRec`](../../include/atomiverse/cpu/molrec.h)
- [`XYZ`](../../include/atomiverse/cpu/recorder/xyz.h)
- [`mulliken_charge` and related helpers](../../include/atomiverse/cpu/observer/mulliken.h)
- [`DensityGrid`](../../include/atomiverse/cpu/observer/density.h)

## When To Leave The API Reference

If you need a conceptual explanation of what the software does, go back to the
[User Guide](../user-guide/index.md).

If you need to understand why the code is structured this way or where to place
a new feature, switch to the [Developer Manual](../developer/index.md).

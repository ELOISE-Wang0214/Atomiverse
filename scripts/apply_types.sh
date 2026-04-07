#!/bin/bash
# Mechanical refactoring: Replace Real/Int with CUDA mixed-precision types
# in all CUDA backend files EXCEPT frame.cuh and types.cuh.
set -euo pipefail

ROOT="/Users/roykid/work/molcrafts/Atomiverse"

##############################################################################
# PART 1: PUBLIC HEADERS (include/atomiverse/cuda/**/*.cuh)
# These are CUDA-free API headers: Real -> double, Int -> int for members/params
##############################################################################

# List of public headers to process (exclude frame.cuh, types.cuh)
PUBLIC_HEADERS=(
  "include/atomiverse/cuda/neighbor_list.cuh"
  "include/atomiverse/cuda/potentials/lj.cuh"
  "include/atomiverse/cuda/potentials/harmonic_bond.cuh"
  "include/atomiverse/cuda/potentials/harmonic_angle.cuh"
  "include/atomiverse/cuda/potentials/multi_harmonic_dihedral.cuh"
  "include/atomiverse/cuda/md/velocity_verlet.cuh"
  "include/atomiverse/cuda/md/init.cuh"
  "include/atomiverse/cuda/md/kinetic_energy.cuh"
  "include/atomiverse/cuda/md/nose_hoover.cuh"
  "include/atomiverse/cuda/md/barostat.cuh"
  "include/atomiverse/cuda/md/deform.cuh"
  "include/atomiverse/cuda/md/md.cuh"
  "include/atomiverse/cuda/minimizer/sd.cuh"
  "include/atomiverse/cuda/minimizer/cg.cuh"
  "include/atomiverse/cuda/observer/stress.cuh"
  "include/atomiverse/cuda/observer/rdf.cuh"
  "include/atomiverse/cuda/observer/mulliken.cuh"
  "include/atomiverse/cuda/recorder/thermo.cuh"
  "include/atomiverse/cuda/recorder/recorder.cuh"
  "include/atomiverse/cuda/recorder/density_recorder.cuh"
  "include/atomiverse/cuda/scf/hamiltonian.cuh"
  "include/atomiverse/cuda/scf/eigen_solver.cuh"
  "include/atomiverse/cuda/scf/density_builder.cuh"
  "include/atomiverse/cuda/scf/convergence_check.cuh"
  "include/atomiverse/cuda/scf/rhf_fock.cuh"
  "include/atomiverse/cuda/scf/rhf_density.cuh"
  "include/atomiverse/cuda/scf/scf.cuh"
  "include/atomiverse/cuda/scf/device_integrals.cuh"
  "include/atomiverse/cuda/component.cuh"
  "include/atomiverse/cuda/context.cuh"
  "include/atomiverse/cuda/transfer.cuh"
  "include/atomiverse/cuda/backend.cuh"
  "include/atomiverse/cuda/backend_traits.cuh"
  "include/atomiverse/cuda/error.cuh"
  "include/atomiverse/cuda/launch.cuh"
  "include/atomiverse/cuda/timer.cuh"
  "include/atomiverse/cuda/async_bridge.cuh"
)

# Kernel declaration headers need cuda types, not just double
KERNEL_HEADERS=(
  "include/atomiverse/cuda/potentials/lj_kernel.cuh"
  "include/atomiverse/cuda/potentials/harmonic_bond_kernel.cuh"
)

##############################################################################
# PART 2: SOURCE FILES (src/cuda/**/*.cu and src/cuda/**/*.cuh)
##############################################################################

SOURCE_FILES=(
  "src/cuda/transfer.cu"
  "src/cuda/driver.cu"
  "src/cuda/async_bridge.cu"
  "src/cuda/neighbor_list.cu"
  "src/cuda/potentials/lj.cu"
  "src/cuda/potentials/harmonic_bond.cu"
  "src/cuda/potentials/harmonic_angle.cu"
  "src/cuda/potentials/multi_harmonic_dihedral.cu"
  "src/cuda/md/velocity_verlet.cu"
  "src/cuda/md/init.cu"
  "src/cuda/md/kinetic_energy.cu"
  "src/cuda/md/nose_hoover.cu"
  "src/cuda/md/deform.cu"
  "src/cuda/md/barostat.cu"
  "src/cuda/minimizer/sd.cu"
  "src/cuda/observer/stress.cu"
  "src/cuda/scf/hamiltonian.cu"
  "src/cuda/scf/eigen_solver.cu"
  "src/cuda/scf/density_builder.cu"
  "src/cuda/scf/convergence_check.cu"
  "src/cuda/scf/driver.cu"
  "src/cuda/scf/rhf_fock.cu"
  "src/cuda/scf/rhf_density.cu"
  "src/cuda/scf/device_integrals.cu"
  "src/cuda/md/neighborlist/sort_and_cluster.cu"
  "src/cuda/md/neighborlist/setup.cu"
  "src/cuda/md/neighborlist/grid.cu"
  "src/cuda/md/neighborlist/pair_search.cu"
  "src/cuda/md/neighborlist/fused_search.cu"
  "src/cuda/md/neighborlist/export_pairs.cu"
)

SOURCE_HEADERS=(
  "src/cuda/md/neighborlist/nlist_data.cuh"
  "src/cuda/md/neighborlist/nlist_impl.cuh"
  "src/cuda/md/neighborlist/nlist_error.cuh"
  "src/cuda/md/neighborlist/nlist_memory.cuh"
  "src/cuda/md/neighborlist/flat_pair_list.cuh"
  "src/cuda/md/neighborlist/morton.cuh"
  "src/cuda/md/neighborlist/triclinic.cuh"
  "src/cuda/scf/cublas_dispatch.cuh"
)

echo "=== Processing public headers ==="
for f in "${PUBLIC_HEADERS[@]}"; do
  fp="$ROOT/$f"
  [ -f "$fp" ] || { echo "SKIP (missing): $f"; continue; }
  echo "  $f"
done

echo "=== Processing kernel headers ==="
for f in "${KERNEL_HEADERS[@]}"; do
  fp="$ROOT/$f"
  [ -f "$fp" ] || { echo "SKIP (missing): $f"; continue; }
  echo "  $f"
done

echo "=== Processing source files ==="
for f in "${SOURCE_FILES[@]}"; do
  fp="$ROOT/$f"
  [ -f "$fp" ] || { echo "SKIP (missing): $f"; continue; }
  echo "  $f"
done

echo "=== Processing source headers ==="
for f in "${SOURCE_HEADERS[@]}"; do
  fp="$ROOT/$f"
  [ -f "$fp" ] || { echo "SKIP (missing): $f"; continue; }
  echo "  $f"
done

echo "All files exist. Ready for transformation."

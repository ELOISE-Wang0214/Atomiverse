#!/usr/bin/env python3
"""
Transform CUDA backend files to use mixed-precision types.

Categories:
1. Public API headers (CUDA-free): Real -> double, Int -> int
2. Kernel declaration headers: Real -> cuda types
3. Source .cu files: Real -> work_t/real_t/prec_t depending on context
4. Internal .cuh in src/: template-based, mostly already parameterized

Strategy: Per-file manual replacement with context awareness.
"""

import os
import re
import sys

ROOT = "/Users/roykid/work/molcrafts/Atomiverse"

def read_file(path):
    with open(os.path.join(ROOT, path), 'r') as f:
        return f.read()

def write_file(path, content):
    with open(os.path.join(ROOT, path), 'w') as f:
        f.write(content)

def transform_public_header(path):
    """Transform CUDA-free API headers: Real -> double, Int -> int for members."""
    content = read_file(path)
    original = content

    # Replace include of core/real.h -- these headers don't need cuda types
    # They just use double for their own members
    content = content.replace('#include "atomiverse/core/real.h"', '')
    content = content.replace('#include "atomiverse/core/int.h"', '')

    # For headers that include component.cuh (which includes frame.cuh which
    # still includes core/real.h), the Real type is still available transitively.
    # But we want to replace Real usage in the header's own API with double.

    # Replace Real in function signatures, member declarations, and return types
    # Be careful not to replace inside strings, comments referencing Real type name

    # Pattern: standalone Real (not part of a longer word)
    # But only in declarations, not in doxygen comments describing the type
    content = re.sub(r'\bReal\b', 'double', content)
    content = re.sub(r'\bInt\b(?!\*)', 'int', content)

    # Fix over-replacements in comments/strings - restore Real in doxygen type references
    # Actually, since these are mechanical, let's just do them all.
    # The comment references to Real as a type are fine as double too.

    if content != original:
        write_file(path, content)
        print(f"  CHANGED: {path}")
    else:
        print(f"  UNCHANGED: {path}")

def transform_kernel_header(path):
    """Transform kernel declaration headers (lj_kernel.cuh, harmonic_bond_kernel.cuh)."""
    content = read_file(path)
    original = content

    # Replace includes
    content = content.replace('#include "atomiverse/core/real.h"',
                              '#include "atomiverse/cuda/types.cuh"')
    content = content.replace('#include "atomiverse/core/int.h"', '')

    # Add using declarations after namespace opening
    # These headers declare __global__ kernels so they need the cuda types
    # The kernel pointers need to be real_t/prec_t/work_t

    # For lj_kernel.cuh:
    if 'lj_kernel' in path:
        # kMinPairDist: constexpr work_t (computation constant)
        content = content.replace('inline constexpr Real kMinPairDist = Real(1e-10);',
                                  'inline constexpr cuda::work_t kMinPairDist = cuda::work_t(1e-10);')

        # lj_force kernel signature: fx,fy,fz are real_t (per-atom force),
        # energy is prec_t, pair data is real_t, scalar params are work_t
        content = content.replace(
            '''__global__ void lj_force(
    Real* __restrict__ fx, Real* __restrict__ fy, Real* __restrict__ fz,
    Real* __restrict__ energy,
    const Int* __restrict__ pi, const Int* __restrict__ pj,
    const Real* __restrict__ pdx, const Real* __restrict__ pdy,
    const Real* __restrict__ pdz, const Real* __restrict__ pdist,
    Real epsilon, Real sigma, Real cutoff, Real v_shift,
    const Int* __restrict__ n_pairs_ptr);''',
            '''__global__ void lj_force(
    cuda::real_t* __restrict__ fx, cuda::real_t* __restrict__ fy, cuda::real_t* __restrict__ fz,
    cuda::prec_t* __restrict__ energy,
    const int* __restrict__ pi, const int* __restrict__ pj,
    const cuda::real_t* __restrict__ pdx, const cuda::real_t* __restrict__ pdy,
    const cuda::real_t* __restrict__ pdz, const cuda::real_t* __restrict__ pdist,
    cuda::work_t epsilon, cuda::work_t sigma, cuda::work_t cutoff, cuda::work_t v_shift,
    const int* __restrict__ n_pairs_ptr);''')

        # lj_shift inline function
        content = content.replace(
            '''inline Real lj_shift(Real epsilon, Real sigma, Real cutoff) {
    Real sr = sigma / cutoff;
    Real sr2 = sr * sr;
    Real sr6 = sr2 * sr2 * sr2;
    Real sr12 = sr6 * sr6;
    return Real(4) * epsilon * (sr12 - sr6);
}''',
            '''inline double lj_shift(double epsilon, double sigma, double cutoff) {
    double sr = sigma / cutoff;
    double sr2 = sr * sr;
    double sr6 = sr2 * sr2 * sr2;
    double sr12 = sr6 * sr6;
    return 4.0 * epsilon * (sr12 - sr6);
}''')

    elif 'harmonic_bond_kernel' in path:
        # harmonic_bond_k kernel: positions are real_t, forces are real_t,
        # energy is prec_t, K/r0 are work_t (per-type params), box is prec_t
        content = content.replace(
            '''__global__ void harmonic_bond_k(
    const Real* __restrict__ x, const Real* __restrict__ y, const Real* __restrict__ z,
    Real* __restrict__ fx, Real* __restrict__ fy, Real* __restrict__ fz,
    Real* __restrict__ energy,
    const Int* __restrict__ bi, const Int* __restrict__ bj, const Int* __restrict__ bt,
    const Real* __restrict__ K, const Real* __restrict__ r0,
    const Real* __restrict__ box, Int n_bonds);''',
            '''__global__ void harmonic_bond_k(
    const cuda::real_t* __restrict__ x, const cuda::real_t* __restrict__ y, const cuda::real_t* __restrict__ z,
    cuda::real_t* __restrict__ fx, cuda::real_t* __restrict__ fy, cuda::real_t* __restrict__ fz,
    cuda::prec_t* __restrict__ energy,
    const int* __restrict__ bi, const int* __restrict__ bj, const int* __restrict__ bt,
    const cuda::work_t* __restrict__ K, const cuda::work_t* __restrict__ r0,
    const cuda::prec_t* __restrict__ box, int n_bonds);''')

    if content != original:
        write_file(path, content)
        print(f"  CHANGED: {path}")
    else:
        print(f"  UNCHANGED: {path}")

def transform_source_file(path):
    """Transform .cu source files with mixed-precision types."""
    content = read_file(path)
    original = content

    basename = os.path.basename(path)

    # Remove 'using Real = ATV::Real;' and 'using Int = ATV::Int;'
    content = content.replace('using Real = ATV::Real;\n', '')
    content = content.replace('using Int = ATV::Int;\n', '')

    # Add using declarations for cuda types after includes, before namespace
    # We need: using real_t = cuda::real_t; etc. within the ATV namespace
    # Actually, since these are .cu files inside namespace ATV, we can use
    # cuda::real_t directly, or add using declarations.

    # For .cu files in src/cuda/ that are inside namespace ATV:
    # The approach is to add 'using namespace cuda;' style declarations
    # or use the qualified names. Let me use qualified names in kernel code
    # and add type aliases at file scope.

    # The cleanest approach: add a block of using declarations after the includes,
    # before any namespace. Then replace Real/Int throughout.

    # First, let's handle the specific files that need special treatment.

    if content != original:
        write_file(path, content)
        print(f"  CHANGED: {path}")
    else:
        print(f"  UNCHANGED: {path}")

# Process files
print("=== Public headers ===")
headers = [
    "include/atomiverse/cuda/potentials/harmonic_bond.cuh",
    "include/atomiverse/cuda/potentials/harmonic_angle.cuh",
    "include/atomiverse/cuda/potentials/multi_harmonic_dihedral.cuh",
    "include/atomiverse/cuda/md/velocity_verlet.cuh",
    "include/atomiverse/cuda/md/init.cuh",
    "include/atomiverse/cuda/md/kinetic_energy.cuh",
    "include/atomiverse/cuda/md/nose_hoover.cuh",
    "include/atomiverse/cuda/md/barostat.cuh",
    "include/atomiverse/cuda/md/deform.cuh",
    "include/atomiverse/cuda/minimizer/sd.cuh",
    "include/atomiverse/cuda/scf/hamiltonian.cuh",
    "include/atomiverse/cuda/scf/eigen_solver.cuh",
    "include/atomiverse/cuda/scf/convergence_check.cuh",
    "include/atomiverse/cuda/scf/rhf_fock.cuh",
    "include/atomiverse/cuda/scf/rhf_density.cuh",
    "include/atomiverse/cuda/scf/device_integrals.cuh",
    "include/atomiverse/cuda/observer/rdf.cuh",
    "include/atomiverse/cuda/observer/mulliken.cuh",
    "include/atomiverse/cuda/recorder/thermo.cuh",
    "include/atomiverse/cuda/recorder/recorder.cuh",
    "include/atomiverse/cuda/recorder/density_recorder.cuh",
]

for h in headers:
    transform_public_header(h)

print("\n=== Kernel headers ===")
for h in ["include/atomiverse/cuda/potentials/lj_kernel.cuh",
          "include/atomiverse/cuda/potentials/harmonic_bond_kernel.cuh"]:
    transform_kernel_header(h)

print("\nDone. Source files need manual per-file treatment.")

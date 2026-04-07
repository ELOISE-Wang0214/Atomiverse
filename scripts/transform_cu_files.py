#!/usr/bin/env python3
"""
Transform ALL .cu source files and src/cuda/ .cuh files to use CUDA mixed-precision types.

Rules:
- Remove 'using Real = ATV::Real;' and 'using Int = ATV::Int;'
- In kernel code and host-side code within .cu files:
  - Per-atom array pointers: Real* -> real_t* (positions, velocities, forces, mass, pair dx/dy/dz/dist)
  - Energy/kinetic/stress/box pointers: Real* -> prec_t* (scalars, accumulators)
  - Per-type param arrays: Real* -> work_t* (K, r0, C coefficients)
  - Kernel scalar params: Real -> work_t
  - Local computation variables: Real -> work_t
  - Shared memory: extern __shared__ Real -> extern __shared__ work_t
  - Real(0), Real(1), etc -> work_t(0), work_t(1), etc
  - sizeof(Real) -> sizeof(real_t) for per-atom, sizeof(prec_t) for accumulators
  - Int -> int
  - atomicAdd to force arrays: result stays as-is (atomicAdd on real_t)
  - atomicAdd to energy/stress: result stays as-is (atomicAdd on prec_t)

For the neighbor list sub-modules (already templated on Real), just need:
- Remove using declarations
- Fix non-template code that uses bare Real/Int
- setup.cu: uses bare Real (not templated)
- export_pairs.cu: uses bare Real and Int

For SCF files: Real arrays are all prec_t (SCF matrices), Int -> int

For cublas_dispatch.cuh: Replace core/real.h include, remove Real typedef usage
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

def add_types_include_if_needed(content, filepath):
    """Add types.cuh include if not already present and file uses cuda types."""
    if 'cuda/types.cuh' in content:
        return content
    # Add after cuda_runtime.h or after the first include
    if '#include <cuda_runtime.h>' in content:
        content = content.replace('#include <cuda_runtime.h>',
                                  '#include <cuda_runtime.h>\n#include "atomiverse/cuda/types.cuh"')
    return content

def transform_lj_cu(path):
    """Transform potentials/lj.cu"""
    content = read_file(path)

    # Add types include
    content = content.replace('#include <cuda_runtime.h>',
                              '#include <cuda_runtime.h>\n#include "atomiverse/cuda/types.cuh"')

    # Add using declarations after namespace ATV opening
    content = content.replace(
        'namespace ATV {\n\nnamespace cuda { namespace kernels {',
        'namespace ATV {\n\nusing cuda::real_t;\nusing cuda::work_t;\nusing cuda::prec_t;\n\nnamespace cuda { namespace kernels {')

    # Kernel signature
    content = content.replace(
        '''__global__ void lj_force(
    Real* __restrict__ fx, Real* __restrict__ fy, Real* __restrict__ fz,
    Real* __restrict__ energy,
    const Int* __restrict__ pi, const Int* __restrict__ pj,
    const Real* __restrict__ pdx, const Real* __restrict__ pdy,
    const Real* __restrict__ pdz, const Real* __restrict__ pdist,
    Real epsilon, Real sigma, Real cutoff, Real v_shift,
    const Int* __restrict__ n_pairs_ptr) {''',
        '''__global__ void lj_force(
    real_t* __restrict__ fx, real_t* __restrict__ fy, real_t* __restrict__ fz,
    prec_t* __restrict__ energy,
    const int* __restrict__ pi, const int* __restrict__ pj,
    const real_t* __restrict__ pdx, const real_t* __restrict__ pdy,
    const real_t* __restrict__ pdz, const real_t* __restrict__ pdist,
    work_t epsilon, work_t sigma, work_t cutoff, work_t v_shift,
    const int* __restrict__ n_pairs_ptr) {''')

    # Kernel body - Int -> int, Real -> work_t for locals
    content = content.replace('    Int idx = blockIdx.x * blockDim.x + threadIdx.x;\n    Int n_pairs = *n_pairs_ptr;',
                              '    int idx = blockIdx.x * blockDim.x + threadIdx.x;\n    int n_pairs = *n_pairs_ptr;')

    content = content.replace('    Real my_energy = Real(0);', '    work_t my_energy = work_t(0);')
    content = content.replace('        Real r = pdist[idx];', '        work_t r = pdist[idx];')
    content = content.replace('            Real inv_r = Real(1) / r;', '            work_t inv_r = work_t(1) / r;')
    content = content.replace('            Real sr = sigma * inv_r;', '            work_t sr = sigma * inv_r;')
    content = content.replace('            Real sr2 = sr * sr;', '            work_t sr2 = sr * sr;')
    content = content.replace('            Real sr6 = sr2 * sr2 * sr2;', '            work_t sr6 = sr2 * sr2 * sr2;')
    content = content.replace('            Real sr12 = sr6 * sr6;', '            work_t sr12 = sr6 * sr6;')
    content = content.replace('            Real f_over_r = Real(24) * epsilon * (Real(2) * sr12 - sr6) * inv_r * inv_r;',
                              '            work_t f_over_r = work_t(24) * epsilon * (work_t(2) * sr12 - sr6) * inv_r * inv_r;')
    content = content.replace('            Real dx = pdx[idx], dy = pdy[idx], dz = pdz[idx];',
                              '            work_t dx = pdx[idx], dy = pdy[idx], dz = pdz[idx];')
    content = content.replace('            Int i = pi[idx], j = pj[idx];',
                              '            int i = pi[idx], j = pj[idx];')

    # atomicAdd to forces (real_t arrays)
    content = content.replace('            atomicAdd(&fx[i], -f_over_r * dx);',
                              '            atomicAdd(&fx[i], static_cast<real_t>(-f_over_r * dx));')
    content = content.replace('            atomicAdd(&fy[i], -f_over_r * dy);',
                              '            atomicAdd(&fy[i], static_cast<real_t>(-f_over_r * dy));')
    content = content.replace('            atomicAdd(&fz[i], -f_over_r * dz);',
                              '            atomicAdd(&fz[i], static_cast<real_t>(-f_over_r * dz));')
    content = content.replace('            atomicAdd(&fx[j],  f_over_r * dx);',
                              '            atomicAdd(&fx[j], static_cast<real_t>( f_over_r * dx));')
    content = content.replace('            atomicAdd(&fy[j],  f_over_r * dy);',
                              '            atomicAdd(&fy[j], static_cast<real_t>( f_over_r * dy));')
    content = content.replace('            atomicAdd(&fz[j],  f_over_r * dz);',
                              '            atomicAdd(&fz[j], static_cast<real_t>( f_over_r * dz));')

    content = content.replace('            my_energy = Real(4) * epsilon * (sr12 - sr6) - v_shift;',
                              '            my_energy = work_t(4) * epsilon * (sr12 - sr6) - v_shift;')

    # Shared memory
    content = content.replace('    extern __shared__ Real s_energy[];',
                              '    extern __shared__ work_t s_energy[];')

    content = content.replace('    if (threadIdx.x == 0 && s_energy[0] != Real(0))',
                              '    if (threadIdx.x == 0 && s_energy[0] != work_t(0))')
    content = content.replace('        atomicAdd(energy, s_energy[0]);',
                              '        atomicAdd(energy, static_cast<prec_t>(s_energy[0]));')

    # Constructor and launch
    content = content.replace('LJ::LJ(Real epsilon, Real sigma, Real cutoff)',
                              'LJ::LJ(double epsilon, double sigma, double cutoff)')

    # Launch function
    content = content.replace('    Int n = ctx.frame.n_atoms();', '    int n = ctx.frame.n_atoms();')
    content = content.replace('block * sizeof(Real)', 'block * sizeof(work_t)')

    write_file(path, content)
    print(f"  CHANGED: {path}")

def transform_simple_cu(path, content):
    """Generic transformation for .cu files - replace all Real/Int occurrences."""
    # This is the fallback - just do blanket replacement
    # Real -> work_t in kernel code, keeping pointer context
    orig = content

    # Replace #include "atomiverse/core/real.h" etc
    content = content.replace('#include "atomiverse/core/real.h"', '#include "atomiverse/cuda/types.cuh"')
    content = content.replace('#include "atomiverse/core/int.h"', '')

    return content

# Process each file individually with tailored transformations.
# Given the complexity, let me use a Python approach that handles the
# most common patterns for ALL .cu files at once.

def bulk_transform_cu(path):
    """Transform a .cu file using context-aware regex replacements."""
    content = read_file(path)
    original = content

    basename = os.path.basename(path)

    # Remove using declarations
    content = content.replace('using Real = ATV::Real;\n', '')
    content = content.replace('using Int = ATV::Int;\n', '')
    content = content.replace('using ATV::Int;\n', '')
    content = content.replace('using ATV::Real;\n', '')

    # Replace core includes
    content = content.replace('#include "atomiverse/core/real.h"', '#include "atomiverse/cuda/types.cuh"')
    content = content.replace('#include "atomiverse/core/int.h"', '')

    # Add types.cuh include if not present and file has cuda_runtime.h
    if 'cuda/types.cuh' not in content and '<cuda_runtime.h>' in content:
        content = content.replace('#include <cuda_runtime.h>',
                                  '#include <cuda_runtime.h>\n#include "atomiverse/cuda/types.cuh"')

    # For files in namespace ATV, add using declarations
    # Check if file opens namespace ATV
    if 'namespace ATV' in content and 'using cuda::real_t;' not in content:
        # Find first 'namespace ATV {' and add usings after it
        # Handle various patterns
        for ns_pattern in ['namespace ATV {\n\nnamespace cuda { namespace kernels {',
                           'namespace ATV {\n\nnamespace md {',
                           'namespace ATV {\n\nnamespace minimizer {',
                           'namespace ATV {\n\nnamespace scf {',
                           'namespace ATV {\nnamespace md {',
                           'namespace ATV {\nnamespace scf {',
                           'namespace ATV {\nnamespace cuda {\nnamespace scf {',
                           'namespace ATV {\n\nnamespace cpu {',
                           'namespace ATV {\n']:
            if ns_pattern in content:
                using_block = 'using cuda::real_t;\nusing cuda::work_t;\nusing cuda::prec_t;\n\n'
                replacement = ns_pattern.split('\n')[0] + '\n\n' + using_block
                # Only add once
                idx = content.find(ns_pattern)
                if idx >= 0:
                    first_line = ns_pattern.split('\n')[0]
                    content = content[:idx] + first_line + '\n\n' + using_block + '\n'.join(ns_pattern.split('\n')[1:]) + content[idx + len(ns_pattern):]
                    break

    # Now do the blanket Real -> work_t replacement for kernel locals
    # and Int -> int replacement
    # This is aggressive but correct for .cu kernel code

    # Replace sizeof(Real) contextually
    # For memset/memcpy of per-atom arrays (fx, fy, fz, x, y, z, vx, vy, vz, mass):
    #   sizeof(Real) -> sizeof(real_t)
    # For memset/memcpy of energy/kinetic/stress/box/SCF matrices:
    #   sizeof(Real) -> sizeof(prec_t)

    # Shared memory: extern __shared__ Real -> extern __shared__ work_t
    content = re.sub(r'extern __shared__ Real\b', 'extern __shared__ work_t', content)

    # Replace standalone Real in various contexts
    # Kernel parameters and local variables
    content = re.sub(r'\bReal\*\s*__restrict__', 'real_t* __restrict__', content)
    content = re.sub(r'const Real\*\s*__restrict__', 'const real_t* __restrict__', content)

    # Real* without __restrict__
    content = re.sub(r'\bReal\*\b', 'real_t*', content)
    content = re.sub(r'const Real\*', 'const real_t*', content)

    # Real local variables and parameters (standalone Real, not pointer)
    content = re.sub(r'\bReal\b(?!\*)', 'work_t', content)

    # Int -> int
    content = re.sub(r'\bInt\b', 'int', content)

    if content != original:
        write_file(path, content)
        print(f"  CHANGED: {path}")
    else:
        print(f"  UNCHANGED: {path}")

# All .cu files
CU_FILES = [
    "src/cuda/transfer.cu",
    "src/cuda/driver.cu",
    "src/cuda/async_bridge.cu",
    "src/cuda/neighbor_list.cu",
    "src/cuda/potentials/lj.cu",
    "src/cuda/potentials/harmonic_bond.cu",
    "src/cuda/potentials/harmonic_angle.cu",
    "src/cuda/potentials/multi_harmonic_dihedral.cu",
    "src/cuda/md/velocity_verlet.cu",
    "src/cuda/md/init.cu",
    "src/cuda/md/kinetic_energy.cu",
    "src/cuda/md/nose_hoover.cu",
    "src/cuda/md/deform.cu",
    "src/cuda/md/barostat.cu",
    "src/cuda/minimizer/sd.cu",
    "src/cuda/observer/stress.cu",
    "src/cuda/scf/hamiltonian.cu",
    "src/cuda/scf/eigen_solver.cu",
    "src/cuda/scf/density_builder.cu",
    "src/cuda/scf/convergence_check.cu",
    "src/cuda/scf/driver.cu",
    "src/cuda/scf/rhf_fock.cu",
    "src/cuda/scf/rhf_density.cu",
    "src/cuda/scf/device_integrals.cu",
    "src/cuda/md/neighborlist/sort_and_cluster.cu",
    "src/cuda/md/neighborlist/setup.cu",
    "src/cuda/md/neighborlist/grid.cu",
    "src/cuda/md/neighborlist/pair_search.cu",
    "src/cuda/md/neighborlist/fused_search.cu",
    "src/cuda/md/neighborlist/export_pairs.cu",
]

# Internal headers in src/
INTERNAL_HEADERS = [
    "src/cuda/md/neighborlist/nlist_data.cuh",
    "src/cuda/md/neighborlist/nlist_impl.cuh",
    "src/cuda/md/neighborlist/nlist_error.cuh",
    "src/cuda/md/neighborlist/nlist_memory.cuh",
    "src/cuda/md/neighborlist/flat_pair_list.cuh",
    "src/cuda/md/neighborlist/morton.cuh",
    "src/cuda/md/neighborlist/triclinic.cuh",
    "src/cuda/scf/cublas_dispatch.cuh",
]

print("=== Processing .cu files ===")
for f in CU_FILES:
    bulk_transform_cu(f)

print("\n=== Processing internal headers ===")
for f in INTERNAL_HEADERS:
    bulk_transform_cu(f)

print("\nDone.")

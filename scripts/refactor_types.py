#!/usr/bin/env python3
"""
Mechanical refactoring: Replace Real/Int with CUDA mixed-precision types.

Headers (.cuh) that are CUDA-free (component API headers):
  - Constructor params, member scalars, return types: Real -> double
  - Include swap: core/real.h -> cuda/types.cuh
  - Remove core/int.h include
  - Int -> int for member fields (where not already int)

Headers (.cuh) that contain kernel declarations (lj_kernel.cuh, harmonic_bond_kernel.cuh):
  - Kernel pointer params use cuda types (real_t, prec_t)
  - These need cuda/types.cuh

Source files (.cu):
  - Kernel code: Real -> work_t for locals/computation
  - Per-atom pointers: Real* -> real_t*
  - Energy/kinetic/stress/box pointers: Real* -> prec_t*
  - Shared memory: Real -> work_t
  - sizeof(Real) -> sizeof(real_t) or sizeof(prec_t) as appropriate
  - atomicAdd to forces -> real_t cast
  - atomicAdd to energy/stress/kinetic -> prec_t cast
  - Remove using Real = ATV::Real etc.
  - Int -> int
"""

import re
import sys
import os

def process_file(filepath, content):
    """Process a single file and return new content."""
    basename = os.path.basename(filepath)

    # Skip excluded files
    if basename in ('frame.cuh', 'types.cuh'):
        return content

    return content

if __name__ == '__main__':
    print("Use individual edit commands instead")

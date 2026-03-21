---
name: av-debug
description: Troubleshoot build, test, or runtime failures with structured diagnosis. Read-only — never writes code.
argument-hint: "<error message or symptom>"
user-invocable: true
---

# Debug (diagnosis only — never writes code)

Read CLAUDE.md for build commands and architecture context.

Use this skill to **understand** a failure. To **fix** it, use `/av-fix` instead.

## Procedure

### 1. Classify the failure
- **Build failure**: CMake config error, compiler error, linker error
- **Test failure**: assertion failure, crash, timeout, runtime error
- **Runtime failure**: illegal memory access, kernel launch failure, segfault, NaN/inf results

### 2. Gather information
```bash
cmake --build build -j$(nproc) 2>&1 | tail -100
ctest --test-dir build --output-on-failure 2>&1
```

### 3. Diagnosis by type

**Build failure**: Check CMakeLists.txt for missing sources. Check includes (facade headers must not include cuda_runtime.h or xtensor headers). Check Real/Int type usage.

**Test failure**: Read the failing test and the component under test. Check lifecycle compliance (build allocates, launch computes, destroy frees). For CUDA: check grid/block sizing. For CPU: check xtensor shape mismatches.

**CUDA runtime error**: Check kernel launch config (grid > 0, block <= 1024). Check device pointer validity (allocated in build?). Check async operations match stream. Check cudaGetLastError() after kernel launches.

**CPU runtime error**: Check null/dangling pointer access. Check xtensor shape/dimension mismatches. Check out-of-bounds access on SoA arrays.

**NaN/inf results**: Check force computation (division by zero at r=0?). Check unit conversions (kFtm2v). Check time step size. Check initial conditions.

### 4. Report (no edits)
- Root cause
- Fix recommendation (what to change, not the change itself)
- Preventive measure (test to add)

# libxc Integration & Native XC Fix

**Status:** draft
**Priority:** high
**Depends on:** none

## Problem

Atomiverse UKS-LSDA gives systematically lower energies than PySCF by ~0.1 Ha.
UHF (same integrals, no grid) matches PySCF to float32 precision (~1e-7 Ha),
proving the integrals are correct. The error is entirely in the numerical XC integration pipeline.

### Root cause analysis

Three independent bugs compound:

| # | Component | File | Bug | Impact |
|---|-----------|------|-----|--------|
| 1 | Angular grid | `grid.cpp:160-207` | Abandoned real Lebedev-110 tabulation; uses product grid (10 GL × 11 uniform φ). Product rule has wrong weight distribution for spherical harmonics. | ~0.05 Ha on atoms |
| 2 | Radial grid | `grid.cpp:235-254` | Becke mapping `r = R(1+x)/(1-x)` with Chebyshev points is a valid scheme but differs from PySCF's Mura-Knowles / Treutler-Ahlrichs. More importantly, the number of radial points (50) may be too few for this particular mapping. | ~0.02 Ha |
| 3 | Bragg-Slater radii | `grid.cpp:25-36` | Only H/He/Li defined; Z≥4 falls back to 1.0 Å. He uses 0.25 Å (estimated, not standard). | Affects Becke partition weights for multi-atom systems |

The XC functional itself (VWN5 parameters, Slater exchange coefficient) appears correct on paper,
but has never been validated point-by-point against a reference implementation.

## Iron rules (added to CLAUDE.md)

1. **No fallback values.** Lookup tables must cover all supported elements explicitly. Unsupported input must `throw std::invalid_argument`, never silently return a default.
2. **libxc as XC ground truth.** Native XC output must match libxc to 1e-12 at every test point. Any discrepancy is a bug in native code.

## Strategy

Integrate **libxc** as an optional dependency gated by `ATV_USE_LIBXC`.
Use it as ground truth at two levels:

1. **Point-level XC validation**: feed identical (ρ_α, ρ_β, σ_αα, σ_αβ, σ_ββ) to both native and libxc, assert numerically same output (exc, vrho, vsigma). Fix native code until tests pass.
2. **Grid-level validation**: replace the grid with PySCF-compatible Lebedev + Treutler-Ahlrichs. Compare E_xc on the same density.
3. **End-to-end validation**: UKS total energy must match PySCF within 1e-5 Ha for atoms, 1e-4 Ha for molecules.

## Phases

### Phase 1: libxc CMake integration

**Files:** `CMakeLists.txt`

- Add `option(ATV_USE_LIBXC "Use libxc for XC functionals" OFF)`
- `find_package(Libxc)` or `FetchContent` if not found
- Define `ATV_HAS_LIBXC` compile definition when available
- Link `Libxc::xc` to atomiverse target
- Gate all libxc code behind `#ifdef ATV_HAS_LIBXC`

### Phase 2: libxc XC evaluator wrapper

**New files:**
- `include/atomiverse/cpu/scf/xc_libxc.h`
- `src/cpu/scf/xc_libxc.cpp`

Provide a function with the same signature as `evaluate_xc()`:

```cpp
#ifdef ATV_HAS_LIBXC
XCOutput evaluate_xc_libxc(XCType type,
                           Real rho_a, Real rho_b,
                           Real sigma_aa, Real sigma_ab, Real sigma_bb);
#endif
```

Implementation:
- LSDA → `XC_LDA_X` (Slater) + `XC_LDA_C_VWN` (VWN5, functional id 7, NOT VWN3=id 5)
- PBE → `XC_GGA_X_PBE` + `XC_GGA_C_PBE`
- Call `xc_lda_exc_vxc()` / `xc_gga_exc_vxc()` with spin-polarized inputs
- Convert libxc output convention to Atomiverse's `XCOutput` struct

**Key libxc conventions to handle:**
- libxc `exc` is energy per particle (per electron), Atomiverse uses same
- libxc `vrho[0]` = ∂(ρ·exc)/∂ρ_α, `vrho[1]` = ∂(ρ·exc)/∂ρ_β
- libxc `vsigma[0]` = ∂(ρ·exc)/∂σ_αα, `[1]` = ∂/∂σ_αβ, `[2]` = ∂/∂σ_ββ
- Verify sign conventions match

### Phase 3: Point-level XC tests — fix native implementation

**New file:** `tests/test_xc_vs_libxc.cpp`

Test plan:
```
FOR each functional in {LSDA, PBE}:
  FOR each test point in {
    unpolarized (ρ_α=ρ_β=0.5),
    polarized (ρ_α=0.8, ρ_β=0.2),
    fully polarized (ρ_α=1.0, ρ_β=0),
    low density (ρ_α=1e-6, ρ_β=1e-8),
    high density (ρ_α=100, ρ_β=50),
    GGA: various σ values
  }:
    native = evaluate_xc(type, ρ_α, ρ_β, ...)
    libxc  = evaluate_xc_libxc(type, ρ_α, ρ_β, ...)
    ASSERT |native.exc - libxc.exc| < 1e-12
    ASSERT |native.vrho_alpha - libxc.vrho_alpha| < 1e-12
    ...
```

**Expected fixes needed in `xc_functional.cpp`:**
- Verify `kCx` value (Slater constant) matches libxc's `XC_LDA_X`
- Verify VWN5 correlation gives exact same numbers as libxc `XC_LDA_C_VWN` (id=7)
- Check the spin-interpolation formula: `f(ζ)` and its derivative — common source of subtle bugs
- Check `vrho` convention: is it `∂exc/∂ρ_σ` or `∂(ρ·exc)/∂ρ_σ`? Must match libxc
- PBE: verify κ, μ, β, γ constants against libxc source

### Phase 4: Fix angular grid — real Lebedev tabulation

**File:** `grid.cpp`

Replace the product grid with proper Lebedev quadrature tables.

Source for Lebedev data:
- libxc ships with Lebedev tables, or
- Dmitri Laikov's tables (public domain): 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810
- PySCF's `pyscf/lib/dft/grid_basis.c` has them too

For Phase 4, implement at least: 6, 50, 110, 194, 302, 590-point grids.

Storage format: static const arrays of (x, y, z, w) tuples, generated by a Python script that pulls from Lebedev-Laikov tables. Store in a separate file `lebedev_data.h` (auto-generated, committed).

**Validation test:** `tests/test_grid.cpp`
```
// Spherical harmonic orthogonality: ∫ Y_l^m · Y_l'^m' dΩ = δ_{ll'} δ_{mm'}
// Test for l up to the grid's stated exactness degree
FOR n_angular in {6, 110, 302}:
  grid = lebedev(n_angular)
  FOR l in 0..max_degree(n_angular):
    FOR m in -l..l:
      integral = Σ_i w_i * Y_l^m(θ_i, φ_i) * Y_l^m(θ_i, φ_i)
      ASSERT |integral - 1.0| < 1e-14
```

### Phase 5: Fix radial grid — Mura-Knowles or Treutler-Ahlrichs

**File:** `grid.cpp`

Option A: **Mura-Knowles** (MK, simpler, PySCF default for non-pruned grids):
```
x_i = (i + 0.5) / n   for i = 0, ..., n-1
r_i = -R * ln(1 - x_i^3)
w_i = R * 3 * x_i^2 / (1 - x_i^3) / n * r_i^2
```

Option B: **Treutler-Ahlrichs** (TA, used by PySCF with pruning):
```
x_i = Gauss-Chebyshev-II points on [-1, 1]
r_i = R / ln(2) * (1 + x_i)^α * ln(2 / (1 - x_i))
```

Recommend **Option A (Mura-Knowles)** for simplicity and PySCF compatibility. Keep Becke mapping as fallback behind a runtime option.

**Validation test:**
```
// Integrate r^n * exp(-α*r^2) * 4π r^2 dr analytically and compare
// Known: ∫_0^∞ r^{2+n} exp(-αr²) dr = Γ((n+3)/2) / (2 α^{(n+3)/2})
FOR n in {0, 2, 4, 6}:
  FOR α in {0.1, 1.0, 10.0}:
    numerical = Σ_i w_i * r_i^n * exp(-α * r_i^2)
    exact = gamma((n+3)/2) / (2 * α^((n+3)/2))
    ASSERT |numerical - exact| / exact < 1e-10
```

### Phase 6: Complete Bragg-Slater radii table

**File:** `grid.cpp`

Replace the 3-element table with the complete Bragg-Slater set (Z=1..118).
Source: J.C. Slater, JCP 41, 3199 (1964), extended by Cordero et al.

At minimum, cover Z=1..36 (first two periods + transition metals).

Values (Angstrom): H=0.25, He=**0.31**, Li=1.45, Be=1.05, B=0.85, C=0.70, N=0.65, O=0.60, F=0.50, Ne=**0.38**, Na=1.80, Mg=1.50, Al=1.25, Si=1.10, P=1.00, S=1.00, Cl=1.00, Ar=**0.71**, ...

Note: He=0.31 (not 0.25). The current 0.25 is wrong.

### Phase 7: End-to-end UKS validation

**File:** `tests/test_uks.cpp` (update), `python/examples/03_uks.py` (update)

Validation targets (with libxc + fixed grid):
| System | Expected |dE| vs PySCF |
|--------|--------------------------|
| H atom | < 1e-6 Ha |
| He atom | < 1e-5 Ha |
| Li atom | < 1e-4 Ha |
| H2 eq | < 1e-5 Ha |
| H2 PES scan | < 1e-4 Ha at all points |

Update `03_uks.py` threshold from "KNOWN OFFSET" to PASS/FAIL with 1% threshold.

### Phase 8 (optional): Runtime libxc fallback

Add a runtime switch in `uks_solve` to use libxc evaluator instead of native:

```cpp
XCOutput xc;
#ifdef ATV_HAS_LIBXC
if (use_libxc)
    xc = evaluate_xc_libxc(type, rho_a, rho_b, ...);
else
#endif
    xc = evaluate_xc(type, rho_a, rho_b, ...);
```

This allows users to choose libxc for production accuracy or native for educational/debugging purposes.

## File manifest

| File | Action | Phase |
|------|--------|-------|
| `CMakeLists.txt` | Add `ATV_USE_LIBXC` option, find/fetch libxc | 1 |
| `include/atomiverse/cpu/scf/xc_libxc.h` | New: libxc wrapper header | 2 |
| `src/cpu/scf/xc_libxc.cpp` | New: libxc wrapper implementation | 2 |
| `tests/test_xc_vs_libxc.cpp` | New: point-level XC comparison | 3 |
| `src/cpu/scf/xc_functional.cpp` | Fix: match libxc output exactly | 3 |
| `include/atomiverse/cpu/scf/lebedev_data.h` | New: tabulated Lebedev points | 4 |
| `src/cpu/scf/grid.cpp` | Fix: real Lebedev, Mura-Knowles radial, complete radii | 4-6 |
| `tests/test_grid.cpp` | Update: spherical harmonic + radial accuracy tests | 4-5 |
| `tests/test_uks.cpp` | Update: end-to-end energy comparison | 7 |
| `python/examples/03_uks.py` | Update: change KNOWN OFFSET to PASS/FAIL | 7 |
| `src/cpu/scf/scf_solve.cpp` | Optional: runtime libxc switch | 8 |

## Success criteria

1. `test_xc_vs_libxc`: all XC point comparisons match to 1e-12
2. `test_grid`: spherical harmonic integrals exact to grid's stated degree
3. `test_uks`: He atom energy matches PySCF to 1e-5 Ha
4. `03_uks.py`: all examples PASS with 1% threshold
5. Native XC without libxc gives same results as with libxc (within 1e-12)

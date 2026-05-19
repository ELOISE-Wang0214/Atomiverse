---
slug: libxc-xc-fix-05-uks-validation
criteria:
  - id: ac-001
    summary: UKS-LSDA on 5 systems within per-system Ha tolerance vs PySCF
    type: scientific
    pass_when: |
      Running ctest target test_uks_solve from a clean build:
      |E_atv - E_PySCF_const| <= 1e-6 Ha for H atom,
      <= 1e-5 Ha for He atom,
      <= 1e-4 Ha for Li atom,
      <= 1e-5 Ha for H2 at 0.74 Å,
      and <= 1e-4 Ha at each R in the H2 PES scan
      {0.5, 0.6, 0.74, 0.85, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5} Å,
      with PySCF reference values present as constexpr double constants
      in tests/cpu/scf/test_uks_solve.cpp and a header comment naming
      the PySCF version, basis, grid spec, and conv_tol that produced
      them.
    status: pending
  - id: ac-002
    summary: Full ctest suite passes in the default build configuration
    type: runtime
    pass_when: |
      `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
      && cmake --build build -j$(nproc)
      && ctest --test-dir build --output-on-failure`
      exits 0 with zero failing tests. No special ATV_BUILD_XC_ORACLE
      toggle is required.
    status: pending
  - id: ac-003
    summary: 03_uks.py exits 0 on green and non-zero on red with PASS/FAIL line
    type: runtime
    pass_when: |
      `python python/examples/03_uks.py` prints a line matching
      `Result: PASS` or `Result: FAIL`, and exits with status 0 iff
      atom_max_rel < 1.0% AND h2_conv_max_rel < 1.0% (the converged
      R <= 1.5 Å subset). Any failure of either condition yields
      non-zero exit status and a `Result: FAIL` line. The docstring
      no longer contains the phrase "KNOWN OFFSET".
    status: pending
---

# Acceptance criteria

## ac-001 — UKS-LSDA per-system tolerance vs PySCF (scientific)

Five system families exercised in `test_uks_solve`: H atom (doublet); He atom (singlet); Li atom (doublet); H₂ at 0.74 Å; H₂ PES at the ten R values. PySCF references are `constexpr double` in the test file, reproducible from the header-comment recipe. The H₂ PES scan is value-parameterized — a violation at any R fails the test and names the offending R.

## ac-002 — Full ctest green in default configuration (runtime)

Clean default build exits zero. No libxc dependency: `ATV_BUILD_XC_ORACLE` is orthogonal to this gate.

## ac-003 — Python example PASS/FAIL gate (runtime)

`03_uks.py` no longer documentation-only. Docstring describes a validation gate (no "KNOWN OFFSET" framing). `main()` ends with `sys.exit(0 if passed else 1)`.

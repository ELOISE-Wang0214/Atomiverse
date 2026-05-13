---
slug: pw-dft-si-04-pw-driver
criteria:
  - id: ac-001
    summary: pw::solve_si() declared in header
    type: code
    pass_when: |
      `include/atomiverse/cpu/pw/solve_si.h` exists with `double pw::solve_si(...)` declaration.
    status: pending
  - id: ac-002
    summary: test_pw_si target when ATV_BUILD_PW=ON
    type: code
    pass_when: |
      cmake -DATV_BUILD_PW=ON with FFTW3 produces `test_pw_si` executable target.
    status: pending
  - id: ac-003
    summary: ATV_BUILD_PW=OFF excludes PW
    type: code
    pass_when: |
      cmake -DATV_BUILD_PW=OFF succeeds, no test_pw_si target, no PW symbols in lib.
    status: pending
  - id: ac-004
    summary: ATV_BUILD_PW=ON without FFTW3 fails clearly
    type: code
    pass_when: |
      cmake fails with message containing "FFTW3".
    status: pending
  - id: ac-005
    summary: test_pw_si compiles, links, runs
    type: runtime
    pass_when: |
      ctest -R test_pw_si exits 0.
    status: pending
  - id: ac-006
    summary: Si total energy within 0.5 Ry of -107.8 Ry
    type: scientific
    pass_when: |
      pw::solve_si() returns E such that |E + 107.8| < 0.5 Ry.
    status: pending
  - id: ac-007
    summary: Doxygen docstring with units for solve_si()
    type: docs
    pass_when: |
      @param and @return documented with energy units (Ry).
    status: pending
---

---
slug: regressions-harness
criteria:
  - id: ac-001
    summary: regressions/ builds behind its own switch and never runs in default ctest
    type: code
    pass_when: |
      cmake -DATV_BUILD_REGRESSIONS=ON configures; `cmake --build build --target
      regressions` builds AND runs every registered case. With the switch OFF the
      regressions tree is not configured. `ctest` (default) registers ZERO
      regression cases — grep of the ctest test list shows no reg_* entry. This
      preserves CLAUDE.md's "tests/ is unit-only and if it is registered it runs".
    status: verified
    note: |
      ctest -N on a REGRESSIONS=ON build reports Total Tests: 0; regressions/CMakeLists.txt contains no add_test(). `--target regressions` builds and runs every case.

  - id: ac-002
    summary: Regression cases are out-of-library consumers of the public API only
    type: code
    pass_when: |
      Every regressions/**/reg_*.cpp includes only public `atomiverse/...` headers
      and links only the `atomiverse` target. A grep for `internal/` or any
      relative-path include into src/ over regressions/ returns nothing. A case
      that cannot be written against the public API is reported as an API gap, not
      worked around by reaching into internals.
    status: verified
    note: |
      reg_lj_argon.cpp includes 8 public atomiverse/ headers and links only the atomiverse target; grep finds no internal/ or ../src reach-in.

  - id: ac-003
    summary: Every case cites a literature reference, not a self-generated golden value
    type: science
    pass_when: |
      Each regressions/<domain>/<case>/ carries a REFERENCE.md naming the source
      (book/paper + section/equation/page), the expected value, and the physical
      justification for the tolerance. No expected value is a frozen snapshot of
      Atomiverse's own current output. Reviewed by reading REFERENCE.md.
    status: verified
    note: |
      REFERENCE.md cites Rahman 1964 (Phys. Rev. 136, A405) for the model and Frenkel&Smit / Allen&Tildesley / Hairer-Lubich-Wanner for the symplectic-conservation criterion. The tolerance is an envelope from theory, not a frozen snapshot of our output.

  - id: ac-004
    summary: Seed case — LJ argon NVE energy conservation — passes and is fast
    type: scientific
    pass_when: |
      regressions/md/lj_argon/reg_lj_argon.cpp runs LJ argon (Rahman 1964
      parameters: eps/kB = 120 K, sigma = 3.4 Angstrom) under velocity Verlet NVE
      and asserts total-energy drift stays below the documented bound over the run;
      exits non-zero on violation and prints got/want/tol. Wall time <= 60 s
      single-threaded CPU.
    status: verified
    note: |
      PASS. N=216 Rahman argon, dt=5 fs, 10 ps: worst |dE|/(N*eps) = 1.795e-04 vs bound 1e-2; second-half drift (1.707e-04) does not exceed the first (1.795e-04) => bounded and non-secular. Wall 0.11 s single-threaded (limit 60 s).

  - id: ac-005
    summary: The three-way tests / benchmarks / regressions split is documented
    type: code
    pass_when: |
      CLAUDE.md's Repository layout lists regressions/ with its one-line purpose,
      and the tests-vs-benchmarks rule is restated as a three-way split:
      tests/ = unit only, WHOLE SUITE RUNS IN SECONDS (the fast-iteration loop),
      one function/method per test, E2E STRICTLY FORBIDDEN;
      benchmarks/ = external-ENGINE parity + performance;
      regressions/ = out-of-library end-to-end vs the LITERATURE, <= 60 s,
      and EVERY delivered feature ships exactly one.
      regressions/README.md carries the hard rules.
    status: verified
    note: |
      CLAUDE.md Repository layout now lists regressions/, and Conventions carries the 'Three test tiers, and they never mix' rule (tests/ = unit, seconds, e2e STRICTLY FORBIDDEN). regressions/README.md carries the six hard rules.

  - id: ac-006
    summary: Regression cases are flat main()s — no single-use helper functions
    type: code
    pass_when: |
      Per the user directive "只运行一次的代码不要抽出来作为单独的函数": each
      reg_*.cpp is one file whose setup / run / check reads top-to-bottom inside
      main(). It does NOT carve one-shot code into setup_frame() / run_md() /
      check_energy() helpers that each have exactly one call site. Reviewed by
      reading the files; a helper is only acceptable with >= 2 call sites.
    status: verified
    note: |
      reg_lj_argon.cpp is one flat main(): setup, run, and checks read top to bottom, with zero single-use helper functions.
---

# Acceptance — regressions-harness

Gate for every downstream breadth spec: none of `manybody-potentials`,
`md-constraints`, `path-driver-neb`, `mc-gcmc-npt` may reach `code-complete`
without landing its own passing case under `regressions/`.

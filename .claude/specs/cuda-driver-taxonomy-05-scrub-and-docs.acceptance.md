---
slug: cuda-driver-taxonomy-05-scrub-and-docs
criteria:
  - id: ac-001
    summary: Dangling CUDA scf::Phase instantiation removed; CUDA tree free of scf::Phase
    type: code
    pass_when: |
      src/cuda/driver.cu contains no `Driver<scf::Phase` instantiation,
      and `grep -rn 'scf::Phase' include/atomiverse/cuda/ src/cuda/`
      returns zero matches.
    status: pending
  - id: ac-002
    summary: ATV::scf::Phase enum kept and its doc-comment re-pointed to cpu::pw
    type: code
    pass_when: |
      include/atomiverse/facade/stage.h still declares `namespace ATV::scf`
      with `enum class Phase` and the `phase_names(Phase)` overload, and its
      doc-comment no longer contains the strings "Legacy CUDA-only" or
      "will retire", instead naming the CPU plane-wave SCF (`cpu::pw`) as the
      owner.
    status: pending
  - id: ac-003
    summary: CPU plane-wave scf::Phase instantiation left intact
    type: code
    pass_when: |
      src/cpu/driver.cpp still contains
      `template class Driver<scf::Phase, cpu::device_tag>;`.
    status: pending
  - id: ac-004
    summary: CLAUDE.md Driver-taxonomy table shows migrated CUDA shape
    type: code
    pass_when: |
      The "Driver taxonomy" table in CLAUDE.md lists the CUDA rows as the
      three-category shape (cuda::Evaluate / cuda::md::MD /
      cuda::minimizer::Minimize) with the cuda::scf::SCFSolver inner solver,
      and no CUDA row contains the strings "legacy", "pending migration", or
      "algorithm-as-driver".
    status: pending
  - id: ac-005
    summary: Full CUDA build and CUDA test suite green after the scrub
    type: runtime
    evaluator_hint: requires CUDA toolkit
    pass_when: |
      `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON` then
      `cmake --build build -j$(nproc)` succeed, and
      `ctest --test-dir build -R cuda --output-on-failure` reports all
      CUDA tests passing (requires a CUDA toolkit).
    status: pending
---

# Acceptance criteria

- **ac-001** — The mechanical core of the link. The single dead explicit
  instantiation in `src/cuda/driver.cu` is gone, and the CUDA tree (headers +
  sources) carries no surviving `scf::Phase` token, including `@see` and usage-
  snippet doc-comments re-pointed to the `cuda::scf::SCFSolver` API. The grep gate
  is the RED test from the first Task — failing before the scrub, passing after.

- **ac-002** — Guards the corrected scope: `ATV::scf::Phase` is **not** deleted
  (its CPU plane-wave owner still needs it), but its doc-comment is corrected from
  the now-false "Legacy CUDA-only … will retire" to name `cpu::pw` as the owner.

- **ac-003** — Regression guard that the CPU plane-wave SCF instantiation was not
  collaterally removed while scrubbing the CUDA one.

- **ac-004** — Documentation truthfulness: CLAUDE.md's Driver taxonomy table
  reflects the shipped migrated CUDA shape with no stale "legacy/pending" wording.
  A grep handle on the table makes this third-party-checkable.

- **ac-005** — End-to-end no-behavior-change guarantee: the full CUDA build links
  without the dead instantiation and the entire CUDA test suite still passes.
  `runtime` type because it requires building and executing on a CUDA toolkit.

---
title: CUDA driver taxonomy — link 05: scrub dangling scf::Phase and document the migration
status: approved
created: 2026-05-22
follows: cuda-driver-taxonomy-04-scf
---

# CUDA driver taxonomy — link 05: scrub dangling scf::Phase and document the migration

## Summary

This is the fifth and final link of the chain that migrated the CUDA backend onto
the three-category driver taxonomy already shipped on CPU. By the time this link
runs, links 01–04 have landed: CUDA is fully namespaced, `cuda::md::MD` /
`cuda::minimizer::Minimize` / `cuda::Evaluate` are the three category drivers, and
`cuda::scf::SCF` has been demoted to the `cuda::scf::SCFSolver` inner solver. This
link removes the now-dangling CUDA references to the legacy `ATV::scf::Phase` enum,
re-points the enum's documentation to its one genuine remaining owner, and updates
project documentation (CLAUDE.md, the engineering notes log) so the recorded
architecture matches the shipped reality. There is no behavior change: the build
and the full CUDA test suite must compile and pass exactly as before.

## Domain basis

None. This link contains no equations and no physics — it removes a dead explicit
instantiation, re-points doc-comments, and updates project documentation. No
behavior changes; correctness is the unchanged CUDA build + test suite.

## Design

**Correction to the chain's planning assumption.** The chain plan and this link's
brief both assumed `ATV::scf::Phase` (`facade/stage.h`, namespace `ATV::scf`)
would be referenced by nothing once link 04 demoted CUDA SCF to an inner solver,
and could therefore be deleted outright. A tree-wide grep performed while drafting
this spec shows that is false: the CPU **plane-wave** SCF subsystem is a live,
non-legacy consumer.

- `src/cpu/driver.cpp:12` explicitly instantiates `Driver<scf::Phase, cpu::device_tag>`.
- Seven `cpu::pw::PW_*` components (`pw_poisson.h`, `pw_hamiltonian.h`, `pw_lda.h`,
  `pw_veff.h`, `pw_density.h`, `pw_convergence.h`, `pw_davidson.h`) return
  `DefaultBinding` values keyed on `ATV::scf::Phase::*`.

The Gaussian-basis CPU SCF is indeed an inner solver and does not use the enum —
but the plane-wave SCF runs as a `Driver<scf::Phase>` outer loop. The brief's
description was accurate only for the Gaussian path.

Consequently `ATV::scf::Phase` **must not be deleted** in this link. Deleting it
would break the CPU plane-wave build. Migrating `cpu::pw` off `scf::Phase` onto its
own phase enum / category driver is a genuine architectural refactor with behavior
risk — out of scope for a mechanical hygiene link and out of scope for this chain
(which is a CUDA migration). It is recorded as a follow-up in the notes log.

What this link **does** do:

1. **Scrub the dangling CUDA instantiation.** `src/cuda/driver.cu:118`
   (`template class Driver<scf::Phase, cuda::device_tag>;`) is now dead after
   link 04 — no CUDA driver is parameterized on `scf::Phase` anymore. Remove that
   one line. The CPU instantiation in `src/cpu/driver.cpp:12` stays.

2. **Scrub CUDA-side `scf::Phase` doc references.** Six CUDA SCF headers carry
   `/// @see scf::Phase::*` comments (`scf.cuh`, `convergence_check.cuh`,
   `rhf_fock.cuh`, `density_builder.cuh`, `rhf_density.cuh`, `eigen_solver.cuh`,
   `hamiltonian.cuh`) plus three `solver.on(scf::Phase::end_of_step)` usage
   snippets in doc-comments (`scf.cuh`, `density_recorder.cuh`, `mulliken.cuh`).
   After link 04 these classes are `cuda::scf::internal::*` step helpers, not
   phase-bound components. Re-point the `@see` lines and usage snippets to the
   `cuda::scf::SCFSolver` inner-solver API. This removes every CUDA reference to
   `scf::Phase`.

3. **Re-point the enum doc-comment.** `facade/stage.h`'s `namespace ATV::scf`
   doc-comment currently reads "Legacy CUDA-only enum … This namespace will retire
   once the CUDA SCF migration ships." That is now wrong on both counts: it is not
   CUDA-only (plane-wave CPU SCF owns it) and it is not retiring. Rewrite the
   doc-comment to state that `ATV::scf::Phase` is the phase enum for the CPU
   plane-wave SCF outer driver (`cpu::pw`), and that Gaussian-basis CPU SCF and the
   migrated CUDA SCF are inner solvers that do not use it. The enum body,
   `phase_names()` overload, and namespace are otherwise untouched.

4. **Update CLAUDE.md.** The "Driver taxonomy" table CUDA rows change from
   "legacy — pending migration" / "legacy — algorithm-as-driver" to the migrated
   three-category shape; the Phase-naming line drops the "CUDA-only legacy" note
   and restates `scf::Phase` as the plane-wave SCF enum; the Conventions namespace
   list and the "Implemented components" tables reflect the migrated CUDA shape; a
   short CUDA single-point / MD snippet may be added to the example block.

5. **Append an engineering-notes entry** dated 2026-05-22 to `.claude/notes/notes.md`
   recording the CUDA driver taxonomy migration, mirroring the CPU entry style, and
   noting the two latent bugs fixed in link 01 plus the plane-wave-SCF follow-up
   surfaced here.

No symbols are added or renamed. The only deleted symbol-level artifact is one
explicit-instantiation line.

## Files to create or modify

- `include/atomiverse/facade/stage.h` — rewrite the `namespace ATV::scf`
  doc-comment (lines ~87–93); enum body and `phase_names()` unchanged.
- `src/cuda/driver.cu` — delete the `Driver<scf::Phase, cuda::device_tag>`
  explicit instantiation (line 118).
- `include/atomiverse/cuda/scf/scf.cuh` — re-point `@see scf::Phase` and the
  `solver.on(scf::Phase::end_of_step)` usage snippet to the `SCFSolver` API.
- `include/atomiverse/cuda/scf/convergence_check.cuh` — re-point `@see scf::Phase::convergence_check`.
- `include/atomiverse/cuda/scf/rhf_fock.cuh` — re-point `@see scf::Phase::build_hamiltonian`.
- `include/atomiverse/cuda/scf/density_builder.cuh` — re-point `@see scf::Phase::build_density`.
- `include/atomiverse/cuda/scf/rhf_density.cuh` — re-point `@see scf::Phase::build_density`.
- `include/atomiverse/cuda/scf/eigen_solver.cuh` — re-point `@see scf::Phase::diagonalize`.
- `include/atomiverse/cuda/scf/hamiltonian.cuh` — re-point `@see scf::Phase::build_hamiltonian`.
- `include/atomiverse/cuda/recorder/density_recorder.cuh` — re-point the
  `solver.on(scf::Phase::end_of_step)` usage snippet.
- `include/atomiverse/cuda/observer/mulliken.cuh` — re-point the
  `solver.on(scf::Phase::end_of_step)` usage snippet.
- `CLAUDE.md` — Driver taxonomy table, Phase-naming line, Conventions namespace
  list, Implemented components tables, example block.
- `.claude/notes/notes.md` — append the 2026-05-22 CUDA driver taxonomy entry.

> Note: links 02–04 land the `cuda::scf::internal::` relocation of the SCF step
> helpers; if those headers no longer exist under `cuda/scf/` by the time this
> link runs, the corresponding `@see` re-point items are satisfied by link 04 and
> this link only needs to confirm no `scf::Phase` token survives in the CUDA tree.

## Tasks

- [ ] Write a grep gate verifying `scf::Phase` is absent from CUDA code (tests/cuda/scf_phase_scrub.sh or a CMake/CTest check over include/atomiverse/cuda/ + src/cuda/) — RED before the scrub
- [ ] Delete the `Driver<scf::Phase, cuda::device_tag>` explicit instantiation in src/cuda/driver.cu
- [ ] Re-point all CUDA-side `@see scf::Phase` and `solver.on(scf::Phase::...)` doc references (cuda/scf/*.cuh, cuda/recorder/density_recorder.cuh, cuda/observer/mulliken.cuh) to the cuda::scf::SCFSolver API
- [ ] Rewrite the `namespace ATV::scf` doc-comment in include/atomiverse/facade/stage.h to identify cpu::pw plane-wave SCF as its owner and drop the "Legacy CUDA-only / will retire" wording
- [ ] Update CLAUDE.md: Driver taxonomy table CUDA rows, Phase-naming line, Conventions namespace list, Implemented components tables, example block
- [ ] Append the 2026-05-22 CUDA driver taxonomy entry to .claude/notes/notes.md (mirror the CPU "[2026-05-06] CPU driver taxonomy redesign" style; record the two latent bugs fixed and the plane-wave-SCF follow-up)
- [ ] Run full check + test suite

## Testing strategy

- **Happy path / RED gate:** a grep gate asserts `grep -rn 'scf::Phase'` over
  `include/atomiverse/cuda/` and `src/cuda/` returns no matches. It fails before
  the scrub (RED) and passes after (GREEN). It deliberately does **not** cover
  `include/atomiverse/cpu/` or `src/cpu/`, where `cpu::pw` legitimately keeps
  `scf::Phase`.
- **Build integrity:** the full project build with CUDA enabled
  (`cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON`, `cmake --build build`)
  must compile after removing the dangling instantiation — confirming nothing in
  the CUDA tree still depended on `Driver<scf::Phase, cuda::device_tag>`.
- **CPU plane-wave regression:** the CPU build and any plane-wave SCF tests must
  still compile and pass — confirms the kept `ATV::scf::Phase` enum and the
  `src/cpu/driver.cpp` instantiation remain intact and that the doc-comment rewrite
  did not touch the enum body.
- **Edge case:** `grep -rn 'Driver<scf::Phase, cuda'` over the whole tree returns
  nothing (no other CUDA instantiation site exists).
- **Full CUDA test suite:** `ctest --test-dir build -R cuda --output-on-failure`
  green — no behavior change. Requires a CUDA toolkit.
- **Docs:** CLAUDE.md "Driver taxonomy" table contains no "legacy" / "pending
  migration" / "algorithm-as-driver" wording in its CUDA rows; the notes entry
  dated 2026-05-22 exists in `.claude/notes/notes.md`.

## Out of scope

- **Deleting `ATV::scf::Phase`.** Surfaced during drafting as still owned by the
  CPU plane-wave SCF (`cpu::pw`); deletion would break that build. Kept.
- **Migrating `cpu::pw` plane-wave SCF off `scf::Phase`** onto its own phase enum /
  category driver. A real architectural refactor with behavior risk, unrelated to
  the CUDA migration this chain performs. Recorded as a follow-up in the notes
  entry; not done here.
- **Any remaining CUDA performance work.**
- **The CPU `ResponseDriver` / `vib::Phase` (4th CPU category).** Unrelated and
  untouched.
- **`docs/developer/cpu-internals.md` and `.claude/architecture.md`** — those carry
  stale `scf::Phase` prose too, but they are generated/architectural docs outside
  the CLAUDE.md + notes scope this link owns; left for a docs-sweep link.

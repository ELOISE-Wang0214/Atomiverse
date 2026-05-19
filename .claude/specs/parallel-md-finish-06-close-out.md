---
title: Close out parallel-md-engine — full suite green + retire monolith
slug: parallel-md-finish-06-close-out
status: approved
created: 2026-05-19
parent: parallel-md-engine.md
depends-on:
  - parallel-md-finish-01-profile-csvs
  - parallel-md-finish-02-lammps-reference
  - parallel-md-finish-03-atomiverse-8rank
  - parallel-md-finish-04-scaling-report
  - parallel-md-finish-05-omp-numa
---

# Close out parallel-md-engine — full suite green + retire monolith

## Summary

Final closure step: build the full source tree, run the full `cpu` + `mpi` ctest labels, verify zero regressions vs the pre-Phase-5 test count, flip `parallel-md-engine.md`'s frontmatter `status` from `in-progress` → `retired`, update the INDEX entry, and write the chain-completion summary. Closes parent ac-025 and retires the parent spec. After this lands, all six chain sub-specs can themselves be flipped to `retired` in a single follow-up cleanup commit.

## Domain basis

No physics. Pure orchestration + documentation.

## Design

Touched entities:

- `.claude/specs/parallel-md-engine.md` — flip `status: in-progress` → `status: retired` in frontmatter; add a closure note at the top of the body pointing to the chain.
- `.claude/specs/INDEX.md` — change the parallel-md-engine line to `[retired — see parallel-md-finish-NN]` and list the 6 chain sub-spec entries indented underneath.
- The chain sub-spec files themselves — frontmatter `status` should already be at `code-complete` from prior `/mol:impl` runs; this sub-spec flips them all to `retired` once `ac-001` passes (full suite green).

No source code changes.

## Files to create or modify

- `.claude/specs/parallel-md-engine.md` (modify frontmatter + closure note)
- `.claude/specs/INDEX.md` (modify parallel-md line, indent chain children)
- `.claude/specs/parallel-md-finish-01-profile-csvs.md` (modify frontmatter)
- `.claude/specs/parallel-md-finish-02-lammps-reference.md` (modify frontmatter)
- `.claude/specs/parallel-md-finish-03-atomiverse-8rank.md` (modify frontmatter)
- `.claude/specs/parallel-md-finish-04-scaling-report.md` (modify frontmatter)
- `.claude/specs/parallel-md-finish-05-omp-numa.md` (modify frontmatter)
- `.claude/specs/parallel-md-finish-06-close-out.md` (modify frontmatter — self)

## Tasks

- [ ] Verify all five upstream chain sub-specs have `status: code-complete` in their frontmatter; abort if any is still at `approved` or `in-progress` (means the corresponding `/mol:impl` run did not finish).
- [ ] Build full tree fresh: `cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_TESTS=ON -DATV_BUILD_BENCHMARKS=ON -DATV_BUILD_RHF=ON -DATV_USE_MPI=ON && cmake --build build -j`.
- [ ] Run the full suite: `ctest --test-dir build -L "cpu|mpi" --output-on-failure`. Verify all tests pass and the test count is ≥ the pre-Phase-5 baseline recorded in the parent spec status snapshot (Phase 1–4 exit: 20 MPI tests + the existing CPU suite). Two new tests must appear in the count: `test_nose_hoover_t_ramp` and `test_md_no_alloc` (the latter must report SKIPPED on this Darwin host, which is still a pass).
- [ ] Flip parent frontmatter: in `.claude/specs/parallel-md-engine.md`, change `status: in-progress` → `status: retired`, update `last-updated:` to today, and add a brief one-paragraph closure note under the existing "Closure note (2026-05-19)" pointing to the SHA of the chain-completion commit.
- [ ] Flip each chain sub-spec's frontmatter `status:` to `retired` once they verify; do this in the same commit as the parent flip.
- [ ] Update `INDEX.md`: replace the parallel-md-engine line with `[retired — see parallel-md-finish-NN chain]`; add the six chain entries as a nested block.
- [ ] Commit: `chore(parallel-md): retire parent spec + chain after full ctest green (closes parent ac-025)`.

## Testing strategy

- **Happy path:** ctest exits 0 with no regressions. All frontmatter flips succeed. INDEX entry is human-readable and links the chain.
- **Edge case:** if a chain sub-spec is still at `code-complete` rather than `retired`, that's still acceptable for this sub-spec's closure — but follow up with a manual retire-pass in a subsequent commit.
- **Edge case:** if any test in the `cpu|mpi` labels fails, this sub-spec halts. Triage the failure first (likely a regression introduced by chain spec 04's bench runs); do not flip to `retired` until green.

## Out of scope

- Re-running benchmarks (sub-spec 04's job).
- Touching `msu-pe-validation.md` (independent spec).
- 32-rank or cluster work (parent Phase 6 deferred).

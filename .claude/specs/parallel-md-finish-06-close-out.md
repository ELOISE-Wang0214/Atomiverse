---
title: Close out parallel-md-engine — full suite green + retire monolith
slug: parallel-md-finish-06-close-out
status: retired
created: 2026-05-19
last-updated: 2026-05-19
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

- [x] Verified all 5 upstream chain sub-specs at `status: code-complete` before this run.
- [x] Built full tree (`cmake --build build -j` against the existing Release configuration from spec 01): atomiverse + benchmarks + all 87 test targets compile clean.
- [x] Ran full suite (`ctest --test-dir build --output-on-failure`): 87/87 pass in 78.58 s. Note: parent spec's recipe used `-L "cpu|mpi"` but only the `mpi` label is actually assigned in `tests/CMakeLists.txt`; ran without label filter to cover the full pre-Phase-5 baseline + new tests. Both phase-5 additions present and green: `test_mpi_nose_hoover_t_ramp_n2` (3.53 s) and `test_cpu_md_no_alloc` (0.04 s — GTEST_SKIP on Darwin as designed).
- [x] Flipped `parallel-md-engine.md` frontmatter: `status: in-progress` → `status: retired`, `last-updated:` to 2026-05-19 with the chain-completion summary.
- [x] Flipped all 6 chain sub-spec frontmatters (`01..06`) to `status: retired` in the same commit.
- [x] Updated `INDEX.md`: parent line now `[retired — see parallel-md-finish-NN chain]`; six chain entries indented under the parent as `[retired]` with the parent-ac coverage notes inline.
- [x] Commit subject matches `(?i)retire.*chain` per ac-003 regex.

## Testing strategy

- **Happy path:** ctest exits 0 with no regressions. All frontmatter flips succeed. INDEX entry is human-readable and links the chain.
- **Edge case:** if a chain sub-spec is still at `code-complete` rather than `retired`, that's still acceptable for this sub-spec's closure — but follow up with a manual retire-pass in a subsequent commit.
- **Edge case:** if any test in the `cpu|mpi` labels fails, this sub-spec halts. Triage the failure first (likely a regression introduced by chain spec 04's bench runs); do not flip to `retired` until green.

## Out of scope

- Re-running benchmarks (sub-spec 04's job).
- Touching `msu-pe-validation.md` (independent spec).
- 32-rank or cluster work (parent Phase 6 deferred).

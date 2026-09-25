# Classical Model Support Notes

## Purpose

This reference owns implementation-scope notes for the EAM and Stillinger–Weber
examples in [Many-Body Potentials](../ch16-many-body-potentials.md). Physical
derivations and exercises belong in that tutorial; backend availability and
parameter-loading restrictions belong here and in the
[roadmap](../../roadmap/index.md).

## Evidence status

The following restrictions were inherited from the earlier tutorial. They
have not been revalidated against a pinned engine revision in this workspace.
They are a checklist for source inspection and execution, not a statement
of current support or a commitment to a future implementation.

- The earlier examples described single-element parameterizations and no
  alloy cross-term implementation. Verify the selected element and parameter
  set; a lookup constructor does not imply a general alloy-file reader.
- The earlier notes described single-rank execution and rejection under MPI.
  Verify the current distributed algorithm rather than assuming that a
  restriction or a later improvement still applies.
- The earlier notes reported no CUDA counterpart and no Tersoff, MEAM, or
  ReaxFF implementation in that documented path. Check the roadmap and
  matching source before interpreting these historical absences as current.

## What to verify

For EAM, forces depend on environment derivatives at both endpoints. In a
distributed algorithm, identify where those quantities are completed and
communicated, not only where coordinate and force halos are exchanged.
The necessary communication depends on the decomposition and algorithm;
a blanket assertion about one halo width is not a substitute for that analysis.

For angular models, verify complete neighbor information around each central
atom, triplet counting, forces on all participating atoms, and stress if the
cell will be changed. Record single-rank and multi-rank comparisons separately.

Attach an exact commit, compiler and backend settings, parameter identifier,
complete input, command, expected quantity, and tolerance to every new support
claim. Distinguish a source declaration, a passing local test, and an external
numerical comparison. Update this note when the implementation evidence changes;
do not duplicate it as an undated limitation in the tutorial.

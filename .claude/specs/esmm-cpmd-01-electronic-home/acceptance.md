---
spec: esmm-cpmd-01-electronic-home
created: 2026-07-26
criteria:
  - id: ac-001
    summary: HF-only build compiles the electronic/ cluster
    type: build
    pass_when: |
      cmake configure+build with -DATV_BUILD_HF=ON -DATV_BUILD_GPW=OFF
      compiles src/cpu/electronic/{cpmd,xl_bomd,aspc_predictor}.cpp and
      src/cpu/linalg/orbital_tangent.cpp with no reference to gpw/internal.
    status: pending
  - id: ac-002
    summary: Old propagator headers deleted, zero aliases
    type: code
    pass_when: |
      include/atomiverse/cpu/gpw/{cpmd.h,xl_bomd.h} and
      include/atomiverse/cpu/scf/aspc_predictor.h do not exist; grep over
      include/ src/ tests/ regressions/ finds no gpw::Cpmd, gpw::XlBomd,
      scf::ASPCPredictor, or includes of the deleted paths.
    status: pending
  - id: ac-003
    summary: Namespace is ATV::cpu::electronic and linalg primitives renamed
    type: code
    pass_when: |
      All five propagator classes live in namespace ATV::cpu::electronic;
      cpu/linalg/orbital_tangent.h declares band_gradient / s_tangent /
      lowdin_orthonormalize (+ _k twins); gpw/internal/ot.h retains only the
      OT parameterization, preconditioners, and fixed-H minimizers.
    status: pending
  - id: ac-004
    summary: Whole suite green with unchanged numerics
    type: runtime
    pass_when: |
      Full ctest passes; test_gpw_cpmd, test_gpw_xl_bomd, test_cpu_gpw_ot,
      test_scf_evaluator, test_accelerators pass; the 15 cpmd/sgcp/xl-bomd
      regressions under regressions/dynamics/ pass with their existing pins
      untouched (no tolerance edits in this spec).
    status: pending
  - id: ac-005
    summary: CMake placement is HF-gated propagators + unconditional linalg
    type: code
    pass_when: |
      The three electronic/*.cpp are listed in the ATV_BUILD_HF source block;
      orbital_tangent.cpp is in the unconditional CPU source list; nothing
      electronic/ remains in the ATV_BUILD_GPW block.
    status: pending
  - id: ac-006
    summary: GPW local Angstrom-to-Bohr constant removed
    type: code
    pass_when: |
      src/cpu/gpw/gpw_evaluator.cpp no longer defines a local kAng2Bohr;
      conversions use units::kAngstromInBohr.
    status: pending
out_of_scope:
  - Behavior changes in any propagator
  - New public members on GPWEvaluator or MeanFieldMethod
  - CUDA mirroring
---

# Acceptance — `esmm-cpmd-01-electronic-home`

**Done** means the five electronic-state propagators live in
`cpu/electronic/` (namespace `electronic::`), the generic orbital primitives
live in `cpu/linalg/orbital_tangent.h`, an HF-only build compiles them, and
every existing test and regression passes bit-for-bit.

## AC-001 — layering proof

The whole point of the move: SCF can consume the propagators without GPW.

## AC-002 — no fossils

House rule: no backward compatibility, no forwarding shims.

## AC-003 — names speak function

`electronic::Cpmd`, `linalg::band_gradient` — not origin.

## AC-004 — pure-refactor gate

Bit-for-bit green is the only accepted evidence for a move.

## AC-005 — build-graph honesty

Propagators are HF-tier infrastructure, not GPW-private.

## AC-006 — drive-by debt

Named debt #6 from the campaign plan closes here.

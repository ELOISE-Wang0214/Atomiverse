---
title: Fix Atomiverse PE force-field constants to match the MSU reference
slug: pe-ff-migrate-dreiding-ua
status: approved
created: 2026-05-19
last-updated: 2026-05-19
parent: msu-pe-validation.md
---

# Fix Atomiverse PE force-field constants to match the MSU reference

## Summary

The MSU reference deck `benchmarks/md/reference/polyethylene/lammps_msu_paper.in`
(Tschopp 2010, downloaded from CAVS ICME) is the ground truth for the
Dreiding-UA polyethylene force field. The Atomiverse PE code in
`benchmarks/md/polyethylene_tension.cpp`, `benchmarks/md/pe_profile_single.cpp`,
and 13 PE-related tests **uses wrong constants** — OPLS-AA-like
(θ₀ = 112.7°, k_θ = 58.35 kcal/mol/rad², OPLS dihedral coefficients) that
do not match the MSU paper. This spec replaces those wrong constants
with the MSU paper values (θ₀ = 109.5°, k_θ = 60, K₁..K₄ = 1.73, −4.49,
0.776, 6.99) in 15 source files. **Pure code change** — no LAMMPS
re-runs, no edits to the MSU reference files (those are ground truth and
immutable).

Outcome: Atomiverse and `lammps_msu_paper.in` run the same physics.
`msu-pe-validation.md` Phase D (quantitative ±5% LAMMPS comparison)
becomes runnable.

## Domain basis

Force-field parameters verbatim from
`benchmarks/md/reference/polyethylene/lammps_msu_paper.in`:

```
bond_coeff      1  350     1.53           # k [kcal/mol/Å²], r0 [Å]
angle_coeff     1   60   109.5             # k [kcal/mol/rad²], θ0 [°]
dihedral_coeff  1  1.73 -4.49 0.776 6.99 0.0   # multi/harmonic
pair_coeff      1 1  0.112   4.01           # ε [kcal/mol], σ [Å]
```

Compared to the wrong constants currently in Atomiverse:

| Term | Atomiverse (wrong) | MSU paper (ground truth) |
|---|---|---|
| `bond r0` | 1.529 Å | 1.53 Å (sub-Å rounding; effectively identical) |
| `angle k` | 58.35 kcal/mol/rad² | **60 kcal/mol/rad²** |
| `angle θ0` | 112.7° | **109.5°** (sp³ tetrahedral) |
| `dihedral K1..K4` | 2.0058, −2.5128, 0.2710, 6.2900 | **1.73, −4.49, 0.776, 6.99** |
| `LJ ε, σ` | 0.112, 4.01 | 0.112, 4.01 ✓ |

Reference: Tschopp et al. 2010, Hossain–Tschopp–Horstemeyer 2010
*Polymer*. Citation in `SOURCES.md`.

## Design

This spec is **literal value replacement across 15 source files** and a
status flip on the parent spec. No reference data is regenerated; no
LAMMPS runs are required. The MSU reference files
(`lammps_msu_paper.in`, `msu_pe.data`) are ground truth and must not be
edited.

| File | What changes |
|---|---|
| `benchmarks/md/polyethylene_tension.cpp` | `kAngleK_kcal` 58.35 → 60, `kAngleTheta0` 112.7° → 109.5°, `kDihedralC` array → 1.73, −4.49, 0.776, 6.99, 0.0. |
| `benchmarks/md/pe_profile_single.cpp` | Same three constants (near-duplicate constants block). |
| `tests/cpu/md/test_pe_single_chain.cpp`, `test_pe_multi_chain.cpp`, `test_pe_nve.cpp`, `test_pe_nvt.cpp`, `test_lj_exclusions.cpp` | Same three constants. Tests use formula-based oracles (V = k·Δθ², NVE drift, T equipartition); no numerical-reference recalibration. |
| `tests/mpi/test_pe_*.cpp` (8 files), `test_nose_hoover_t_ramp.cpp` | Same three constants. Multi-rank tests compare 1-rank vs N-rank, so absolute values don't matter; constants need only be consistent. |
| `benchmarks/md/reference/polyethylene/SOURCES.md` | Flip the "Atomiverse FF discrepancy" subsection to "resolved" with the commit SHA. |
| `.claude/specs/msu-pe-validation.md` | Status snapshot Phase D row: `🟡 blocked on FF migration` → `🟡 ready to run`. |

No new files, no header churn, no API surface change. The architectural
audit listed in earlier drafts of this spec **was already completed** as
a side effect of `pe-potential-abstraction` —
`.claude/notes/potential-abstraction-audit.md` covers it; no re-audit
needed here.

### Implementation notes

1. **Single-pass replacement.** A targeted regex/sed pass on the three
   literal numbers (`58.35`, `112.7`, plus the four OPLS dihedral
   coefficients) across the 15 files is mechanically safe — these
   literals do not appear elsewhere in the codebase. Hand-check the
   diff before committing.
2. **Test recalibration cost: zero.** Atomiverse PE tests use
   formula-based oracles (V_angle = k_θ·Δθ² at small displacement, NVE
   drift `< 1e-4`, NVT temperature within `±10 K` of setpoint). No test
   embeds numerical reference values from a specific LAMMPS run that
   would need updating. Expected outcome: every existing PE test still
   passes.
3. **No LAMMPS run required.** The MSU reference deck and data file
   live in this repo unchanged. Quantitative LAMMPS-vs-Atomiverse
   comparison is `msu-pe-validation.md` Phase D's job; this spec only
   fixes the code.

## Files to create or modify

Modified only — no new files:

```
benchmarks/md/polyethylene_tension.cpp
benchmarks/md/pe_profile_single.cpp
tests/cpu/md/test_lj_exclusions.cpp
tests/cpu/md/test_pe_multi_chain.cpp
tests/cpu/md/test_pe_nve.cpp
tests/cpu/md/test_pe_nvt.cpp
tests/cpu/md/test_pe_single_chain.cpp
tests/mpi/test_nose_hoover_t_ramp.cpp
tests/mpi/test_pe_chain_crossing.cpp
tests/mpi/test_pe_chain_drifting.cpp
tests/mpi/test_pe_chains_8rank.cpp
tests/mpi/test_pe_chains_multirank.cpp
tests/mpi/test_pe_full_multirank.cpp
tests/mpi/test_pe_msu_pipeline.cpp
tests/mpi/test_pe_nvt_multirank.cpp
benchmarks/md/reference/polyethylene/SOURCES.md        (resolved annotation)
.claude/specs/msu-pe-validation.md                      (status table only)
```

## Tasks

- [ ] Write a tiny RED test `tests/cpu/md/test_pe_ff_constants.cpp`
      that pin-points the Dreiding-UA values (asserts `kAngleK_kcal == 60.0`,
      `kAngleTheta0 ≈ 109.5° in radians`, `kDihedralC[0..3] == {1.73,
      −4.49, 0.776, 6.99}`). Confirm RED.
- [ ] Replace the three constant values in `benchmarks/md/polyethylene_tension.cpp`.
- [ ] Replace the three constant values in `benchmarks/md/pe_profile_single.cpp`.
- [ ] Replace the three constant values in each of the 13
      `tests/{cpu,mpi}/test_pe_*.cpp` + `test_nose_hoover_t_ramp.cpp`
      + `test_lj_exclusions.cpp` files. Confirm the new pin test is GREEN.
- [ ] Run full `ctest --test-dir build --output-on-failure` and verify
      no regression in the 89-test baseline (87 pre-spec + the post-
      pe-potential-abstraction `test_pair_lj_cut` + `test_pair_morse_smoke`).
      Total after this spec: **90 tests**.
- [ ] Update `benchmarks/md/reference/polyethylene/SOURCES.md`
      "Atomiverse FF discrepancy" subsection to "resolved" with the
      commit SHA.
- [ ] Update `.claude/specs/msu-pe-validation.md` status snapshot:
      flip Phase D from `🟡 blocked on FF migration` to `🟡 ready to run`.

## Testing strategy

**Happy path:**

- Every existing PE test passes unchanged (formula-based oracles adapt
  to the new constants).
- The new `test_pe_ff_constants.cpp` pin test guards against future
  drift back to OPLS values.
- Full ctest count: 89 → 90 (one new pin test).

**Edge cases:**

- `test_pe_nve.cpp`'s NVE drift bound (`< 1e-4` relative) was tuned on
  the old OPLS-like constants. A 3.2° shift in θ₀ stiffens the angle
  oscillation modestly. If the drift bound is now violated, the
  tolerance was over-tight on the wrong physics; document a one-line
  tolerance adjustment with rationale — do not paper over a real bug.
- `test_mpi_pe_msu_pipeline_n8.cpp` runs the full 5-stage protocol at 8
  ranks; this is the most sensitive integration test. Failure here =
  something deeper than a constant swap, abort the spec.

**Domain validation** (`$META.science.required: true`):

- After migration, an all-trans C₂₀H₄₂ single-chain energy
  decomposition under Atomiverse should agree with LAMMPS under
  `lammps_msu_paper.in` constants on the same geometry to **1e-6 eV
  absolute per term** for E_bond, E_angle, E_dihedral. This is
  `msu-pe-validation.md` Phase A.2's reference test; **this spec
  enables it but does not run it** (it's downstream).
- Quantitative LAMMPS-vs-Atomiverse paper-scale comparison
  (stress-strain modulus / yield / hardening within ±5%) is
  `msu-pe-validation.md` Phase D's job. This spec **unblocks** Phase D
  by making the physics match.

**Performance:**

- No kernel-timing change expected — same number of bonds / angles /
  dihedrals per step, only the coefficient values differ. Any drift
  larger than ±1% on `LJ::compute` / `BondHarmonic::compute` etc. is
  suspect and worth investigating.

## Out of scope

- **Editing the MSU reference files** (`lammps_msu_paper.in`,
  `msu_pe.data`). They are ground truth; ground truth flows MSU →
  Atomiverse, never the other way.
- **Paper-scale Atomiverse run on `msu_pe.data`** (10 chains × 1000
  monomers). That's `msu-pe-validation.md` Phase D's deliverable; this
  spec only unblocks it.
- **Tools to read `msu_pe.data` directly into a `cpu::Frame`** via the
  `molrs` LAMMPS data reader — tracked separately under
  `msu-pe-validation` Phase C.
- **Refactoring the duplicated FF constants into a shared header.**
  The 15 files all declare `kAngleK_kcal`, `kAngleTheta0`, etc. in
  their own scope. A single `include/atomiverse/benchmarks/pe_ff_constants.h`
  would prevent the next FF migration from being this painful.
  Tracked as a future `pe-ff-constants-extract` spec; not done here to
  keep this change mechanically simple and reviewable.
- **Architecture audit on the potential layer.** Already done as part
  of `pe-potential-abstraction` —
  `.claude/notes/potential-abstraction-audit.md` is the deliverable.

---
slug: gpw-edl-aimd
criteria:
  - id: ac-001
    summary: ExtXYZ loads 1214-atom Frame via molrs bridge with correct box and Z
    type: code
    pass_when: |
      cpu::load_extxyz on final_frame.xyz returns a Frame with 1214 atoms, box
      matrix = diag(13.128, 10.582, 47.054) Å within 1e-6, and per-atom type set
      to atomic number Z for all five elements {Zn,O,H,Na,Cl}. No hand-rolled C++
      parser exists (read goes through xyz_read_first_frame in molrs-cxxapi ->
      cpu::materialize).
    status: verified

  - id: ac-002
    summary: Missing Lattice header throws named exception, no fallback box
    type: code
    pass_when: |
      load_extxyz on an ExtXYZ file lacking a Lattice="..." header throws
      std::runtime_error naming the missing field; no default/identity box is
      silently substituted.
    status: verified

  - id: ac-003
    summary: In-tree GTH-PBE + DZVP-MOLOPT-SR-GTH tables for the five EDL elements
    type: code
    pass_when: |
      tests/cpu/gpw/test_cpu_gpw_data_table.cpp passes: the in-tree constexpr
      tables (cpu/gpw/basis_data/{gth_pbe,molopt_dzvp}.h) supply GTH-PBE-q1 (H),
      q6 (O), q12 (Zn, l=0,1,2 incl. the d channel), q9 (Na, semicore), q7 (Cl)
      and DZVP-MOLOPT-SR-GTH for all five elements. Every r_loc / C_k / r_l /
      h^l_ij and every contraction coefficient matches its published source
      (Krack TCA 114, 145 (2005); VandeVondele-Hutter JCP 127, 114105 (2007))
      to the digits printed there, with a provenance comment naming the source.
      Projector counts are derived from the table, not hardcoded at the call
      site. An unsupported element or channel THROWS (no fallback value).
    status: pending
    note: |
      REPLACES the original ac-003, which asserted the counts were "read from the
      CP2K GTH_POTENTIALS entry". Reading a third-party data file at test time is
      the banned dependency this refactor removes. Same numbers, in-tree source.

  - id: ac-004
    summary: No third-party file is read at build or test time; CP2K gate is gone
    type: code
    pass_when: |
      The ATV_CP2K_DATA_DIR CMake option and its tests/CMakeLists.txt gate are
      deleted; no test source references $ATV_CP2K_DATA_DIR, GTH_POTENTIALS or
      BASIS_MOLOPT as a path to read. A grep gate over tests/ + src/ + include/
      finds no third-party package path. The full GPW suite configures, builds and
      passes with NO external data directory present anywhere on the machine.
    status: pending
    note: |
      This is the refactor's binding acceptance signal. Before it, 14 of the GPW
      tests (forces, evaluator, aimd, geoopt, diis, heavy_elements + all parity)
      only registered when ATV_CP2K_DATA_DIR pointed at a CP2K install — without
      it GPW had 12 pure-math tests, zero real-molecule physics and zero force
      validation.

  - id: ac-005
    summary: Analytic GPW forces match the 6N+1 finite-difference oracle
    type: scientific
    pass_when: |
      Analytic gpw_forces agrees with the central 6N+1 FD oracle (step 1e-4 Bohr)
      component-wise to <= 5e-5 Ha/Bohr on H2 AND on ZnO (which exercises the l=2
      d-projector derivative, the l=3 f-basis derivative and the GGA grid term),
      running on the in-tree data tables.
    status: pending
    note: |
      Was `verified` (H2 2e-6, ZnO 2.7e-5 Ha/Bohr) — but that run read CP2K's data
      files. Reset to pending: it must be re-proven on the in-tree tables, which is
      also the sharpest check that those tables are transcribed correctly. This is
      the reference-free correctness gate for forces and needs no external engine.

  - id: ac-006
    summary: Single-point energies match pinned literature/reference anchors
    type: scientific
    pass_when: |
      benchmarks/qc/gpw/validation passes: H2, Ar, H2O, ZnO and NaCl GPW
      single-point total energies match their pinned constexpr anchor values
      within 1e-4 Ha (ZnO 0.10 mHa, NaCl 0.03 mHa, H2O 0.43 mHa were the measured
      residuals; the ~4e-4 Ha H2O gap is the known single-grid-vs-multigrid
      discretization residual). Anchors are hardcoded constants with a provenance
      comment — nothing is invoked at build or test time. E2E: benchmarks only,
      never ctest.
    status: pending
    note: |
      REPLACES the original ac-004 (CP2K single-point parity) AND the original
      ac-006 (CP2K force parity). The force-parity criterion is DELETED outright:
      the FD oracle (ac-005) is the term-complete correctness gate, so a CP2K force
      reference added a forbidden dependency and zero information.

  - id: ac-007
    summary: GPWEvaluator uses analytic forces (one SCF + one gradient pass)
    type: code
    pass_when: |
      GPWEvaluator reports capabilities() including Capability::Force and populates
      frame forces via gpw_forces (analytic) with exactly one SCF + one gradient
      pass per launch (no 6N+1 re-solve), on the in-tree data tables.
    status: pending
    note: Was verified; reset with ac-005 — same reason (proven only on CP2K data).

  - id: ac-008
    summary: FixCSVR reproduces the canonical kinetic-energy distribution
    type: scientific
    pass_when: |
      tests/cpu/md/test_cpu_md_csvr.cpp passes: on a cheap classical MD system at
      target T with τ=TIMECON, the sampled kinetic-energy distribution matches the
      canonical (Maxwell-Boltzmann) distribution within the test's statistical
      tolerance, and N_f correctly excludes COM degrees of freedom.
    status: verified
    note: CP2K-free (classical Ar); unaffected by the refactor.

  - id: ac-009
    summary: FixCSVR conserved-quantity diagnostic is drift-free
    type: runtime
    pass_when: |
      The CSVR effective-energy bookkeeping term H-tilde shows no secular drift
      over the test trajectory (drift below the test's tolerance), confirming Eq.15
      bookkeeping is accumulated correctly.
    status: verified
    note: CP2K-free; unaffected by the refactor.

  - id: ac-010
    summary: Capstone — 1214-atom EDL loads and runs a short GPW-CSVR NVT trajectory
    type: scientific
    pass_when: |
      benchmarks/qc/gpw/validation: final_frame.xyz loads into a cpu::Frame; a short
      GPW-AIMD NVT trajectory (PBE, five in-tree GTH-PBE PPs, DZVP-MOLOPT-SR-GTH,
      CSVR τ=20 fs, dt 0.5 fs, T 330 K, a few steps) launches and completes, AND the
      reference-free physics gates hold: analytic forces agree with the FD oracle on
      an EDL fragment, and the CSVR conserved quantity H-tilde shows no secular drift
      over the trajectory. E2E: benchmarks only, never ctest.
    status: pending
    note: |
      REPLACES the original ac-010, whose bar was per-step CP2K parity. The bar is now
      "launches + completes + passes the reference-free gates" — which is what actually
      proves the engine, and is reachable without an external engine. Blocked on ac-003
      /ac-004 (data) and on the M5 gth_matrix_element atom screening (currently O(N³):
      the unscreened atom loop is 81% of runtime at N=24 and projects to hundreds of
      hours/SCF at N=1214).

  - id: ac-011
    summary: Full-system wall-time and memory profile documented
    type: performance
    pass_when: |
      benchmarks/qc/gpw/validation produces a report recording per-step wall-time and
      peak memory for the full 1214-atom system after the gth_matrix_element atom
      screening + OpenMP parallelization of collocate/FFT/KS/force hot paths; the report
      exists and the run completes a single step within the node's memory.
    status: pending
    note: Depends on ac-010.

  - id: ac-012
    summary: New public symbols carry Doxygen with units
    type: docs
    pass_when: |
      cpu::load_extxyz, cpu::md::FixCSVR, the gpw_forces public entry points and the
      new in-tree data-table accessors each have @brief, @param (with eV/Å/fs units
      where applicable), and @return; the doc-coverage gate over the new/changed public
      headers passes.
    status: pending

  - id: ac-013
    summary: DFT-D3(BJ) energy and force (incl. CN-derivative) — reference-free FD gate
    type: scientific
    pass_when: |
      tests/cpu/gpw/test_cpu_gpw_dftd3.cpp passes: D3(BJ) dispersion energy uses the
      in-tree parameter table (cpu/gpw/internal/dftd3_params.h; PBE s6=1.0, a1=0.4289,
      s8=0.7875, a2=4.4407, cited to Grimme 2010 / Grimme-Ehrlich-Goerigk 2011) and the
      analytic force — INCLUDING the second-pass CN-derivative term — agrees with a
      central finite-difference force check component-wise within tolerance. Out of the
      AIMD launch critical path.
    status: pending
    note: |
      REPLACES the original ac-013, which required agreement with "CP2K's PBE-D3(BJ)
      numbers ... verified vs CP2K's dftd3.dat". The params are already in-tree and the
      FD check is what actually catches a missing CN-derivative term.
---

# Acceptance criteria — GPW EDL AIMD (refactored 2026-07-11)

## What changed and why

The original criteria could not be met **as written**, and not because the physics was
too hard. Four of them (ac-004, ac-006, ac-010, ac-013) demanded agreement with an
external CP2K run, and — more seriously — **every GPW test that touched a real molecule
read its pseudopotential and basis data out of a CP2K installation at test time**
(`$ATV_CP2K_DATA_DIR/{GTH_POTENTIALS,BASIS_MOLOPT}`, 14 tests behind that gate). With
the gate off, which is how CI must run, GPW had 12 pure-math tests, **no real-molecule
physics and no force validation at all**. That is a banned third-party dependency
sitting on the critical path of the whole domain.

The refactor keeps the physics and moves the provenance in-tree: GTH-PBE and
DZVP-MOLOPT-SR-GTH are **published tables** (Krack 2005; VandeVondele-Hutter 2007 — the
CP2K data files are transcriptions of those papers), so they become constexpr tables
under `cpu/gpw/basis_data/`, exactly as `cpu/scf/basis_data/pople_631g.h` already pins
its numbers. Previously-measured CP2K energies survive as **pinned anchor constants**
with provenance comments; nothing is invoked at build or test time. "Match CP2K exactly"
still binds — **the numbers are matched, the dependency is not.**

## The gates that actually prove the engine

Ranked by correctness value, and note that the top two need **no external engine at all**:

1. **ac-005 — analytic forces vs the 6N+1 FD oracle.** Term-complete and reference-free:
   a wrong Pulay, projector-derivative or GGA-grid term cannot survive it. This is why the
   deleted CP2K force-parity criterion carried zero information.
2. **ac-010 / ac-009 — energy conservation and drift-free CSVR H̃.** Reference-free
   trajectory-level physics.
3. **ac-003 / ac-004 — in-tree data + no third-party read.** The structural gate; ac-005
   doubles as the sharpest check that the tables were transcribed correctly.
4. **ac-006 — pinned single-point anchors.** Regression protection, not discovery.

ac-005 and ac-007 were reset from `verified` to `pending`: they were proven, but only on
CP2K-sourced data. They must be re-proven on the in-tree tables.

Placement follows the repo rule: `tests/` is unit-only (ms-scale, one function); every
converged-SCF / AIMD / anchor comparison lives in `benchmarks/qc/gpw/validation/` and runs
at acceptance time, never under ctest.

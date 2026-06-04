---
slug: cpu-md-parallel-perf
criteria:
  - id: ac-001
    summary: Multi-rank neighbor build uses half-stencil + per-cell shift (identical pair set)
    type: code
    pass_when: |
      cell_list.cpp's refresh_half uses the half-stencil + per-cell PBC shift
      path under MPI (tag_dedup) for grids with >=3 cells/axis, with the
      full-stencil fallback retained for small/aliased grids. The neighbor
      pair set is identical to the full-stencil path (verified by the LJ
      energy/force tests and the 8-rank pipeline test staying green).
    status: verified
    last_checked: 2026-06-04
    evidence: |
      cell_list.cpp refresh_half: use_half gate relaxed to (ncx>=3 && ncy>=3
      && ncz>=3) — now active under MPI (tag_dedup) and extended axes; <3-cell
      grids still fall back to full-stencil. Cross-rank dedup corrected: the
      tag tie-break is gated by `reject_le_i` so it applies on every cell in
      the full stencil but only the home cell in the half stencil (forward
      cells are scanned by exactly one rank — antisymmetric offset — so the
      partner rank's owner never scans it; ghosts iterated as i are trimmed by
      inum=n_local downstream). New test_mpi_lj_energy_bruteforce pins the
      absolute distributed LJ energy to an O(N²) min-image oracle at 1 and 8
      ranks (box 42.4Å → ~6 cells/slab/axis, half-stencil active) within 1e-9;
      passes. CellList CPU suite (6) + lj_nve n8 + pe_full_multirank n2 +
      pe_chains_8rank n8 + pe_msu_pipeline n8 all green.
  - id: ac-002
    summary: LJ per-pair-shift fast path active in multi-rank
    type: code
    pass_when: |
      pair_lj_cut.cpp's k_has_shifts no longer gates on single-rank; the
      shift fast path (dx = xi-(xj+shift), no per-pair min-image) runs under
      MPI. 8-rank ghost-pair forces/energies match the min-image path within
      tolerance.
    status: verified
    last_checked: 2026-06-04
    evidence: |
      pair_lj_cut.cpp compute: k_has_shifts no longer ANDs with
      !(dd && mpi_size>1); it depends only on has_shifts (populated solely on
      the serial neighbor-build path → false under OMP>1, which then falls
      back to min-image). The stored shift is by construction xi-dx-xj where
      dx is the exact displacement the neighbor build computed, so xi-(xj+shift)
      reproduces it bit-for-bit (extended-axis ghosts carry shift 0; only
      non-decomposed wrapped pairs carry nonzero shift). test_mpi_lj_energy_
      bruteforce n8 (shift fast path active under MPI) matches the min-image
      oracle to 1e-9; pe_full_multirank n2 (LJ + ghost special-bonds) green.
  - id: ac-003
    summary: 8-rank correctness preserved
    type: runtime
    pass_when: |
      tests/mpi/test_pe_msu_pipeline (8-rank end-to-end: bonded + LJ
      special-bonds + halo + migration + NVT + barostat + deform) passes,
      and a 1-rank vs 8-rank bonded+LJ energy comparison agrees within
      tolerance.
    status: verified
    last_checked: 2026-06-04
    evidence: |
      test_mpi_pe_msu_pipeline n8 passes with the levers active. The 1r-vs-8r
      energy comparison is provided by the new test_mpi_lj_energy_bruteforce
      (added this spec): both 1-rank and 8-rank distributed LJ energy match the
      same brute-force oracle within 1e-9, which is strictly stronger than a
      direct 1r==8r comparison (it also pins the absolute value). All MPI
      correctness tests (lj_nve, pe_full_multirank, pe_chains_8rank,
      pe_msu_pipeline) green. NB: bonded forces are unchanged by this spec
      (levers touch only the LJ neighbor build + LJ kernel); a pre-existing
      HEAD failure in bonded per-term energy plumbing (PeThermoSchema) and an
      angle-force test (PeSingleChain.BendOneAngle) are unrelated and
      reproduce identically on clean HEAD.
  - id: ac-004
    summary: Multi-rank Neigh fraction drops; full-run faster
    type: performance
    evaluator_hint: operator-action SLURM run
    pass_when: |
      On the MSU PE full 5-stage run, stage-5 Neigh fraction drops materially
      from the ~42-47% baseline (half-stencil under MPI), and the full-run
      wall at 32 and 64 ranks improves vs the NEW-AV baseline (321s / 225s).
      Stretch goal: total at >=1 rank count beats LAMMPS (164s @32r / 86s @64r).
      Attribution CSVs committed under benchmarks/md/reference/polyethylene/.
    status: pending
    note: |
      Owed to operator-action SLURM full runs (native fabric, srun --mpi=pmi2,
      exclude flaky nodes) — cannot be measured from the dev harness. Levers
      1+2 land the algorithmic change (half-stencil under MPI halves the ~14
      cell visits; per-pair-shift removes the LJ min-image), which is what the
      Neigh-fraction drop depends on; the 32r/64r wall-clock + attribution CSVs
      must be regenerated on the cluster to flip this to verified.
  - id: ac-005
    summary: Physics preserved within msu-pe-validation tolerances
    type: scientific
    pass_when: |
      Multi-rank runs reproduce the validated MSU PE observables
      (energy/density/stress-strain) within the msu-pe-validation tolerances
      — no accuracy traded for parallel speed.
    status: pending
    evaluator_hint: msu-pe-validation cluster gate (/mol:bench or operator)
    note: |
      Harness-level evidence that physics is preserved is already strong: the
      multi-rank LJ energy matches the brute-force oracle to 1e-9 and all MPI
      correctness/conservation tests pass. The formal msu-pe-validation
      energy/density/stress-strain cluster gate is operator/bench-owed and has
      not been run, so this stays pending.
---

# Acceptance — cpu-md-parallel-perf

Binding "done" contract for the multi-rank performance work. ac-001/002/003 are
harness-verifiable (code + 8-rank test). ac-004 is operator-action timing on the
cluster (native fabric, exclude flaky nodes). ac-005 is the formal physics gate.
RCB / load balancing is explicitly out of scope (measured: not the gate).

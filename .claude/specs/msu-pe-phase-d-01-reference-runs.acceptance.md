---
slug: msu-pe-phase-d-01-reference-runs
criteria:
  - id: ac-001
    summary: deck-fidelity audit asserts lammps.in physics equals lammps_msu_paper.in
    type: code
    pass_when: |
      tools/audit_lammps_deck.sh exists and, run against the in-tree
      benchmarks/md/reference/polyethylene/lammps.in and
      lammps_msu_paper.in, exits 0 — every fix, ensemble, integrator,
      timestep, run step-count and force-field coefficient matches,
      with only the non-physics output hooks (fix print / dump /
      compute rdf / log / write_restart / reset_timestep / thermo_style)
      differing; a deliberately mutated lammps.in (changed timestep)
      makes it exit non-zero naming the offending line.
    status: verified
    last_checked: 2026-05-22
  - id: ac-002
    summary: canonicalize_msu_output.py converts raw engine output to the canonical schema
    type: code
    evaluator_hint: pytest tools/test_canonicalize_msu_output.py
    pass_when: |
      tools/canonicalize_msu_output.py exists with engine ∈ {lammps,
      atomiverse} and kind ∈ {thermo,stress-strain} subcommands; the
      pytest suite tools/test_canonicalize_msu_output.py passes,
      covering thermo→`step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz`,
      stress-strain→`strain,sigma_xx`, the space-separated def1/def2
      LAMMPS format, header-only input, and a malformed-row loud failure.
    status: verified
    last_checked: 2026-05-22
  - id: ac-003
    summary: stress_strain_lammps_8r.csv generated from the existing 8-rank paper run
    type: code
    pass_when: |
      benchmarks/md/reference/polyethylene/stress_strain_lammps_8r.csv
      exists, has header `strain,sigma_xx`, was produced by
      canonicalize_msu_output.py from
      runs/lammps_n8_paper/PE_nc10_cl1000.def1.txt; its first strain
      row is ≈ 0, the strain column increases monotonically toward
      ε ≈ 1.0, and sigma_xx is in GPa (tension-positive).
    status: verified
    last_checked: 2026-05-22
  - id: ac-004
    summary: run_msu_lammps.sh uses lmp_cuda CPU mode and shared canonicalization
    type: code
    pass_when: |
      tools/run_msu_lammps.sh defaults its LAMMPS binary to lmp_cuda
      run in CPU mode (no GPU package args), keeps LMP_BIN override,
      and routes raw output to the canonical CSVs via
      canonicalize_msu_output.py rather than an inline copy; grep of
      the script confirms lmp_cuda and the canonicalize call.
    status: verified
    last_checked: 2026-05-22
  - id: ac-005
    summary: run_msu_pdc.sh submits PDC Dardel main-partition sbatch jobs
    type: code
    pass_when: |
      tools/run_msu_pdc.sh exists, accepts an engine (lammps|atomiverse)
      and rank count, requires an -A <project> account flag (exits
      non-zero with a clear error when omitted, no default account),
      targets the `main` partition with a 24h wall-time cap, and
      submits via sbatch wrapping run_msu_lammps.sh / run_msu_atomiverse.sh.
    status: verified
    last_checked: 2026-05-22
  - id: ac-006
    summary: run scripts pass a shortened-protocol smoke test producing schema-correct CSVs
    type: runtime
    pass_when: |
      A shortened-protocol smoke run of run_msu_lammps.sh and
      run_msu_atomiverse.sh (reduced run step counts via the extra-args
      path, on a login node) completes and produces thermo_<engine>_<N>r.csv
      with header `step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz` and
      stress_strain_<engine>_<N>r.csv with header `strain,sigma_xx`.
    status: pending
  - id: ac-007
    summary: SOURCES.md and README.md match the finalized scripts and step count
    type: docs
    pass_when: |
      SOURCES.md § Phase D.1 and README.md state the lmp_cuda CPU-mode
      launch, the PDC Dardel main-partition sbatch procedure, the
      331800-step total (correcting the prior ~360000 figure), and the
      canonicalize_msu_output.py canonicalization step.
    status: pending
  - id: ac-008
    summary: four full PDC reference runs land the eight canonical CSVs and stage-4 artifacts
    type: runtime
    pass_when: |
      After the operator runs run_msu_pdc.sh for lammps and atomiverse
      at 1 and 8 ranks on PDC Dardel, the eight canonical CSVs
      (thermo_{lammps,atomiverse}_{1,8}r.csv,
      stress_strain_{lammps,atomiverse}_{1,8}r.csv) plus the stage-4
      thermo log, finalized RDF curve and trajectory XYZ for both
      engines at both rank counts are committed under
      benchmarks/md/reference/polyethylene/.
    status: pending
---

# Acceptance criteria

- **ac-001** (code) — Turns the `SOURCES.md` § Phase D.1 prose claim ("`lammps.in` physics ≡ `lammps_msu_paper.in`") into a machine-verifiable gate. Verified by running the audit script on the in-tree decks and on a mutated copy.
- **ac-002** (code) — The canonicalization converter is pure tooling with deterministic column mappings; its pytest suite verifies every engine/kind combination, the space-separated `def1`/`def2` LAMMPS format, and edge cases independently of any cluster run.
- **ac-003** (code) — A real, already-completed 8-rank LAMMPS paper-deck run is brought into the canonical schema *now*, giving `msu-pe-phase-d-04`'s ac-002 load-bearing gate its first real LAMMPS-side reference data without a re-run.
- **ac-004** (code) — The LAMMPS launcher is finalized to the user's `lmp_cuda` binary in CPU mode and to the shared canonicalization path. Grep-verifiable.
- **ac-005** (code) — The PDC Dardel submission wrapper exists with a mandatory account flag (no fallback per project convention) and `main`-partition / 24h-cap awareness. Grep-verifiable on the script.
- **ac-006** (runtime) — A shortened-protocol smoke run proves the run scripts wire end-to-end and emit correctly-schema'd CSVs before any 24h cluster commitment. Runs on a login node, not CI.
- **ac-007** (docs) — The reference-directory docs must not drift from the finalized scripts; the prior ~360000-step figure is corrected to the actual 331800.
- **ac-008** (runtime) — The deferred operator gate, parked at code-complete (same shape as `msu-pe-phase-d-04` ac-002). The four ~24h PDC Dardel reference runs and the commit of the eight canonical CSVs plus the stage-4 thermo / RDF / trajectory artifacts. Not a CI step — an operator action that unblocks the downstream `msu-pe-align-03` / `d-03` / `d-04` runtime gates.

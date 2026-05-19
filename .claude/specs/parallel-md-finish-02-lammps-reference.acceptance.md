---
slug: parallel-md-finish-02-lammps-reference
criteria:
  - id: ac-001
    summary: lammps.in is self-contained (no external data file required)
    type: code
    pass_when: |
      `grep -E '^\s*read_data' benchmarks/md/reference/polyethylene/lammps.in`
      returns zero matches. `grep -E '^\s*(lattice|create_box|create_atoms|create_bonds)'`
      returns ≥ 4 matches. Running
      `/opt/homebrew/bin/lmp_mpi -in benchmarks/md/reference/polyethylene/lammps.in
      -log /tmp/lammps_smoke.log` from the directory exits 0.
    status: pending

  - id: ac-002
    summary: thermo_reference.csv committed with documented header and data
    type: code
    pass_when: |
      `benchmarks/md/reference/polyethylene/thermo_reference.csv`
      exists, first line is exactly
      `step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz`, and the file
      has ≥ 11 lines (header + ≥ 10 data rows). The last
      data row's `T` value is finite and within [80, 150]
      (end of stage 4 cooling).
    status: pending

  - id: ac-003
    summary: stress_strain_reference.csv committed with documented header and data
    type: code
    pass_when: |
      `benchmarks/md/reference/polyethylene/stress_strain_reference.csv`
      exists, first line is exactly `strain,sigma_xx`, and
      the file has ≥ 11 lines (header + ≥ 10 data rows).
      The last data row's `strain` value is ≥ 0.3
      (proving the tension stage actually ran).
    status: pending

  - id: ac-004
    summary: README.md version + hardware blocks filled with real values
    type: docs
    pass_when: |
      `benchmarks/md/reference/polyethylene/README.md` contains
      no `<TBD` markers in the "Required versions" subsection.
      `grep -E 'lammps[- _]?(stable|patch)'` returns ≥ 1 match
      (regex from parent ac-005). The README has a subsection
      whose heading matches `(?i)hardware|cpu|fingerprint` and
      contains the substring of `sysctl -n machdep.cpu.brand_string`
      output. (Closes parent ac-005 README portion.)
    status: pending
---

# Acceptance criteria

- **ac-001** is the load-bearing change: by removing the `read_data` dependency, the LAMMPS pipeline becomes runnable from a clean checkout without any moltemplate / MSU-wiki download. The smoke-run requirement proves the procedural geometry block actually compiles and runs through all 5 stages.
- **ac-002 / ac-003** verify that the two CSVs needed by parent ac-005 land with the documented column schema. End-of-stage-4 T ∈ [80, 150] K confirms NVT cooling ran; final strain ≥ 0.3 confirms tension actually deformed the box.
- **ac-004** closes the version + hardware annotation requirement from parent ac-005. The `lammps[-_ ]?(stable|patch)` regex is inherited verbatim from the parent.

Together these four close parent ac-005 in full.

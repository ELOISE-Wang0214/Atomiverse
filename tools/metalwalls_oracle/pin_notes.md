# Pin notes — metalwalls-method-loop

## four-site-tip4p-foothold (2026-07-23)

- **Oracle:** LAMMPS `lj/cut/tip4p/cut` (units metal)
- **Binary:** `/home/jicli594/.local/x86_64/bin/lmp` — 30 Mar 2026 Development
- **Inputs:** `tools/metalwalls_oracle/inputs/lammps/four-site-tip4p-foothold/`
- **Command:** `lmp -in in.tip4p_dimer_metal -log log.metal`
- **Pin:** PE = 0.52924228894955327451 eV; forces in REFERENCE.md
- **Conversion:** metal units match ATV eV; ε_OO = 0.1550 kcal/mol × `kKcalPerMolInEv`

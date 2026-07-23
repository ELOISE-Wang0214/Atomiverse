# Pin notes — metalwalls-method-loop

## four-site-tip4p-foothold (2026-07-23)

- **Oracle:** LAMMPS `lj/cut/tip4p/cut` (units metal)
- **Binary:** `/home/jicli594/.local/x86_64/bin/lmp` — 30 Mar 2026 Development
- **Inputs:** `tools/metalwalls_oracle/inputs/lammps/four-site-tip4p-foothold/`
- **Command:** `lmp -in in.tip4p_dimer_metal -log log.metal`
- **Pin:** PE = 0.52924228894955327451 eV; forces in REFERENCE.md
- **Conversion:** metal units match ATV eV; ε_OO = 0.1550 kcal/mol × `kKcalPerMolInEv`

## born-repulsion-foothold (2026-07-24)

- **Oracle:** LAMMPS `pair_style born` (units metal, shift yes)
- **Binary:** `/home/jicli594/.local/x86_64/bin/lmp` — 30 Mar 2026 Development
- **Inputs:** `tools/metalwalls_oracle/inputs/lammps/in.born_dimer`
- **Command:** `lmp -in in.born_dimer`
- **Pin:** PE = 509.3671885909988 eV; F1x = -201.72730961087780 eV/Å
- **Formula:** A exp((σ−r)/ρ) − C/r⁶ + D/r⁸ (ρ softness length)
## fumi-tosi-nacl (2026-07-24)

- **Oracle:** LAMMPS `hybrid/overlay born + coul/cut` (units metal; `pair_modify pair born shift yes`)
- **Inputs:** `tools/metalwalls_oracle/inputs/lammps/in.fumi_tosi_nacl4`
- **Command:** `lmp -in in.fumi_tosi_nacl4`
- **Pin:** PE = 54.7707332302457 eV; forces in REFERENCE.md
- **ATV:** multi-type `cpu::PairBorn` + `cpu::PairCoulCut` on 4-ion NaCl square

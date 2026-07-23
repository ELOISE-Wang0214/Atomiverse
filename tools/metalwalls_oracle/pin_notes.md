# Pin notes — metalwalls-method-loop

## four-site-tip4p-foothold (2026-07-23)

- **Oracle:** LAMMPS `lj/cut/tip4p/cut` (units metal)
- **Binary:** `/home/jicli594/.local/x86_64/bin/lmp` — 30 Mar 2026 Development
- **Inputs:** `tools/metalwalls_oracle/inputs/lammps/four-site-tip4p-foothold/`
- **Command:** `lmp -in in.tip4p_dimer_metal -log log.metal`
- **Pin:** PE = 0.52924228894955327451 eV; forces in REFERENCE.md
- **Conversion:** metal units match ATV eV; ε_OO = 0.1550 kcal/mol × `kKcalPerMolInEv`

## thomas-fermi-electrode (2026-07-24)

- **Oracle:** MetalWalls Thomas–Fermi model (Scalfi et al. JCP 153, 174704 (2020))
- **Why not LAMMPS:** ELECTRODE package is perfect conductor only (no \(l_{\mathrm{TF}}\))
- **Binary:** no `mw` on agent host → offline discrete pin (same algebra as ConstPotBO + TF diagonal)
- **Script:** `tools/metalwalls_oracle/inputs/metalwalls/thomas-fermi-electrode/pin_gold_cap.py`
- **Command:** `python3 .../pin_gold_cap.py`
- **Pin:** empty Au capacitor toy, \(\Delta\Psi=1\) V, \(\eta=1/0.56\) Å⁻¹, \(d=4/a_0^3\),
  \(R=40\): \(Q_{\mathrm{top}}(l_{\mathrm{TF}}=0)=4.620945546104\times 10^{-3}\,e\),
  \(Q_{\mathrm{top}}(l_{\mathrm{TF}}=2\,\mathrm{Å})=3.934601758679\times 10^{-3}\,e\)
- **perf_exempt:** no MetalWalls twin binary

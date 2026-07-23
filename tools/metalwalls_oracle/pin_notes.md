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
## dipole-electrode-joint (2026-07-24)

- **Oracle:** MetalWalls multipole + Gaussian electrode BO form (JCP 2022 / JOSS 2020)
- **Why not LAMMPS:** no SCF PIM ⊕ const-pot joint twin (ELECTRODE non-polarizable; Drude ≠ PIM)
- **Binary:** no `mw` on agent host → offline algebraic pin
- **Script:** `tools/metalwalls_oracle/inputs/metalwalls/dipole-electrode-joint/pin_joint.py`
- **Command:** `python3 .../pin_joint.py`
- **Pin (joint α>0):** Q_top = +7.3285008958883610e-02 e; μ_cat,z = 5.8951916148990660e-02 e·Å; μ_an,z = 8.7432701391015569e-02 e·Å
- **Pin (α=0):** Q_top = +7.5339627935767609e-02 e (matches ConstPotBO)
- **perf_exempt:** no MetalWalls twin binary

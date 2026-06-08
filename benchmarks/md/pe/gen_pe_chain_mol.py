#!/usr/bin/env python3
"""Generate a LAMMPS molecule template for ONE united-atom PE zig-zag chain,
matching benchmarks/md/polyethylene_tension.cpp make_pe_melt_local() exactly,
plus print the lattice geometry so a deck can `create_atoms mol` an nx*ny*nz
melt with NO giant data file (scalable to any size).

Chain geometry (identical to make_pe_melt_local):
  R0=1.53, theta0=109.5deg, mass(CH2)=14.027
  dz = R0*sqrt((1+cos theta0)/2),  dx = sqrt(R0^2 - dz^2)
  atom m (0..L-1): x = (+/-) dx/2 alternating, y = 0, z = m*dz
  per-cell box: spacing x spacing x cz, cz = (L-1)*dz + gap_z
Topology: linear chain -> (L-1) bonds, (L-2) angles, (L-3) dihedrals, all type 1.
"""
import math, argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, default=18)        # atoms per chain
    ap.add_argument("--spacing", type=float, default=5.0)
    ap.add_argument("--gap-z", type=float, default=4.0)
    ap.add_argument("--out", default="pe_chain.mol")
    a = ap.parse_args()
    L = a.L
    theta0 = math.radians(109.5)
    R0 = 1.53
    dz = R0 * math.sqrt((1.0 + math.cos(theta0)) / 2.0)
    dx = math.sqrt(R0 * R0 - dz * dz)
    cz = (L - 1) * dz + a.gap_z

    coords = []
    for m in range(L):
        sgn = 1.0 if (m % 2 == 0) else -1.0
        coords.append((sgn * dx / 2.0, 0.0, m * dz))

    # One chain CENTERED in a (spacing x spacing x cz) unit cell, so that
    # `read_data` + `replicate nx ny nz` reproduces make_pe_melt_local's
    # aligned grid exactly — no random rotation, no clashes.
    sp = a.spacing
    cx = cy = sp / 2.0
    cz0 = a.gap_z / 2.0
    placed = [(cx + x, cy + y, cz0 + z) for (x, y, z) in coords]

    lines = []
    lines.append("# United-atom PE zig-zag chain unit cell "
                 "(matches make_pe_melt_local); replicate to build a melt")
    lines.append("")
    lines.append(f"{L} atoms")
    lines.append(f"{L-1} bonds")
    lines.append(f"{L-2} angles")
    lines.append(f"{L-3} dihedrals")
    lines.append("")
    lines.append("1 atom types")
    lines.append("1 bond types")
    lines.append("1 angle types")
    lines.append("1 dihedral types")
    lines.append("")
    lines.append(f"0.0 {sp:.10f} xlo xhi")
    lines.append(f"0.0 {sp:.10f} ylo yhi")
    lines.append(f"0.0 {cz:.10f} zlo zhi")
    lines.append("")
    lines.append("Masses")
    lines.append("")
    lines.append("1 14.027")
    lines.append("")
    lines.append("Atoms # molecular")
    lines.append("")
    for i,(x,y,z) in enumerate(placed, 1):
        # atom-ID molecule-ID type x y z
        lines.append(f"{i} 1 1 {x:.10f} {y:.10f} {z:.10f}")
    lines.append("")
    lines.append("Bonds")
    lines.append("")
    for b in range(1, L):
        lines.append(f"{b} 1 {b} {b+1}")
    lines.append("")
    lines.append("Angles")
    lines.append("")
    for g in range(1, L-1):
        lines.append(f"{g} 1 {g} {g+1} {g+2}")
    lines.append("")
    lines.append("Dihedrals")
    lines.append("")
    for d in range(1, L-2):
        lines.append(f"{d} 1 {d} {d+1} {d+2} {d+3}")
    lines.append("")
    with open(a.out, "w") as f:
        f.write("\n".join(lines))

    print(f"wrote {a.out}: L={L} bonds={L-1} angles={L-2} dihedrals={L-3}")
    print(f"GEOM dz={dz:.10f} dx={dx:.10f} cz={cz:.10f} spacing={a.spacing}")
    print(f"LATTICE a1={a.spacing} a2={a.spacing} a3={cz:.10f}  "
          f"(box = nx*{a.spacing} x ny*{a.spacing} x nz*{cz:.6f} Angstrom)")
    # density check: mass per chain / cell volume
    vol = a.spacing * a.spacing * cz
    m_chain = L * 14.027  # g/mol
    NA = 6.02214076e23
    rho = (m_chain / NA) / (vol * 1e-24)  # g/cm^3
    print(f"DENSITY ~= {rho:.4f} g/cc  (cell vol {vol:.2f} A^3, {L} UA/chain)")

if __name__ == "__main__":
    main()

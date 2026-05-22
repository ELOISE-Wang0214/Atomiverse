#!/usr/bin/env python3
"""Overlay the Atomiverse and LAMMPS MSU-PE uniaxial-tension curves.

Stdlib-only SVG (reuses the plotting helpers in plot_msu_reference.py).
Reads the canonical stress-strain reference CSVs and writes:
  msu_compare_stress_strain.svg

Usage: plot_msu_compare.py [reference_dir]
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_msu_reference import _panel, _svg, _read_csv, _roll, _pct


def _linfit(xs, ys):
    n = len(xs)
    if n < 2:
        return 0.0
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    den = n * sxx - sx * sx
    return (n * sxy - sx * sy) / den if den else 0.0


def _descriptors(eps, sig):
    """Elastic modulus (ε<0.02 least-squares slope) and yield (max σ, ε≤0.1)."""
    el = [(x, y) for x, y in zip(eps, sig) if 0.0 <= x < 0.02]
    E = _linfit([x for x, _ in el], [y for _, y in el])
    yld = max((y for x, y in zip(eps, sig) if x <= 0.1), default=0.0)
    return E, yld


def main(argv):
    ref = argv[1] if len(argv) > 1 else os.path.join(
        "benchmarks", "md", "reference", "polyethylene")
    av = _read_csv(os.path.join(ref, "stress_strain_atomiverse_1r.csv"))
    lm = _read_csv(os.path.join(ref, "stress_strain_lammps_1r.csv"))
    ae, asig = av["strain"], av["sigma_xx"]
    le, lsig = lm["strain"], lm["sigma_xx"]

    aE, aY = _descriptors(ae, asig)
    lE, lY = _descriptors(le, lsig)

    width, height = 870, 590
    ox, oy, pw, ph = 86.0, 66.0, 706.0, 446.0
    xr = (0.0, max(max(ae), max(le)))
    allsig = sorted(asig + lsig)
    lo, hi = _pct(allsig, 0.003), _pct(allsig, 0.997)
    pad = 0.09 * ((hi - lo) or 1.0)

    C_LMP, C_LMP_R = "#1d3557", "#86a8c4"
    C_AV, C_AV_R = "#c1121f", "#e6a99f"
    series = [
        (le, lsig, C_LMP_R, 0.7, 0.30),
        (ae, asig, C_AV_R, 0.7, 0.30),
        (le, _roll(lsig, 21), C_LMP, 2.4, 1.0),
        (ae, _roll(asig, 21), C_AV, 2.4, 1.0),
    ]
    clip, frag = _panel(
        0, ox, oy, pw, ph, series, xr, (lo - pad, hi + pad),
        "engineering strain  ε = Lₓ/Lₓ₀ − 1",
        "tensile stress  σₓₓ  [GPa]",
        "MSU amorphous-PE uniaxial tension — Atomiverse vs LAMMPS")

    head = ('<text x="%.1f" y="26" font-size="16.5" font-weight="bold" '
            'text-anchor="middle">MSU amorphous-PE — our code (Atomiverse) '
            'vs LAMMPS reference</text>' % (width / 2))

    # legend + descriptor box, upper-left inside the panel
    lx, ly = ox + 26, oy + 26
    leg = ['<g font-size="12">']
    leg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
               'stroke-width="2.8"/>' % (lx, ly, lx + 30, ly, C_LMP))
    leg.append('<text x="%.1f" y="%.1f" fill="#222">LAMMPS  (lmp_cuda, '
               'reference)</text>' % (lx + 38, ly + 4))
    leg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
               'stroke-width="2.8"/>' % (lx, ly + 22, lx + 30, ly + 22, C_AV))
    leg.append('<text x="%.1f" y="%.1f" fill="#222">Atomiverse  (CPU, '
               'this run)</text>' % (lx + 38, ly + 26))
    leg.append('<text x="%.1f" y="%.1f" font-size="11.5" fill="#555">'
               'elastic modulus E:  LAMMPS %.2f  ·  Atomiverse %.2f GPa  '
               '(%.1f×)</text>' % (lx, ly + 50, lE, aE,
                                   aE / lE if lE else 0.0))
    leg.append('<text x="%.1f" y="%.1f" font-size="11.5" fill="#555">'
               'yield σ_y:  LAMMPS %.3f  ·  Atomiverse %.3f GPa</text>'
               % (lx, ly + 67, lY, aY))
    leg.append('</g>')

    body = "<defs>%s</defs>\n%s\n%s\n%s" % (clip, head, frag, "\n".join(leg))
    out = os.path.join(ref, "msu_compare_stress_strain.svg")
    with open(out, "w") as fh:
        fh.write(_svg(width, height, body))
    print("wrote", out)
    print("  LAMMPS     E=%.3f GPa  yield=%.3f GPa" % (lE, lY))
    print("  Atomiverse E=%.3f GPa  yield=%.3f GPa" % (aE, aY))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

#!/usr/bin/env python3
"""Plot the MSU-PE LAMMPS reference run as standalone SVG figures.

Stdlib-only (csv, math) — no matplotlib / numpy, so it is safe under the
PDC home-directory inode quota and needs no environment setup.

Reads the canonical reference CSVs under a reference directory and writes:
  msu_lammps_thermo.svg         T / P / E_pot over the 5-stage protocol
  msu_lammps_stress_strain.svg  sigma_xx vs engineering strain (stage 5)

Usage: plot_msu_reference.py [reference_dir]
       (default reference_dir: benchmarks/md/reference/polyethylene)
"""
import csv
import math
import os
import sys

# 5-stage protocol boundaries in the monotonic step column (lammps.in omits
# reset_timestep): 10000 + 50000 + 50000 + 50000 + 171800 = 331800.
STAGE_EDGES = [10000, 60000, 110000, 160000]
STAGE_SPANS = [(0, 10000), (10000, 60000), (60000, 110000),
               (110000, 160000), (160000, 331800)]
STAGE_TAGS = ["S1", "S2", "S3", "S4", "S5"]


# ----------------------------------------------------------------------
# small numeric helpers
# ----------------------------------------------------------------------
def _nice_step(rough):
    if rough <= 0:
        return 1.0
    mag = 10.0 ** math.floor(math.log10(rough))
    for m in (1.0, 2.0, 2.5, 5.0, 10.0):
        if rough <= m * mag:
            return m * mag
    return 10.0 * mag


def _axis(lo, hi, target=6):
    if not (hi > lo):
        hi = lo + 1.0
    step = _nice_step((hi - lo) / target)
    a = math.floor(lo / step) * step
    b = math.ceil(hi / step) * step
    ticks, v = [], a
    while v <= b + step * 1e-6:
        ticks.append(v)
        v += step
    return a, b, ticks


def _pct(sorted_vals, q):
    if not sorted_vals:
        return 0.0
    i = int(round(q * (len(sorted_vals) - 1)))
    return sorted_vals[min(len(sorted_vals) - 1, max(0, i))]


def _roll(ys, w):
    """Centred moving average, window w; window shrinks at the edges."""
    n = len(ys)
    if w < 2 or n == 0:
        return list(ys)
    half = w // 2
    out = []
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        out.append(sum(ys[lo:hi]) / (hi - lo))
    return out


def _fmt(v):
    if v == 0:
        return "0"
    av = abs(v)
    if av >= 1e5 or av < 1e-3:
        return "%.0e" % v
    if av >= 100:
        return "%.0f" % v
    if av >= 1:
        return ("%.2f" % v).rstrip("0").rstrip(".")
    return ("%.3f" % v).rstrip("0").rstrip(".")


def _fmt_k(v):
    if v == 0:
        return "0"
    return "%gk" % (v / 1000.0)


# ----------------------------------------------------------------------
# SVG panel
# ----------------------------------------------------------------------
def _panel(idx, ox, oy, w, h, series, xr, yr, xlabel, ylabel, title,
           xfmt=_fmt, stages=False):
    """SVG fragment for one plot panel.

    series : list of (xs, ys, colour, stroke_width, opacity)
    xr, yr : (lo, hi) data ranges for the axes
    stages : draw the 5-stage boundary guides
    """
    xa, xb, xticks = _axis(*xr)
    ya, yb, yticks = _axis(*yr)
    sw = (xb - xa) or 1.0
    sh = (yb - ya) or 1.0

    def sx(x):
        return ox + (x - xa) / sw * w

    def sy(y):
        return oy + h - (y - ya) / sh * h

    o = []
    o.append('<text x="%.1f" y="%.1f" font-size="15" font-weight="bold" '
             'text-anchor="middle">%s</text>' % (ox + w / 2, oy - 9, title))
    # y grid + ticks
    for t in yticks:
        y = sy(t)
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                 'stroke="#e6e6e6" stroke-width="1"/>'
                 % (ox, y, ox + w, y))
        o.append('<text x="%.1f" y="%.1f" font-size="11" text-anchor="end" '
                 'fill="#444">%s</text>' % (ox - 8, y + 4, _fmt(t)))
    # x ticks
    for t in xticks:
        x = sx(t)
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                 'stroke="#333" stroke-width="1"/>' % (x, oy + h, x, oy + h + 5))
        o.append('<text x="%.1f" y="%.1f" font-size="11" '
                 'text-anchor="middle" fill="#444">%s</text>'
                 % (x, oy + h + 18, xfmt(t)))
    # stage boundary guides
    if stages:
        for e in STAGE_EDGES:
            if xa <= e <= xb:
                x = sx(e)
                o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                         'stroke="#b9b9b9" stroke-width="1" '
                         'stroke-dasharray="4 3"/>' % (x, oy, x, oy + h))
        for (s0, s1), tag in zip(STAGE_SPANS, STAGE_TAGS):
            c = 0.5 * (max(s0, xa) + min(s1, xb))
            if xa <= c <= xb:
                o.append('<text x="%.1f" y="%.1f" font-size="10" '
                         'text-anchor="middle" fill="#999">%s</text>'
                         % (sx(c), oy + 13, tag))
    # data (clipped to the panel)
    o.append('<g clip-path="url(#clip%d)">' % idx)
    for xs, ys, colour, swidth, opac in series:
        pts = " ".join("%.1f,%.1f" % (sx(x), sy(y)) for x, y in zip(xs, ys))
        o.append('<polyline points="%s" fill="none" stroke="%s" '
                 'stroke-width="%s" stroke-opacity="%s" '
                 'stroke-linejoin="round"/>' % (pts, colour, swidth, opac))
    o.append('</g>')
    # frame
    o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" '
             'stroke="#333" stroke-width="1.3"/>' % (ox, oy, w, h))
    # axis labels
    o.append('<text x="%.1f" y="%.1f" font-size="12.5" '
             'text-anchor="middle">%s</text>'
             % (ox + w / 2, oy + h + 38, xlabel))
    o.append('<text x="%.1f" y="%.1f" font-size="12.5" text-anchor="middle" '
             'transform="rotate(-90 %.1f %.1f)">%s</text>'
             % (ox - 48, oy + h / 2, ox - 48, oy + h / 2, ylabel))
    clip = ('<clipPath id="clip%d"><rect x="%.1f" y="%.1f" width="%.1f" '
            'height="%.1f"/></clipPath>' % (idx, ox, oy, w, h))
    return clip, "\n".join(o)


def _svg(width, height, body):
    return ('<?xml version="1.0" encoding="UTF-8"?>\n'
            '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
            'viewBox="0 0 %d %d" font-family="DejaVu Sans, Arial, sans-serif">\n'
            '<rect width="%d" height="%d" fill="white"/>\n%s\n</svg>\n'
            % (width, height, width, height, width, height, body))


def _read_csv(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    head = rows[0]
    cols = {name: [] for name in head}
    for r in rows[1:]:
        if len(r) != len(head):
            continue
        for name, cell in zip(head, r):
            try:
                cols[name].append(float(cell))
            except ValueError:
                cols[name].append(float("nan"))
    return cols


# ----------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------
def plot_thermo(cols, out_path):
    step = cols["step"]
    pw, ph = 820.0, 175.0
    ox, gap, top = 70.0, 52.0, 36.0
    width, height = 920, int(top + 3 * ph + 2 * gap + 52)
    xr = (min(step), max(step))

    specs = [
        ("T", "Temperature  [K]", "#d1495b", "Temperature", False),
        ("P", "Pressure  [atm]", "#3a7ca5", "Pressure", True),
        ("E_pot", "E_pot  [kcal/mol]", "#2a9d8f", "Potential energy", True),
    ]
    clips, frags = [], []
    for i, (key, ylabel, colour, title, smooth) in enumerate(specs):
        ys = cols[key]
        oy = top + i * (ph + gap)
        # Frame the y-axis on stages 2-5 (controlled thermodynamics).
        # Stage 1 is an overlap-relaxation transient (T spikes to ~810 K,
        # E_pot starts near 1e8) that would otherwise squash the panel;
        # it clips cleanly at the panel edge instead.
        eq = sorted(y for s, y in zip(step, ys) if s > 10000 and y == y)
        if not eq:
            eq = sorted(v for v in ys if v == v)
        lo, hi = _pct(eq, 0.005), _pct(eq, 0.995)
        pad = 0.08 * ((hi - lo) or 1.0)
        series = [(step, ys, colour, 0.7, 0.35 if smooth else 0.9)]
        if smooth:
            series.append((step, _roll(ys, 25), colour, 1.8, 1.0))
        last = (i == len(specs) - 1)
        clip, frag = _panel(
            i, ox, oy, pw, ph, series, xr, (lo - pad, hi + pad),
            "MD step" if last else "", ylabel, title,
            xfmt=_fmt_k, stages=True)
        clips.append(clip)
        frags.append(frag)

    head = ('<text x="%.1f" y="22" font-size="16.5" font-weight="bold" '
            'text-anchor="middle">MSU amorphous-PE — LAMMPS reference '
            'thermodynamics (5-stage protocol)</text>' % (width / 2))
    body = "<defs>%s</defs>\n%s\n%s" % ("".join(clips), head,
                                        "\n".join(frags))
    with open(out_path, "w") as fh:
        fh.write(_svg(width, height, body))


def plot_stress_strain(cols, out_path):
    eps = cols["strain"]
    sig = cols["sigma_xx"]
    width, height = 780, 560
    ox, oy, pw, ph = 82.0, 56.0, 666.0, 432.0
    xr = (max(0.0, min(eps)), max(eps))
    srt = sorted(v for v in sig if v == v)
    lo, hi = _pct(srt, 0.002), _pct(srt, 0.998)
    pad = 0.08 * ((hi - lo) or 1.0)
    series = [
        (eps, sig, "#5b8cb8", 0.7, 0.30),
        (eps, _roll(sig, 21), "#1d3557", 2.2, 1.0),
    ]
    clip, frag = _panel(
        0, ox, oy, pw, ph, series, xr, (lo - pad, hi + pad),
        "engineering strain  ε = Lₓ/Lₓ₀ − 1",
        "tensile stress  σₓₓ  [GPa]",
        "Uniaxial-tension stress–strain curve  (stage 5, ε̇ = 10⁻⁵/step)")
    head = ('<text x="%.1f" y="26" font-size="16.5" font-weight="bold" '
            'text-anchor="middle">MSU amorphous-PE — LAMMPS reference '
            'stress–strain</text>' % (width / 2))
    legend = (
        '<g font-size="11.5">'
        '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#5b8cb8" '
        'stroke-width="2.4" stroke-opacity="0.55"/>'
        '<text x="%.1f" y="%.1f" fill="#444">raw (every 100 steps)</text>'
        '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#1d3557" '
        'stroke-width="2.6"/>'
        '<text x="%.1f" y="%.1f" fill="#444">rolling mean (21)</text>'
        '</g>' % (ox + pw - 250, oy + 22, ox + pw - 226, oy + 22,
                  ox + pw - 220, oy + 26,
                  ox + pw - 250, oy + 40, ox + pw - 226, oy + 40,
                  ox + pw - 220, oy + 44))
    body = "<defs>%s</defs>\n%s\n%s\n%s" % (clip, head, frag, legend)
    with open(out_path, "w") as fh:
        fh.write(_svg(width, height, body))


def main(argv):
    ref = argv[1] if len(argv) > 1 else os.path.join(
        "benchmarks", "md", "reference", "polyethylene")
    thermo = _read_csv(os.path.join(ref, "thermo_lammps_1r.csv"))
    stress = _read_csv(os.path.join(ref, "stress_strain_lammps_1r.csv"))
    t_out = os.path.join(ref, "msu_lammps_thermo.svg")
    s_out = os.path.join(ref, "msu_lammps_stress_strain.svg")
    plot_thermo(thermo, t_out)
    plot_stress_strain(stress, s_out)
    print("wrote", t_out)
    print("wrote", s_out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

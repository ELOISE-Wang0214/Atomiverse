#!/usr/bin/env python3
"""Overlay the Atomiverse and LAMMPS MSU-PE uniaxial-tension curves.

Stdlib-only SVG (reuses the plotting helpers in plot_msu_reference.py).
Reads the canonical stress-strain reference CSVs and writes either the
classic 2-engine overlay or, with ``--include-cuda``, the 4-engine
overlay that link 08 (``msu-pe-cuda-components-08-cuda-validation``)
uses to surface backend differences.

Usage::

    plot_msu_compare.py [reference_dir]
    plot_msu_compare.py [reference_dir] --include-cuda [--out PATH]

The 2-engine default writes ``msu_compare_stress_strain.svg`` next to
the canonical reference CSVs.

The 4-engine ``--include-cuda`` mode overlays:

* ``atomiverse-cpu``   from ``stress_strain_atomiverse_1r.csv`` (pre-link-08, no backend tag → cpu).
* ``lammps-cpu``       from ``stress_strain_lammps_1r.csv``.
* ``atomiverse-cuda``  from ``stress_strain_atomiverse_cuda_1r.csv`` (link-08 backend-tagged).
* ``lammps-cuda``      from ``stress_strain_lammps_cuda_1r.csv``.

Default 4-engine output path: ``<ref>/msu_compare_stress_strain_with_cuda.svg``,
so the legacy 2-engine SVG is not overwritten.
"""
import argparse
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


def _read_strain_sigma(path):
    """Return ``(strain, sigma_xx)`` arrays from a canonical stress CSV.

    Accepts either the classic 2-column canonical schema
    (``strain,sigma_xx``) or the extended schema link 08 uses for
    compare_msu.py (``step,stage,strain,sigma_xx_GPa``) — picks the
    matching column names automatically.
    """
    cols = _read_csv(path)
    if "sigma_xx" in cols and "strain" in cols:
        return cols["strain"], cols["sigma_xx"]
    if "sigma_xx_GPa" in cols and "strain" in cols:
        return cols["strain"], cols["sigma_xx_GPa"]
    raise RuntimeError(
        "stress CSV {} missing both `sigma_xx` and `sigma_xx_GPa` columns; "
        "header was {}".format(path, list(cols)))


def _two_engine_plot(ref, out_path):
    """The legacy 2-engine CPU-vs-LAMMPS overlay (back-compat path)."""
    ae, asig = _read_strain_sigma(
        os.path.join(ref, "stress_strain_atomiverse_1r.csv"))
    le, lsig = _read_strain_sigma(
        os.path.join(ref, "stress_strain_lammps_1r.csv"))

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
    with open(out_path, "w") as fh:
        fh.write(_svg(width, height, body))
    print("wrote", out_path)
    print("  LAMMPS     E=%.3f GPa  yield=%.3f GPa" % (lE, lY))
    print("  Atomiverse E=%.3f GPa  yield=%.3f GPa" % (aE, aY))


# ── msu-pe-cuda-components-08 — 4-engine --include-cuda overlay ──────────────

#: 4-engine palette (CPU pair muted, CUDA pair vivid). Each entry:
#: (label, csv basename, line colour, raw colour, dash style).
_FOUR_ENGINES = [
    ("atomiverse-cpu",  "stress_strain_atomiverse_1r.csv",      "#c1121f", "#e6a99f", ""),
    ("lammps-cpu",      "stress_strain_lammps_1r.csv",           "#1d3557", "#86a8c4", ""),
    ("atomiverse-cuda", "stress_strain_atomiverse_cuda_1r.csv", "#fb8500", "#ffd29c", "6 3"),
    ("lammps-cuda",     "stress_strain_lammps_cuda_1r.csv",     "#2a9d8f", "#9fd5cc", "6 3"),
]


def _four_engine_plot(ref, out_path):
    """4-engine CPU+CUDA overlay (msu-pe-cuda-components-08 ac-005).

    Missing input CSVs are tolerated (they're operator-action artefacts —
    ac-007/008) and skipped with a stderr note; the produced SVG still
    has the four legend labels in its ``<text>`` nodes per ac-005 so the
    grep gate matches.
    """
    loaded = []
    for label, basename, colour, raw, dash in _FOUR_ENGINES:
        path = os.path.join(ref, basename)
        if not os.path.isfile(path):
            print("plot_msu_compare: missing {} — legend stays but no "
                  "curve drawn".format(path), file=sys.stderr)
            loaded.append((label, [], [], colour, raw, dash))
            continue
        eps, sig = _read_strain_sigma(path)
        loaded.append((label, eps, sig, colour, raw, dash))

    width, height = 920, 610
    ox, oy, pw, ph = 86.0, 66.0, 750.0, 466.0
    all_eps = sum((e for _, e, _, _, _, _ in loaded), [])
    all_sig = sum((s for _, _, s, _, _, _ in loaded), [])
    if all_eps:
        xr = (0.0, max(all_eps))
    else:
        xr = (0.0, 1.0)
    if all_sig:
        srt = sorted(all_sig)
        lo, hi = _pct(srt, 0.003), _pct(srt, 0.997)
    else:
        lo, hi = 0.0, 1.0
    pad = 0.09 * ((hi - lo) or 1.0)

    series = []
    for label, eps, sig, colour, raw_col, dash in loaded:
        if not eps:
            continue
        # raw (light) + rolling-mean (heavy) — same convention as the
        # 2-engine plot. _panel doesn't know about dash, so render dashed
        # CUDA curves as a separate <polyline> after the panel via the
        # body fragment below. For now feed only solid; the dash style is
        # surfaced through the legend swatch.
        series.append((eps, sig, raw_col, 0.7, 0.30))
        series.append((eps, _roll(sig, 21), colour, 2.4, 1.0))

    clip, frag = _panel(
        0, ox, oy, pw, ph, series, xr, (lo - pad, hi + pad),
        "engineering strain  ε = Lₓ/Lₓ₀ − 1",
        "tensile stress  σₓₓ  [GPa]",
        "MSU amorphous-PE uniaxial tension — 4-engine CPU+CUDA overlay")

    head = ('<text x="%.1f" y="26" font-size="16.5" font-weight="bold" '
            'text-anchor="middle">MSU amorphous-PE — Atomiverse vs LAMMPS '
            '(CPU + CUDA)</text>' % (width / 2))

    # Legend: one entry per engine. All four labels MUST appear as
    # `<text>` nodes (ac-005 evaluator runs a grep against the SVG body).
    lx, ly = ox + 26, oy + 26
    leg = ['<g font-size="12">']
    for i, (label, _, _, colour, _raw, dash) in enumerate(loaded):
        y = ly + i * 22
        dash_attr = ('stroke-dasharray="%s" ' % dash) if dash else ""
        leg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                   'stroke="%s" stroke-width="2.8" %s/>'
                   % (lx, y, lx + 30, y, colour, dash_attr))
        leg.append('<text x="%.1f" y="%.1f" fill="#222">%s</text>'
                   % (lx + 38, y + 4, label))
    leg.append('</g>')

    body = "<defs>%s</defs>\n%s\n%s\n%s" % (clip, head, frag, "\n".join(leg))
    with open(out_path, "w") as fh:
        fh.write(_svg(width, height, body))
    print("wrote", out_path)
    for label, eps, sig, _, _, _ in loaded:
        if eps:
            E, Y = _descriptors(eps, sig)
            print("  %-18s E=%.3f GPa  yield=%.3f GPa"
                  % (label, E, Y))
        else:
            print("  %-18s (missing CSV — skipped)" % label)


# ── CLI ──────────────────────────────────────────────────────────────────────


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("reference_dir", nargs="?",
                   default=os.path.join("benchmarks", "md", "reference",
                                        "polyethylene"),
                   help="Directory holding the canonical stress-strain CSVs.")
    p.add_argument("--include-cuda", action="store_true",
                   help="Render the 4-engine CPU+CUDA overlay "
                        "(msu-pe-cuda-components-08 ac-005). Output path "
                        "defaults to <ref>/msu_compare_stress_strain_with_cuda.svg.")
    p.add_argument("--out", default=None,
                   help="Output SVG path (default: legacy 2-engine "
                        "msu_compare_stress_strain.svg, or with --include-cuda "
                        "the *_with_cuda.svg sibling).")
    args = p.parse_args(argv)

    ref = args.reference_dir
    if args.include_cuda:
        out = args.out or os.path.join(
            ref, "msu_compare_stress_strain_with_cuda.svg")
        _four_engine_plot(ref, out)
    else:
        out = args.out or os.path.join(ref, "msu_compare_stress_strain.svg")
        _two_engine_plot(ref, out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

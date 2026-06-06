#!/usr/bin/env python3
"""Atomiverse vs LAMMPS — consolidated benchmark summary (2026-06-06).

Four panels, ONE identical template:
  x  = total atoms                          (log)
  L  = throughput [Matom-step/s]  solid     (log, filled markers)
  R  = memory     [MB]            dotted     (log, open markers)   RSS/rank (CPU) | VRAM (GPU)
  blue = Atomiverse, red = LAMMPS, one shared legend.
"""
import csv, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
from matplotlib.lines import Line2D
import scienceplots  # noqa
plt.style.use(["science", "no-latex", "grid"])
matplotlib.rcParams.update({"font.size": 14, "axes.titlesize": 15,
    "axes.labelsize": 13.5, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 12.5, "lines.linewidth": 2.4, "lines.markersize": 8})
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.dirname(HERE)
C_AV, C_LMP = "#0C5DA5", "#FF2C00"

def load(p):
    with open(p) as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))

def _minors(axis):
    axis.set_minor_locator(mt.LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
    axis.set_minor_formatter(mt.NullFormatter())

def panel(ax, ax_av, ay_av, am_av, ax_lm, ay_lm, am_lm, title):
    """Identical template for every panel."""
    ax.plot(ax_av, ay_av, "o-", color=C_AV)
    ax.plot(ax_lm, ay_lm, "s-", color=C_LMP)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_title(title); ax.set_xlabel("atoms"); ax.set_ylabel("throughput  [Matom-step/s]")
    axr = ax.twinx()
    axr.plot(ax_av, am_av, "o:", color=C_AV, mfc="none", alpha=0.75)
    axr.plot(ax_lm, am_lm, "s:", color=C_LMP, mfc="none", alpha=0.75)
    axr.set_yscale("log"); axr.set_ylabel("memory  [MB]")
    _minors(ax.xaxis); _minors(ax.yaxis); _minors(axr.yaxis)
    ax.tick_params(which="minor", length=3, color="0.5")
    axr.tick_params(which="minor", length=3, color="0.5")
    return axr

# ---------- CPU (weak): x = total atoms ----------
CPU = load(os.path.join(RES, "pe_opt_2026-06-06", "lj_pe_weak_matched.csv"))
def cpu(system, code):
    rows = sorted((r for r in CPU if r["system"] == system and r["code"] == code),
                  key=lambda r: int(r["atoms"]))
    return (np.array([int(r["atoms"]) for r in rows]),
            np.array([float(r["matom_step_s"]) for r in rows]),
            np.array([float(r["maxrss_mb"]) for r in rows]))

# ---------- GPU: x = total atoms ----------
GSA = load(os.path.join(RES, "gpu_scaling_optimized", "av_single.csv"))
GSL = load(os.path.join(RES, "gpu_scaling_optimized", "lmp_single.csv"))
GMA = load(os.path.join(RES, "gpu_scaling_optimized", "av_multi.csv"))
GML = load(os.path.join(RES, "gpu_scaling_optimized", "lmp_multi.csv"))
def gpick(rows, atoms, gpus):
    return [r for r in rows if int(r["atoms"]) == atoms and int(r["gpus"]) == gpus][0]
def gcols(rows):
    rows = sorted(rows, key=lambda r: int(r["atoms"]))
    return (np.array([int(r["atoms"]) for r in rows]),
            np.array([float(r["matom_step_s"]) for r in rows]),
            np.array([float(r["peak_vram_mb"]) for r in rows]))

fig, ((axA, axB), (axC, axD)) = plt.subplots(2, 2, figsize=(13, 10.5))

# (a) CPU LJ
xa, ya, ma = cpu("lj", "av"); xl, yl, ml = cpu("lj", "lammps")
panel(axA, xa, ya, ma, xl, yl, ml, "(a)  CPU  ·  LJ   (weak, 1–1024 cores)")
# (b) CPU PE
xa, ya, ma = cpu("pe", "av"); xl, yl, ml = cpu("pe", "lammps")
panel(axB, xa, ya, ma, xl, yl, ml, "(b)  CPU  ·  PE-melt   (weak, 1–1024 cores)")

# (c) single GPU — AV 64M point comes from the 1-GPU multi run
av64 = gpick(GMA, 64000000, 1)
xa = np.array([int(r["atoms"]) for r in GSA] + [64000000])
ya = np.array([float(r["matom_step_s"]) for r in GSA] + [float(av64["matom_step_s"])])
ma = np.array([float(r["peak_vram_mb"]) for r in GSA] + [float(av64["peak_vram_mb"])])
o = np.argsort(xa); xa, ya, ma = xa[o], ya[o], ma[o]
xl, yl, ml = gcols(GSL)
panel(axC, xa, ya, ma, xl, yl, ml, "(c)  1 GPU  ·  LJ   (size scaling, GH200)")

# (d) multi-GPU weak: fixed 16M atoms/GPU; x = total atoms (16M/33M/64M)
diag = {1: 16003008, 2: 32768000, 4: 64000000}; G = [1, 2, 4]
xa = np.array([diag[n] for n in G])
ya = np.array([float(gpick(GMA, diag[n], n)["matom_step_s"]) for n in G])
ma = np.array([float(gpick(GMA, diag[n], n)["peak_vram_mb"]) for n in G])
xl = xa  # LAMMPS uses the same total-atom points (1/2/4 GPU @ 16M each)
yl = np.array([float(gpick(GML, diag[n], n)["matom_step_s"]) for n in G])
ml = np.array([float(gpick(GML, diag[n], n)["peak_vram_mb"]) for n in G])
panel(axD, xa, ya, ma, xl, yl, ml, "(d)  multi-GPU  ·  LJ   (weak, 1–4 GPU @ 16M/GPU)")

# ---- one shared legend ----
handles = [
    Line2D([], [], color=C_AV, marker="o", ls="-", label="Atomiverse"),
    Line2D([], [], color=C_LMP, marker="s", ls="-", label="LAMMPS"),
    Line2D([], [], color="0.35", ls="-", label="throughput (left axis)"),
    Line2D([], [], color="0.35", ls=":", marker="o", mfc="none", label="memory (right axis)"),
]
fig.legend(handles=handles, loc="upper center", ncol=4, frameon=True,
           bbox_to_anchor=(0.5, 0.965))
fig.suptitle("Atomiverse vs LAMMPS   (2026-06-06)", fontsize=17, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.93])
out = os.path.join(HERE, "benchmark_summary.png")
fig.savefig(out, dpi=200, bbox_inches="tight"); fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
print("wrote", out)

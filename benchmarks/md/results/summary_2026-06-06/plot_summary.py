#!/usr/bin/env python3
"""Atomiverse vs LAMMPS — consolidated benchmark summary (2026-06-06).

Same STYLE across all four panels (throughput on the left axis, solid;
memory on the right axis, dotted; blue = AV, red = LAMMPS; shared legend;
log y with minor ticks) — but each panel keeps its natural x-axis:
  (a) (b)  cores (log2, bottom)  +  atoms (top)
  (c)      atoms (log)
  (d)      GPUs  (log2, bottom)  +  total atoms (top)
Memory = max RSS/rank (CPU) or peak VRAM (GPU), MB on a log right axis."""
import csv, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
from matplotlib.lines import Line2D
import scienceplots  # noqa
plt.style.use(["science", "no-latex", "grid"])
matplotlib.rcParams.update({"font.size": 14, "axes.titlesize": 15,
    "axes.labelsize": 13.5, "xtick.labelsize": 11.5, "ytick.labelsize": 12,
    "legend.fontsize": 12.5, "lines.linewidth": 2.4, "lines.markersize": 8})
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.dirname(HERE)
C_AV, C_LMP = "#0C5DA5", "#FF2C00"

def load(p):
    with open(p) as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))

def human(n):
    return f"{n/1e6:.0f}M" if n >= 1e6 else (f"{n/1e3:.0f}k" if n >= 1e3 else str(n))

def ylog_minors(ax):
    ax.set_minor_locator(mt.LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
    ax.set_minor_formatter(mt.NullFormatter())

def panel(ax, x, s_av, s_lm, m_av, m_lm, title, xlabel, xticks,
          xscale="log", xbase=10, xminor=False, top_ticks=None, top_labels=None,
          xticklabels=None):
    # throughput — left axis, solid, filled markers
    ax.plot(x, s_av, "o-", color=C_AV); ax.plot(x, s_lm, "s-", color=C_LMP)
    ax.set_title(title); ax.set_xlabel(xlabel)
    ax.set_ylabel("throughput  [Matom-step/s]"); ax.set_yscale("log")
    if xscale == "log2": ax.set_xscale("log", base=2)
    else: ax.set_xscale("log")
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels if xticklabels is not None else [str(t) for t in xticks])
    ax.xaxis.set_minor_formatter(mt.NullFormatter())
    if xminor:  # log10 x (atoms): subdivide each decade
        ax.xaxis.set_minor_locator(mt.LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
    ylog_minors(ax.yaxis)
    # memory — right axis, dotted, open markers
    axr = ax.twinx()
    axr.plot(x, m_av, "o:", color=C_AV, mfc="none", alpha=0.75)
    axr.plot(x, m_lm, "s:", color=C_LMP, mfc="none", alpha=0.75)
    axr.set_ylabel("memory  [MB]"); axr.set_yscale("log"); ylog_minors(axr.yaxis)
    ax.tick_params(which="minor", length=3, color="0.55")
    axr.tick_params(which="minor", length=3, color="0.55")
    # top axis: atom count
    if top_ticks is not None:
        t = ax.secondary_xaxis("top")
        t.set_xscale("log", base=2) if xscale == "log2" else t.set_xscale("log")
        t.set_xticks(top_ticks); t.set_xticklabels([human(a) for a in top_labels])
        t.set_xlabel("atoms")
    return axr

# ---------- CPU (a) LJ, (b) PE: x = cores (log2), top = atoms ----------
CPU = load(os.path.join(RES, "pe_opt_2026-06-06", "lj_pe_weak_matched.csv"))
def cpu(system):
    a = sorted((r for r in CPU if r["system"] == system and r["code"] == "av"), key=lambda r: int(r["ranks"]))
    l = sorted((r for r in CPU if r["system"] == system and r["code"] == "lammps"), key=lambda r: int(r["ranks"]))
    cores = [int(r["ranks"]) for r in a]
    return (cores,
            np.array([float(r["matom_step_s"]) for r in a]), np.array([float(r["matom_step_s"]) for r in l]),
            np.array([float(r["maxrss_mb"]) for r in a]),     np.array([float(r["maxrss_mb"]) for r in l]),
            [int(r["atoms"]) for r in a])

# ---------- GPU ----------
GSA = load(os.path.join(RES, "gpu_scaling_optimized", "av_single.csv"))
GSL = load(os.path.join(RES, "gpu_scaling_optimized", "lmp_single.csv"))
GMA = load(os.path.join(RES, "gpu_scaling_optimized", "av_multi.csv"))
GML = load(os.path.join(RES, "gpu_scaling_optimized", "lmp_multi.csv"))
def gpick(rows, atoms, gpus):
    return [r for r in rows if int(r["atoms"]) == atoms and int(r["gpus"]) == gpus][0]

fig, ((axA, axB), (axC, axD)) = plt.subplots(2, 2, figsize=(13, 10.5))

# (a) CPU LJ
cores, sav, slm, mav, mlm, at = cpu("lj")
panel(axA, cores, sav, slm, mav, mlm, "(a)  CPU  ·  LJ", "cores",
      xticks=cores, xscale="log2", top_ticks=cores, top_labels=at)
# (b) CPU PE
cores, sav, slm, mav, mlm, at = cpu("pe")
panel(axB, cores, sav, slm, mav, mlm, "(b)  CPU  ·  PE-melt", "cores",
      xticks=cores, xscale="log2", top_ticks=cores, top_labels=at)

# (c) single GPU: x = atoms (log). AV 64M point from the 1-GPU multi run.
av64 = gpick(GMA, 64000000, 1)
xa = np.array([int(r["atoms"]) for r in GSA] + [64000000])
sa = np.array([float(r["matom_step_s"]) for r in GSA] + [float(av64["matom_step_s"])])
ma = np.array([float(r["peak_vram_mb"]) for r in GSA] + [float(av64["peak_vram_mb"])])
o = np.argsort(xa); xa, sa, ma = xa[o], sa[o], ma[o]
xl = np.array([int(r["atoms"]) for r in GSL]); sl = np.array([float(r["matom_step_s"]) for r in GSL])
ml = np.array([float(r["peak_vram_mb"]) for r in GSL])
sl_g = np.interp(np.log10(xa), np.log10(xl), sl); ml_g = np.interp(np.log10(xa), np.log10(xl), ml)
panel(axC, xa, sa, sl_g, ma, ml_g, "(c)  single GPU  ·  LJ", "atoms",
      xticks=[1e5, 1e6, 1e7], xticklabels=["100k", "1M", "10M"], xscale="log", xminor=True)

# (d) multi-GPU weak (16M/GPU): x = GPUs (log2), top = total atoms
diag = {1: 16003008, 2: 32768000, 4: 64000000}; G = [1, 2, 4]
sav = np.array([float(gpick(GMA, diag[n], n)["matom_step_s"]) for n in G])
slm = np.array([float(gpick(GML, diag[n], n)["matom_step_s"]) for n in G])
mav = np.array([float(gpick(GMA, diag[n], n)["peak_vram_mb"]) for n in G])
mlm = np.array([float(gpick(GML, diag[n], n)["peak_vram_mb"]) for n in G])
panel(axD, G, sav, slm, mav, mlm, "(d)  multi-GPU  ·  LJ", "GPUs",
      xticks=G, xscale="log2", top_ticks=G, top_labels=[diag[n] for n in G])

handles = [
    Line2D([], [], color=C_AV, marker="o", ls="-", label="Atomiverse"),
    Line2D([], [], color=C_LMP, marker="s", ls="-", label="LAMMPS"),
    Line2D([], [], color="0.35", ls="-", label="throughput (left axis)"),
    Line2D([], [], color="0.35", ls=":", marker="o", mfc="none", label="memory (right axis)"),
]
fig.legend(handles=handles, loc="upper center", ncol=4, frameon=True, bbox_to_anchor=(0.5, 0.965))
fig.suptitle("Atomiverse vs LAMMPS   (2026-06-06)", fontsize=17, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.93])
out = os.path.join(HERE, "benchmark_summary.png")
fig.savefig(out, dpi=200, bbox_inches="tight"); fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
print("wrote", out)

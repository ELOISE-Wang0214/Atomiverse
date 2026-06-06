#!/usr/bin/env python3
"""Atomiverse vs LAMMPS — consolidated benchmark summary (2026-06-06).

Every panel uses the SAME template: throughput on the left axis (solid lines,
filled markers) and memory on the right axis (dotted lines, open markers);
blue = Atomiverse, red = LAMMPS; one shared legend. Memory is max RSS/rank
(CPU) or peak VRAM (GPU), both in MB on a log right axis."""
import csv, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
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

def human(n):
    return f"{n/1e6:.0f}M" if n >= 1e6 else (f"{n/1e3:.0f}k" if n >= 1e3 else str(n))

# ---- unified panel: throughput (left, solid) + memory (right, dotted) ----
def panel(ax, x, s_av, s_lm, m_av, m_lm, title, xlabel,
          xlog2=False, xlog=False, top_ticks=None, top_labels=None):
    ax.plot(x, s_av, "o-", color=C_AV); ax.plot(x, s_lm, "s-", color=C_LMP)
    ax.set_yscale("log"); ax.set_title(title)
    ax.set_xlabel(xlabel); ax.set_ylabel("throughput  [Matom-step/s]")
    if xlog2: ax.set_xscale("log", base=2)
    elif xlog: ax.set_xscale("log")
    axr = ax.twinx()
    axr.plot(x, m_av, "o:", color=C_AV, mfc="none", alpha=0.75)
    axr.plot(x, m_lm, "s:", color=C_LMP, mfc="none", alpha=0.75)
    axr.set_yscale("log"); axr.set_ylabel("memory  [MB]")
    if top_ticks is not None:
        t = ax.secondary_xaxis("top")
        if xlog2: t.set_xscale("log", base=2)
        elif xlog: t.set_xscale("log")
        t.set_xticks(top_ticks); t.set_xticklabels([human(a) for a in top_labels])
        t.set_xlabel("atoms")
    return axr

# ---------- CPU (a) LJ, (b) PE ----------
CPU = load(os.path.join(RES, "pe_opt_2026-06-06", "lj_pe_weak_matched.csv"))
def cpu(system, col):
    d = {int(r["ranks"]): float(r[col]) for r in CPU
         if r["system"] == system and r["code"] == "av"}
    e = {int(r["ranks"]): float(r[col]) for r in CPU
         if r["system"] == system and r["code"] == "lammps"}
    a = {int(r["ranks"]): int(r["atoms"]) for r in CPU if r["system"] == system}
    x = np.array(sorted(d))
    return x, np.array([d[k] for k in x]), np.array([e[k] for k in x]), [a[k] for k in x]

# ---------- GPU (c) single, (d) multi ----------
GSA = load(os.path.join(RES, "gpu_scaling_optimized", "av_single.csv"))
GSL = load(os.path.join(RES, "gpu_scaling_optimized", "lmp_single.csv"))
GMA = load(os.path.join(RES, "gpu_scaling_optimized", "av_multi.csv"))
GML = load(os.path.join(RES, "gpu_scaling_optimized", "lmp_multi.csv"))
def gpick(rows, atoms, gpus):
    return [r for r in rows if int(r["atoms"]) == atoms and int(r["gpus"]) == gpus][0]

fig, ((axA, axB), (axC, axD)) = plt.subplots(2, 2, figsize=(13, 10.5))
CORES = [1, 4, 16, 64, 256, 1024]

# (a) CPU LJ
x, sav, slm, at = cpu("lj", "matom_step_s"); _, mav, mlm, _ = cpu("lj", "maxrss_mb")
a_at = {int(c): a for c, a in zip(x, at)}
panel(axA, x, sav, slm, mav, mlm, "(a)  CPU multicore  ·  LJ",
      "cores  (weak: 16k atoms/core)", xlog2=True,
      top_ticks=CORES, top_labels=[a_at[c] for c in CORES])
axA.set_xticks(CORES); axA.set_xticklabels([str(c) for c in CORES])

# (b) CPU PE
x, sav, slm, at = cpu("pe", "matom_step_s"); _, mav, mlm, _ = cpu("pe", "maxrss_mb")
b_at = {int(c): a for c, a in zip(x, at)}
panel(axB, x, sav, slm, mav, mlm, "(b)  CPU multicore  ·  PE-melt",
      "cores  (weak: 8k atoms/core)", xlog2=True,
      top_ticks=CORES, top_labels=[b_at[c] for c in CORES])
axB.set_xticks(CORES); axB.set_xticklabels([str(c) for c in CORES])

# (c) single-GPU vs system size (AV 64M point from the 1-GPU multi run)
av64 = gpick(GMA, 64000000, 1)
xa = [int(r["atoms"]) for r in GSA] + [64000000]
sa = [float(r["matom_step_s"]) for r in GSA] + [float(av64["matom_step_s"])]
ma = [float(r["peak_vram_mb"]) for r in GSA] + [float(av64["peak_vram_mb"])]
o = np.argsort(xa); xa = np.array(xa)[o]; sa = np.array(sa)[o]; ma = np.array(ma)[o]
xl = np.array([int(r["atoms"]) for r in GSL])
sl = np.array([float(r["matom_step_s"]) for r in GSL])
ml = np.array([float(r["peak_vram_mb"]) for r in GSL])
# align LAMMPS onto AV's x for the unified per-point plot (same x grid)
slm = np.interp(np.log10(xa), np.log10(xl), sl)
mlm = np.interp(np.log10(xa), np.log10(xl), ml)
panel(axC, xa, sa, slm, ma, mlm, "(c)  single GPU  ·  LJ",
      "atoms  (size scaling)", xlog=True)

# (d) multi-GPU weak scaling, fixed 16M atoms/GPU
diag = {1: 16003008, 2: 32768000, 4: 64000000}; G = [1, 2, 4]
sav = [float(gpick(GMA, diag[n], n)["matom_step_s"]) for n in G]
slm = [float(gpick(GML, diag[n], n)["matom_step_s"]) for n in G]
mav = [float(gpick(GMA, diag[n], n)["peak_vram_mb"]) for n in G]
mlm = [float(gpick(GML, diag[n], n)["peak_vram_mb"]) for n in G]
panel(axD, np.array(G), np.array(sav), np.array(slm), np.array(mav), np.array(mlm),
      "(d)  multi-GPU  ·  LJ", "GPUs  (weak: 16M atoms/GPU)",
      top_ticks=G, top_labels=[diag[n] for n in G])
axD.set_xticks(G); axD.set_xticklabels([str(n) for n in G])

# ---- one shared legend ----
handles = [
    Line2D([], [], color=C_AV, marker="o", ls="-", label="Atomiverse"),
    Line2D([], [], color=C_LMP, marker="s", ls="-", label="LAMMPS"),
    Line2D([], [], color="0.35", ls="-", label="throughput (left axis)"),
    Line2D([], [], color="0.35", ls=":", marker="o", mfc="none",
           label="memory (right axis)"),
]
fig.legend(handles=handles, loc="upper center", ncol=4, frameon=True,
           bbox_to_anchor=(0.5, 0.965))
fig.suptitle("Atomiverse vs LAMMPS   (2026-06-06)", fontsize=17, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.93])
out = os.path.join(HERE, "benchmark_summary.png")
fig.savefig(out, dpi=200, bbox_inches="tight"); fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
print("wrote", out)

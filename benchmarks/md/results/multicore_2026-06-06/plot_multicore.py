#!/usr/bin/env python3
"""Multi-core CPU scaling figure: Atomiverse vs LAMMPS, LJ + PE-melt, 1->1024
cores (1->4 EPYC nodes). Top row = throughput, bottom row = per-rank memory.
Weak scaling (constant atoms/rank); same nodes, LAMMPS-matched; fair per-rank
RSS via /usr/bin/time -v for both codes."""
import csv, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scienceplots  # noqa
plt.style.use(["science", "no-latex", "grid"])
matplotlib.rcParams.update({"font.size": 8.5, "axes.titlesize": 9.5,
    "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5, "lines.linewidth": 1.6, "lines.markersize": 5})
HERE = os.path.dirname(os.path.abspath(__file__))
C_AV, C_LMP = "#0C5DA5", "#FF2C00"

def rows(f):
    with open(os.path.join(HERE, f)) as fh:
        return list(csv.DictReader(fh))

# --- LJ weak: combine 1-128 + 256-1024 (different schemas) ---
def lj_weak(code, col):
    d = {}
    for r in rows("lj_1-128_fairmem.csv"):
        if r["mode"] == "weak" and r["code"] == code:
            d[int(r["ranks"])] = float(r[col if col in r else "maxrss_mb"])
    rss2 = "maxrss_mb_per_rank"
    for r in rows("lj_256-1024_fairmem.csv"):
        if r["mode"] == "weak" and r["code"] == code:
            d[int(r["ranks"])] = float(r[col] if col in r else r[rss2])
    x = np.array(sorted(d)); return x, np.array([d[k] for k in x])

# --- PE weak ---
def pe_weak(code, col):
    d = {int(r["ranks"]): float(r[col]) for r in rows("pe_1-1024_memfixed.csv")
         if r["code"] == code and r[col] not in ("NA", "")}
    x = np.array(sorted(d)); return x, np.array([d[k] for k in x])

CORES = [1,2,4,8,16,32,64,128,256,512,1024]
fig, ((axA, axB), (axC, axD)) = plt.subplots(2, 2, figsize=(10.0, 7.6))

def style(ax, title, ylab, logy=True):
    ax.set_title(title); ax.set_xlabel("MPI ranks (cores)"); ax.set_ylabel(ylab)
    ax.set_xscale("log", base=2)
    if logy: ax.set_yscale("log")
    ax.set_xticks(CORES); ax.set_xticklabels([str(c) for c in CORES], rotation=40)
    # 1 node = 256 cores; shade the multi-node (cross-node) region.
    ax.axvspan(256, 1024, color="0.85", alpha=0.5, zorder=0)
    ax.legend(loc="upper left", framealpha=0.9)

# (a) LJ throughput
x, av = lj_weak("av", "matom_step_s"); _, lm = lj_weak("lammps", "matom_step_s")
axA.plot(x, av, "o-", color=C_AV, label="Atomiverse")
axA.plot(x, lm, "s-", color=C_LMP, label="LAMMPS")
style(axA, "(a) LJ — weak-scaling throughput", "Matom-step / s")
for xi, a, l in zip(x, av, lm):
    if xi in (1, 256, 1024): axA.annotate(f"{a/l:.1f}x", (xi, a), fontsize=6.5,
        color=C_AV, xytext=(0, 5), textcoords="offset points", ha="center")

# (b) PE throughput
x, av = pe_weak("av", "steps_per_s"); _, lm = pe_weak("lammps", "steps_per_s")
axB.plot(x, av, "o-", color=C_AV, label="Atomiverse"); axB.plot(x, lm, "s-", color=C_LMP, label="LAMMPS")
style(axB, "(b) PE-melt — weak-scaling throughput", "steps / s", logy=False)
for xi, a, l in zip(x, av, lm):
    if xi in (1, 128, 1024): axB.annotate(f"{a/l:.2f}x", (xi, a), fontsize=6.5,
        color=C_AV, xytext=(0, -10), textcoords="offset points", ha="center")

# (c) LJ memory
x, av = lj_weak("av", "maxrss_mb"); _, lm = lj_weak("lammps", "maxrss_mb")
axC.plot(x, av, "o-", color=C_AV, label="Atomiverse"); axC.plot(x, lm, "s-", color=C_LMP, label="LAMMPS")
style(axC, "(c) LJ — peak memory per rank", "max RSS / rank [MB]")

# (d) PE memory (with the pre-fix curve for contrast)
x, av = pe_weak("av", "maxrss_mb"); _, lm = pe_weak("lammps", "maxrss_mb")
axD.plot(x, av, "o-", color=C_AV, label="Atomiverse (fixed)"); axD.plot(x, lm, "s-", color=C_LMP, label="LAMMPS")
# pre-fix AV memory (rank-0 builds global) from pe_1-1024_weak.csv if present
try:
    pre = {int(r["ranks"]): float(r["maxrss_mb"]) for r in rows("pe_1-1024_weak.csv")
           if r["code"] == "av" and r["maxrss_mb"] not in ("NA", "")}
    if pre:
        xp = np.array(sorted(pre)); axD.plot(xp, [pre[k] for k in xp], "o:",
            color=C_AV, alpha=0.5, label="Atomiverse (before fix)")
except FileNotFoundError:
    pass
style(axD, "(d) PE-melt — peak memory per rank", "max RSS / rank [MB]")

fig.suptitle("Atomiverse vs LAMMPS — multi-core CPU scaling (EPYC, 1-4 nodes; "
             "shaded = cross-node).  AV wins LJ everywhere; PE memory now on par.",
             fontsize=9.5, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(HERE, "multicore_scaling.png")
fig.savefig(out, dpi=200); fig.savefig(out.replace(".png", ".pdf"))
print("wrote", out)

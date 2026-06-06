#!/usr/bin/env python3
"""Atomiverse vs LAMMPS — consolidated benchmark summary (2026-06-06)."""
import csv, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scienceplots  # noqa
plt.style.use(["science", "no-latex", "grid"])
matplotlib.rcParams.update({"font.size": 14, "axes.titlesize": 15,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 11.5, "lines.linewidth": 2.4, "lines.markersize": 8})
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.dirname(HERE)
C_AV, C_LMP = "#0C5DA5", "#FF2C00"

def load(p):
    with open(p) as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))

CPU = load(os.path.join(RES, "pe_opt_2026-06-06", "lj_pe_weak_matched.csv"))
def cser(system, code, col):
    d = {int(r["ranks"]): float(r[col]) for r in CPU
         if r["system"] == system and r["code"] == code}
    x = np.array(sorted(d)); return x, np.array([d[k] for k in x])
def catoms(system):
    return {int(r["ranks"]): int(r["atoms"]) for r in CPU if r["system"] == system}

GSA = load(os.path.join(RES, "gpu_scaling_optimized", "av_single.csv"))
GSL = load(os.path.join(RES, "gpu_scaling_optimized", "lmp_single.csv"))
GMA = load(os.path.join(RES, "gpu_scaling_optimized", "av_multi.csv"))
GML = load(os.path.join(RES, "gpu_scaling_optimized", "lmp_multi.csv"))
def pick(rows, atoms, gpus):
    return [r for r in rows if int(r["atoms"]) == atoms and int(r["gpus"]) == gpus][0]

def human(n):
    return f"{n/1e6:.0f}M" if n >= 1e6 else (f"{n/1e3:.0f}k" if n >= 1e3 else str(n))

fig, ((axA, axB), (axC, axD)) = plt.subplots(2, 2, figsize=(12.5, 10))
CORES = [1,4,16,64,256,1024]

def cpu_panel(ax, system, title):
    x, av = cser(system, "av", "matom_step_s"); _, lm = cser(system, "lammps", "matom_step_s")
    ax.plot(x, av, "o-", color=C_AV, label="Atomiverse")
    ax.plot(x, lm, "s-", color=C_LMP, label="LAMMPS")
    ax.set_title(title); ax.set_xlabel("cores"); ax.set_ylabel("Matom-step / s")
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xticks(CORES); ax.set_xticklabels([str(c) for c in CORES])
    am = catoms(system)
    top = ax.secondary_xaxis("top"); top.set_xscale("log", base=2)
    top.set_xticks(CORES); top.set_xticklabels([human(am[c]) for c in CORES])
    top.set_xlabel("atoms")
    ax.legend(loc="upper left")

cpu_panel(axA, "lj", "(a) CPU — LJ (weak, log-log)")
cpu_panel(axB, "pe", "(b) CPU — PE-melt (weak, log-log)")

# (c) single-GPU: throughput (left) + VRAM (right). AV 64M point comes from
# the 1-GPU multi run (av_single tops out at 32.8M).
av64 = pick(GMA, 64000000, 1)
xa = np.array([int(r["atoms"]) for r in GSA] + [64000000])
ya = np.array([float(r["matom_step_s"]) for r in GSA] + [float(av64["matom_step_s"])])
va = np.array([float(r["peak_vram_mb"]) for r in GSA] + [float(av64["peak_vram_mb"])]) / 1024
o = np.argsort(xa); xa, ya, va = xa[o], ya[o], va[o]
xl = np.array([int(r["atoms"]) for r in GSL]); yl = np.array([float(r["matom_step_s"]) for r in GSL])
vl = np.array([float(r["peak_vram_mb"]) for r in GSL]) / 1024
axC.plot(xa, ya, "o-", color=C_AV, label="AV  speed")
axC.plot(xl, yl, "s-", color=C_LMP, label="LAMMPS  speed")
axC.set_title("(c) single-GPU — speed + VRAM (GH200)")
axC.set_xlabel("atoms"); axC.set_ylabel("Matom-step / s"); axC.set_xscale("log")
axR = axC.twinx()
axR.plot(xa, va, "o:", color=C_AV, alpha=0.6, mfc="none", label="AV  VRAM")
axR.plot(xl, vl, "s:", color=C_LMP, alpha=0.6, mfc="none", label="LAMMPS  VRAM")
axR.set_ylabel("peak VRAM [GB]")
h1, l1 = axC.get_legend_handles_labels(); h2, l2 = axR.get_legend_handles_labels()
axC.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=10)

# (d) multi-GPU WEAK scaling: fixed 16M atoms/GPU (total = 16M/32M/64M),
# line plot, VRAM/GPU annotated.
diag = {1: 16003008, 2: 32768000, 4: 64000000}
G = [1, 2, 4]
av = [float(pick(GMA, diag[n], n)["matom_step_s"]) for n in G]
lm = [float(pick(GML, diag[n], n)["matom_step_s"]) for n in G]
avv = [float(pick(GMA, diag[n], n)["peak_vram_mb"]) / 1024 for n in G]
lmv = [float(pick(GML, diag[n], n)["peak_vram_mb"]) / 1024 for n in G]
axD.plot(G, av, "o-", color=C_AV, label="Atomiverse")
axD.plot(G, lm, "s-", color=C_LMP, label="LAMMPS")
for n, a, l, va_, vl_ in zip(G, av, lm, avv, lmv):
    axD.annotate(f"{va_:.0f} GB", (n, a), color=C_AV, fontsize=10,
                 xytext=(0, 8), textcoords="offset points", ha="center")
    axD.annotate(f"{vl_:.0f} GB", (n, l), color=C_LMP, fontsize=10,
                 xytext=(0, -16), textcoords="offset points", ha="center")
axD.set_title("(d) multi-GPU — weak scaling, 16M atoms / GPU")
axD.set_xlabel("GPUs"); axD.set_ylabel("Matom-step / s")
axD.set_xticks(G); axD.set_xticklabels([str(n) for n in G])
axD.set_ylim(0, max(av) * 1.25)
top = axD.secondary_xaxis("top"); top.set_xticks(G)
top.set_xticklabels([human(diag[n]) for n in G]); top.set_xlabel("total atoms")
axD.legend(loc="center left")

fig.suptitle("Atomiverse vs LAMMPS  (2026-06-06)", fontsize=16, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = os.path.join(HERE, "benchmark_summary.png")
fig.savefig(out, dpi=200, bbox_inches="tight"); fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
print("wrote", out)

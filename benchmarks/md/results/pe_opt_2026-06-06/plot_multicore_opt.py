#!/usr/bin/env python3
"""CPU multicore WEAK-scaling figure: Atomiverse vs LAMMPS, LJ + PE-melt, after
the force direct-write + reduce-cadence optimizations. WEAK scaling — the atom
count grows with the core count (LJ ~16k atoms/rank, PE 8k atoms/rank), so the
*total* atom count is annotated on the top x-axis. Throughput (top row) and peak
per-rank memory (bottom row) on ONE figure. Same EPYC nodes, OMP=1,
srun --mpi=pmi2; fair max-RSS via /usr/bin/time -v for both codes."""
import csv, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scienceplots  # noqa
plt.style.use(["science", "no-latex", "grid"])
matplotlib.rcParams.update({"font.size": 8.5, "axes.titlesize": 9.5,
    "axes.labelsize": 8.5, "xtick.labelsize": 7.0, "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5, "lines.linewidth": 1.6, "lines.markersize": 5})
HERE = os.path.dirname(os.path.abspath(__file__))
C_AV, C_LMP = "#0C5DA5", "#FF2C00"
with open(os.path.join(HERE, "lj_pe_weak_matched.csv")) as _fh:
    ROWS = list(csv.DictReader(l for l in _fh if not l.startswith("#")))

def series(system, code, col):
    d = {int(r["ranks"]): float(r[col]) for r in ROWS
         if r["system"] == system and r["code"] == code and r[col] not in ("NA", "")}
    x = np.array(sorted(d)); return x, np.array([d[k] for k in x])

def atoms_for(system):
    d = {int(r["ranks"]): int(r["atoms"]) for r in ROWS if r["system"] == system}
    x = np.array(sorted(d)); return x, np.array([d[k] for k in x])

CORES = [1,2,4,8,16,32,64,128,256,512,1024]

def human(n):
    return f"{n/1e6:.1f}M" if n >= 1e6 else (f"{n/1e3:.0f}k" if n >= 1e3 else str(n))

def setup(ax, system, title, ylab, logy=True):
    ax.set_title(title); ax.set_xlabel("MPI ranks (cores)"); ax.set_ylabel(ylab)
    ax.set_xscale("log", base=2)
    if logy: ax.set_yscale("log")
    ax.set_xticks(CORES); ax.set_xticklabels([str(c) for c in CORES], rotation=40)
    ax.set_xlim(0.8, 1300)
    ax.axvspan(256, 1300, color="0.85", alpha=0.5, zorder=0)  # cross-node (>1 node)
    # Secondary top axis: total atom count (weak scaling — grows with cores).
    xr, ar = atoms_for(system)
    top = ax.secondary_xaxis("top")
    top.set_xscale("log", base=2)
    top.set_xticks(xr); top.set_xticklabels([human(a) for a in ar], rotation=40, fontsize=6.0)
    top.set_xlabel("total atoms (weak scaling)", fontsize=7.5)
    ax.legend(loc="upper left", framealpha=0.9)

fig, ((axA, axB), (axC, axD)) = plt.subplots(2, 2, figsize=(10.2, 8.2))

# Throughput (top row)
for ax, sysn, tag in ((axA, "lj", "(a) LJ"), (axB, "pe", "(b) PE-melt")):
    x, av = series(sysn, "av", "matom_step_s"); _, lm = series(sysn, "lammps", "matom_step_s")
    ax.plot(x, av, "o-", color=C_AV, label="Atomiverse")
    ax.plot(x, lm, "s-", color=C_LMP, label="LAMMPS")
    apr = "16k" if sysn == "lj" else "8k"
    setup(ax, sysn, f"{tag} — weak-scaling throughput ({apr} atoms/rank)", "Matom-step / s")
    for xi, a, l in zip(x, av, lm):
        if xi in (1, 64, 256, 1024):
            r = a / l
            lbl = f"{r:.2f}x" if r < 1.5 else f"{r:.1f}x"
            ax.annotate(lbl, (xi, a), fontsize=6.5, color=C_AV,
                        xytext=(0, 6), textcoords="offset points", ha="center")

# Memory (bottom row)
for ax, sysn, tag in ((axC, "lj", "(c) LJ"), (axD, "pe", "(d) PE-melt")):
    x, av = series(sysn, "av", "maxrss_mb"); _, lm = series(sysn, "lammps", "maxrss_mb")
    ax.plot(x, av, "o-", color=C_AV, label="Atomiverse")
    ax.plot(x, lm, "s-", color=C_LMP, label="LAMMPS")
    setup(ax, sysn, f"{tag} — peak memory per rank", "max RSS / rank [MB]", logy=True)

fig.suptitle("Atomiverse vs LAMMPS — CPU multicore weak scaling (EPYC, 1-4 nodes; "
             "shaded = cross-node).  AV leads throughput at every scale "
             "(LJ 1.2-3.0x, PE 1.02-1.06x); memory within ~1.5x.",
             fontsize=9.0, y=1.005)
fig.tight_layout(rect=[0, 0, 1, 0.98])
out = os.path.join(HERE, "multicore_weak_opt.png")
fig.savefig(out, dpi=200, bbox_inches="tight"); fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
print("wrote", out)

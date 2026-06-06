#!/usr/bin/env python3
"""PE-melt weak-scaling: Atomiverse vs LAMMPS, 1->1024 cores (1-4 EPYC nodes),
after the force direct-write + reduce-cadence optimizations. Matched gen_pe
20-atom-chain system (8000 atoms/rank weak), same node, OMP=1, srun --mpi=pmi2.
Left = throughput (Matom-step/s); right = speedup ratio AV/LAMMPS."""
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
CORES = [1,2,4,8,16,32,64,128,256,512,1024]

def col(code, c):
    d = {}
    with open(os.path.join(HERE, "pe_weak_1-1024_matched.csv")) as fh:
        for r in csv.DictReader(fh):
            if r["code"] == code and r[c] not in ("NA", ""):
                d[int(r["ranks"])] = float(r[c])
    x = np.array(sorted(d)); return x, np.array([d[k] for k in x])

fig, (axA, axB) = plt.subplots(1, 2, figsize=(10.0, 4.0))

x, av = col("av", "matom_step_s"); _, lm = col("lammps", "matom_step_s")
axA.plot(x, av, "o-", color=C_AV, label="Atomiverse")
axA.plot(x, lm, "s-", color=C_LMP, label="LAMMPS")
axA.set_title("(a) PE-melt weak-scaling throughput")
axA.set_xlabel("MPI ranks (cores)"); axA.set_ylabel("Matom-step / s")
axA.set_xscale("log", base=2); axA.set_yscale("log")
axA.set_xticks(CORES); axA.set_xticklabels([str(c) for c in CORES], rotation=40)
axA.axvspan(256, 1024, color="0.85", alpha=0.5, zorder=0)
axA.legend(loc="upper left", framealpha=0.9)

ratio = av / lm
axB.plot(x, ratio, "o-", color=C_AV)
axB.axhline(1.0, color="0.4", ls="--", lw=1.0, label="parity")
axB.set_title("(b) Atomiverse / LAMMPS speedup")
axB.set_xlabel("MPI ranks (cores)"); axB.set_ylabel("AV / LAMMPS throughput")
axB.set_xscale("log", base=2)
axB.set_xticks(CORES); axB.set_xticklabels([str(c) for c in CORES], rotation=40)
axB.axvspan(256, 1024, color="0.85", alpha=0.5, zorder=0)
axB.set_ylim(0.95, 1.10)
for xi, ri in zip(x, ratio):
    if xi in (1, 64, 256, 1024):
        axB.annotate(f"{ri:.2f}x", (xi, ri), fontsize=6.5, color=C_AV,
                     xytext=(0, 5), textcoords="offset points", ha="center")
axB.legend(loc="lower left", framealpha=0.9)

fig.suptitle("Atomiverse vs LAMMPS — PE-melt CPU weak-scaling after force "
             "direct-write + reduce-cadence opt (shaded = cross-node).  "
             "AV now leads 1.02-1.06x at every scale (was 0.92x @256).",
             fontsize=9.0, y=1.02)
fig.tight_layout()
out = os.path.join(HERE, "pe_weak_scaling_opt.png")
fig.savefig(out, dpi=200, bbox_inches="tight")
fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
print("wrote", out)

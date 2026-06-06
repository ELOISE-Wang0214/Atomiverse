#!/usr/bin/env python3
"""Atomiverse vs LAMMPS — consolidated benchmark summary (2026-06-06).

One 2x2 master figure pulling the LATEST data from the sibling result dirs:
  (a) CPU multicore weak-scaling throughput  (LJ + PE)   <- ../pe_opt_2026-06-06
  (b) CPU peak per-rank memory               (LJ + PE)   <- ../pe_opt_2026-06-06
  (c) single-GPU throughput vs system size   (LJ)        <- ../gpu_scaling_optimized
  (d) multi-GPU throughput @ 64M atoms       (LJ)        <- ../gpu_scaling_optimized
CPU: EPYC, OMP=1, srun --mpi=pmi2, fair /usr/bin/time -v RSS for both codes.
GPU: GH200, AV single-process multi-GPU (NVLink) vs LAMMPS-KOKKOS/CUDA."""
import csv, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scienceplots  # noqa
plt.style.use(["science", "no-latex", "grid"])
matplotlib.rcParams.update({"font.size": 8.5, "axes.titlesize": 9.5,
    "axes.labelsize": 8.5, "xtick.labelsize": 7.0, "ytick.labelsize": 7.5,
    "legend.fontsize": 7.0, "lines.linewidth": 1.6, "lines.markersize": 5})
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
C_AV, C_LMP = "#0C5DA5", "#FF2C00"

def load(path, predicate=lambda r: True):
    with open(path) as fh:
        return [r for r in csv.DictReader(l for l in fh if not l.startswith("#"))
                if predicate(r)]

# ---- CPU multicore (LJ + PE) ----
CPU = load(os.path.join(RES, "pe_opt_2026-06-06", "lj_pe_weak_matched.csv"))
def cpu_series(system, code, col):
    d = {int(r["ranks"]): float(r[col]) for r in CPU
         if r["system"] == system and r["code"] == code and r[col] not in ("NA", "")}
    x = np.array(sorted(d)); return x, np.array([d[k] for k in x])

# ---- GPU single + multi (LJ) ----
GS_AV  = load(os.path.join(RES, "gpu_scaling_optimized", "av_single.csv"))
GS_LMP = load(os.path.join(RES, "gpu_scaling_optimized", "lmp_single.csv"))
GM_AV  = load(os.path.join(RES, "gpu_scaling_optimized", "av_multi.csv"))
GM_LMP = load(os.path.join(RES, "gpu_scaling_optimized", "lmp_multi.csv"))

fig, ((axA, axB), (axC, axD)) = plt.subplots(2, 2, figsize=(10.4, 8.4))
CORES = [1,2,4,8,16,32,64,128,256,512,1024]

# (a) CPU multicore weak throughput — LJ + PE
for sysn, mk, ls in (("lj", "o", "-"), ("pe", "^", "--")):
    x, av = cpu_series(sysn, "av", "matom_step_s")
    _, lm = cpu_series(sysn, "lammps", "matom_step_s")
    axA.plot(x, av, mk + ls, color=C_AV, label=f"AV ({sysn.upper()})")
    axA.plot(x, lm, mk + ls, color=C_LMP, mfc="none", label=f"LAMMPS ({sysn.upper()})")
axA.set_title("(a) CPU multicore — weak-scaling throughput")
axA.set_xlabel("MPI ranks (cores)"); axA.set_ylabel("Matom-step / s")
axA.set_xscale("log", base=2); axA.set_yscale("log")
axA.set_xticks(CORES); axA.set_xticklabels([str(c) for c in CORES], rotation=40)
axA.axvspan(256, 1100, color="0.85", alpha=0.5, zorder=0)
axA.legend(loc="upper left", ncol=2)

# (b) CPU peak memory per rank — LJ + PE
for sysn, mk, ls in (("lj", "o", "-"), ("pe", "^", "--")):
    x, av = cpu_series(sysn, "av", "maxrss_mb")
    _, lm = cpu_series(sysn, "lammps", "maxrss_mb")
    axB.plot(x, av, mk + ls, color=C_AV, label=f"AV ({sysn.upper()})")
    axB.plot(x, lm, mk + ls, color=C_LMP, mfc="none", label=f"LAMMPS ({sysn.upper()})")
axB.set_title("(b) CPU multicore — peak memory / rank")
axB.set_xlabel("MPI ranks (cores)"); axB.set_ylabel("max RSS / rank [MB]")
axB.set_xscale("log", base=2); axB.set_yscale("log")
axB.set_xticks(CORES); axB.set_xticklabels([str(c) for c in CORES], rotation=40)
axB.axvspan(256, 1100, color="0.85", alpha=0.5, zorder=0)
axB.legend(loc="upper left", ncol=2)

# (c) single-GPU throughput vs system size
xa = np.array([int(r["atoms"]) for r in GS_AV]);  ya = np.array([float(r["matom_step_s"]) for r in GS_AV])
xl = np.array([int(r["atoms"]) for r in GS_LMP]); yl = np.array([float(r["matom_step_s"]) for r in GS_LMP])
axC.plot(xa, ya, "o-", color=C_AV, label="Atomiverse")
axC.plot(xl, yl, "s-", color=C_LMP, label="LAMMPS-KOKKOS")
axC.set_title("(c) single-GPU throughput vs system size (GH200)")
axC.set_xlabel("atoms"); axC.set_ylabel("Matom-step / s")
axC.set_xscale("log")
for a, ya_, yl_ in zip(xa, ya, np.interp(xa, xl, yl)):
    if a in (1000000, 8000000, 32768000):
        axC.annotate(f"{ya_/yl_:.2f}x", (a, ya_), fontsize=6.5, color=C_AV,
                     xytext=(0, 6), textcoords="offset points", ha="center")
axC.legend(loc="lower right")

# (d) multi-GPU throughput @ 64M atoms
def at(rows, n): return float([r for r in rows if int(r["atoms"]) == 64000000
                               and int(r["gpus"]) == n][0]["matom_step_s"])
G = [1, 2, 4]
av = [at(GM_AV, n) for n in G]; lm = [at(GM_LMP, n) for n in G]
w = 0.35; xi = np.arange(len(G))
axD.bar(xi - w/2, av, w, color=C_AV, label="Atomiverse")
axD.bar(xi + w/2, lm, w, color=C_LMP, label="LAMMPS-KOKKOS")
for i, (a, l) in enumerate(zip(av, lm)):
    axD.text(i - w/2, a, f"{a/l:.1f}x", ha="center", va="bottom", fontsize=6.5, color=C_AV)
axD.set_title("(d) multi-GPU throughput @ 64M atoms (GH200, NVLink)")
axD.set_xlabel("GPUs"); axD.set_ylabel("Matom-step / s")
axD.set_xticks(xi); axD.set_xticklabels([str(g) for g in G])
axD.legend(loc="upper left")

fig.suptitle("Atomiverse vs LAMMPS — benchmark summary (2026-06-06).  "
             "CPU: AV leads at every core count (LJ 1.2-3.0x, PE 1.02-1.06x).  "
             "GPU: AV leads at large N single-GPU (~1.15x) and multi-GPU (up to 2.2x @4).",
             fontsize=8.8, y=1.005)
fig.tight_layout(rect=[0, 0, 1, 0.98])
out = os.path.join(HERE, "benchmark_summary.png")
fig.savefig(out, dpi=200, bbox_inches="tight"); fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
print("wrote", out)

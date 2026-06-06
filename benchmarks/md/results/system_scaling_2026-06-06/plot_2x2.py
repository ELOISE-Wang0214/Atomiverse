#!/usr/bin/env python3
"""Publication-grade 2x2 system benchmark: Atomiverse vs LAMMPS, one LJ system
(rho*=0.712, rc=2 sigma) across CPU and GPU.
  (a) CPU weak scaling   (16k atoms/core)  — speed + memory vs cores
  (b) CPU strong scaling (1M atoms fixed)  — speed + memory vs cores
  (c) 1x GH200           — speed + VRAM vs atom count
  (d) multi-GH200 strong scaling @ Pareto N — speed vs GPUs
CPU = double precision (EPYC, MPICH); GPU = mixed precision (KOKKOS/CUDA vs AV CUDA)."""
import csv, os, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scienceplots  # noqa
plt.style.use(["science", "no-latex", "grid"])
matplotlib.rcParams.update({"font.size": 8.5, "axes.titlesize": 9,
    "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7, "lines.linewidth": 1.5, "lines.markersize": 5})
# All CSVs are co-located with this script in the results directory.
SB = GB = os.path.dirname(os.path.abspath(__file__))
C_AV, C_LMP = "#0C5DA5", "#FF2C00"

def rows(path):
    try:
        with open(path) as f:
            return [r for r in csv.DictReader(f) if r.get("status", "ok") == "ok"]
    except FileNotFoundError:
        return []

def alab(x):
    x = float(x)
    return f"{x/1e3:g}k" if x < 1e6 else f"{x/1e6:g}M"

cpu = rows(f"{SB}/cpu_scaling.csv")
def cpu_xy(mode, code, col):
    d = {int(r["cores"]): float(r[col]) for r in cpu
         if r["mode"] == mode and r["code"] == code and r.get(col) not in ("","NA",None) and float(r[col])>0}
    x = np.array(sorted(d)); return x, np.array([d[k] for k in x])

avs, lms = rows(f"{GB}/av_single.csv"), rows(f"{GB}/lmp_single.csv")
avm, lmm = rows(f"{GB}/av_multi.csv"), rows(f"{GB}/lmp_multi.csv")
def gpu_single(rs, col):
    d = {int(r["atoms"]): float(r[col]) for r in rs if r.get(col) not in ("","NA",None) and float(r[col])>0}
    x = np.array(sorted(d)); return x, np.array([d[k] for k in x])

fig, axes = plt.subplots(2, 2, figsize=(10.2, 8.0))
(axA, axB), (axC, axD) = axes

def speed_mem_panel(ax, title, xlabel, series, xticks=None, xlog=True,
                    memcol=None, xtl=None, speedlog=False, memlog=True):
    """series: list of (x, speed, mem_or_None, color, label)."""
    for x, sp, mem, c, lab in series:
        ax.plot(x, sp, "o-", color=c, label=lab)
    ax.set_title(title); ax.set_xlabel(xlabel)
    ax.set_ylabel("throughput [Matom-step / s]" + (" (log)" if speedlog else ""))
    if xlog: ax.set_xscale("log")
    if speedlog: ax.set_yscale("log")
    if xticks is not None:
        ax.set_xticks(xticks); ax.set_xticklabels(xtl or [str(t) for t in xticks])
    axm = None
    if any(m is not None for _, _, m, _, _ in series):
        axm = ax.twinx(); axm.grid(False)
        for x, sp, mem, c, lab in series:
            if mem is not None:
                axm.plot(x, mem, "s--", color=c, alpha=0.45, ms=3.5, lw=1.0)
        if memlog: axm.set_yscale("log")
        axm.set_ylabel(memcol + (" (log, dashed)" if memlog else " (dashed)"))
    return axm

# (a) CPU weak
cores = sorted({int(r["cores"]) for r in cpu if r["mode"] == "weak"}) or [1]
aw_s = cpu_xy("weak", "av", "matom_step_s"); lw_s = cpu_xy("weak", "lammps", "matom_step_s")
aw_m = cpu_xy("weak", "av", "peak_rss_kb"); lw_m = cpu_xy("weak", "lammps", "peak_rss_kb")
speed_mem_panel(axA, "(a) CPU weak scaling — 16k atoms/core (EPYC)", "MPI ranks (cores)",
    [(aw_s[0], aw_s[1], aw_m[1]/1e6, C_AV, "Atomiverse"),
     (lw_s[0], lw_s[1], lw_m[1]/1e6, C_LMP, "LAMMPS")],
    xticks=cores, xtl=[str(c) for c in cores], memcol="peak RSS [GB]",
    speedlog=True)
axA.legend(loc="upper left", framealpha=0.9)

# (b) CPU strong
as_s = cpu_xy("strong", "av", "matom_step_s"); ls_s = cpu_xy("strong", "lammps", "matom_step_s")
as_m = cpu_xy("strong", "av", "peak_rss_kb"); ls_m = cpu_xy("strong", "lammps", "peak_rss_kb")
speed_mem_panel(axB, "(b) CPU strong scaling — 1M atoms fixed (EPYC)", "MPI ranks (cores)",
    [(as_s[0], as_s[1], as_m[1]/1e6, C_AV, "Atomiverse"),
     (ls_s[0], ls_s[1], ls_m[1]/1e6, C_LMP, "LAMMPS")],
    xticks=cores, xtl=[str(c) for c in cores], memcol="peak RSS [GB]",
    speedlog=True)
axB.legend(loc="upper left", framealpha=0.9)

# (c) GPU single vs atoms
ca = gpu_single(avs, "matom_step_s"); cl = gpu_single(lms, "matom_step_s")
cav = gpu_single(avs, "peak_vram_mb"); clv = gpu_single(lms, "peak_vram_mb")
xt = ca[0] if len(ca) else [1]
speed_mem_panel(axC, "(c) Single GH200 — vs system size (mixed prec.)", "atom count",
    [(ca[0], ca[1], cav[1]/1e3, C_AV, "Atomiverse CUDA"),
     (cl[0], cl[1], clv[1]/1e3, C_LMP, "LAMMPS KOKKOS")],
    xticks=xt, xtl=[alab(x) for x in xt], memcol="peak VRAM [GB]",
    speedlog=True)
for t in axC.get_xticklabels(): t.set_rotation(35)
axC.legend(loc="upper left", framealpha=0.9)

# (d) GPU multi strong scaling @ Pareto N (max single-GPU AV throughput)
mcand = sorted({int(r["atoms"]) for r in avm})
aspd = {int(r["atoms"]): float(r["matom_step_s"]) for r in avs}
PN = max(mcand, key=lambda x: aspd.get(x, 0)) if mcand else 0
def scal(rs, n):
    d = {int(r["gpus"]): float(r["matom_step_s"]) for r in rs
         if int(r["atoms"]) == n and r["matom_step_s"] not in ("", "NA")}
    g = np.array(sorted(d)); return g, np.array([d[k] for k in g])
ga, sa = scal(avm, PN); gl, sl = scal(lmm, PN)
axD.plot(ga, sa, "o-", color=C_AV, label="Atomiverse (CUDA P2P)")
axD.plot(gl, sl, "s-", color=C_LMP, label="LAMMPS KOKKOS (MPI)")
if len(gl): axD.plot(gl, sl[0]*gl, ":", color="0.5", lw=1.0, label="ideal (LAMMPS 1-GPU x N)")
axD.set_title(f"(d) Multi-GH200 strong scaling @ {alab(PN)} atoms"); axD.set_xlabel("GPUs")
axD.set_ylabel("throughput [Matom-step / s]"); axD.set_xticks([1,2,4])
axD.legend(loc="upper left", framealpha=0.9)

fig.suptitle("Atomiverse vs LAMMPS — one LJ system (rho*=0.71, rc=2 sigma) on GH200 node "
             "(CPU: EPYC double · GPU: mixed)", fontsize=9.5, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = f"{SB}/sysbench_2x2.png"
fig.savefig(out, dpi=200); fig.savefig(out.replace(".png", ".pdf"))
print("Pareto N =", PN); print("wrote", out)

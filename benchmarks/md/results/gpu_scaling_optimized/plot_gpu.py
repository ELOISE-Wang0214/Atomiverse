#!/usr/bin/env python3
"""GPU LJ benchmark — optimized Atomiverse CUDA vs LAMMPS KOKKOS/CUDA on GH200.

Produces TWO separate figures, each overlaying throughput (solid, left axis) and
peak VRAM (dashed, right axis) for BOTH codes on one plot:
  gpu_single.png  — single GH200, vs atom count.
  gpu_multi.png   — multi-GH200 WEAK scaling (fixed per-card load), vs GPU count.
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    import scienceplots  # noqa
    plt.style.use(["science", "no-latex", "grid"])
except Exception:
    plt.style.use("seaborn-v0_8-whitegrid")
matplotlib.rcParams.update({"font.size": 9, "axes.titlesize": 10.5,
    "axes.labelsize": 9.5, "legend.fontsize": 8, "lines.linewidth": 1.7,
    "lines.markersize": 6})
HERE = os.path.dirname(os.path.abspath(__file__))
C_AV, C_LMP = "#0C5DA5", "#FF2C00"   # blue = Atomiverse, red = LAMMPS

def rows(path):
    try:
        with open(os.path.join(HERE, path)) as f:
            return [r for r in csv.DictReader(f) if r.get("status", "ok") == "ok"]
    except FileNotFoundError:
        return []

def col(rs, key, val):
    d = {float(r[key]): float(r[val]) for r in rs
         if r.get(val) not in ("", "NA", None) and float(r[val]) > 0}
    x = np.array(sorted(d)); return x, np.array([d[k] for k in x])

def alab(x):
    x = float(x)
    return f"{x/1e3:g}k" if x < 1e6 else f"{x/1e6:g}M"

# ───────────────────────── figure 1: single GH200 ─────────────────────────
avs, lms = rows("av_single.csv"), rows("lmp_single.csv")
fig, ax = plt.subplots(figsize=(7.2, 5.0))
ax2 = ax.twinx(); ax2.grid(False)

ax_x, ax_t = col(avs, "atoms", "matom_step_s")
lx_x, lx_t = col(lms, "atoms", "matom_step_s")
l1, = ax.plot(ax_x, ax_t, "o-", color=C_AV)
l2, = ax.plot(lx_x, lx_t, "s-", color=C_LMP)
av_vx, av_v = col(avs, "atoms", "peak_vram_mb")
lv_vx, lv_v = col(lms, "atoms", "peak_vram_mb")
l3, = ax2.plot(av_vx, av_v/1e3, "o--", color=C_AV, alpha=0.5, ms=4, lw=1.1)
l4, = ax2.plot(lv_vx, lv_v/1e3, "s--", color=C_LMP, alpha=0.5, ms=4, lw=1.1)

ax.set_xscale("log"); ax.set_yscale("log"); ax2.set_yscale("log")
ax.set_xlabel("atom count")
ax.set_ylabel("throughput [Matom-step/s]  (solid)")
ax2.set_ylabel("peak VRAM [GB]  (dashed)")
ax.set_title("Single GH200 — LJ throughput & VRAM vs system size (mixed precision)")
xt = ax_x if len(ax_x) else lx_x
ax.set_xticks(xt); ax.set_xticklabels([alab(x) for x in xt], rotation=35)
ax.legend([l1, l2, l3, l4],
          ["Atomiverse throughput", "LAMMPS throughput",
           "Atomiverse VRAM", "LAMMPS VRAM"],
          loc="upper left", framealpha=0.92, ncol=2)
fig.tight_layout()
fig.savefig(f"{HERE}/gpu_single.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{HERE}/gpu_single.pdf", bbox_inches="tight")
plt.close(fig)

# ───────────────────── figure 2: multi GH200 weak scaling ─────────────────
avw, lmw = rows("av_weak.csv"), rows("lmp_weak.csv")
per_card = int(float(avw[0]["atoms_per_gpu"])) if avw else 0
fig, ax = plt.subplots(figsize=(7.2, 5.0))
ax2 = ax.twinx(); ax2.grid(False)

gx_a, t_a = col(avw, "gpus", "matom_step_s")     # total throughput
gx_l, t_l = col(lmw, "gpus", "matom_step_s")
l1, = ax.plot(gx_a, t_a, "o-", color=C_AV)
l2, = ax.plot(gx_l, t_l, "s-", color=C_LMP)
lid = None
if len(gx_a):
    lid, = ax.plot(gx_a, t_a[0]*gx_a, ":", color="0.5", lw=1.2)
# per-GPU VRAM: AV peak_vram_mb is already per-GPU (busiest card); LAMMPS
# peak_vram_mb is summed over the GPUs used, so divide by the GPU count.
va_x, va = col(avw, "gpus", "peak_vram_mb")
vl_x, vl_raw = col(lmw, "gpus", "peak_vram_mb")
vl = vl_raw / vl_x if len(vl_x) else vl_raw
l3, = ax2.plot(va_x, va/1e3, "o--", color=C_AV, alpha=0.5, ms=4, lw=1.1)
l4, = ax2.plot(vl_x, vl/1e3, "s--", color=C_LMP, alpha=0.5, ms=4, lw=1.1)

ax.set_xlabel("GPUs  (weak scaling: total atoms = GPUs × %s)" % alab(per_card))
ax.set_ylabel("total throughput [Matom-step/s]  (solid)")
ax2.set_ylabel("peak VRAM per GPU [GB]  (dashed)")
ax.set_title("Multi-GH200 weak scaling — %s atoms/GPU (NVLink, mixed precision)"
             % alab(per_card))
ax.set_xticks([1, 2, 4])
handles = [l1, l2] + ([lid] if lid is not None else []) + [l3, l4]
labels = ["Atomiverse throughput", "LAMMPS throughput"] \
    + (["ideal (AV 1-GPU × N)"] if lid is not None else []) \
    + ["Atomiverse VRAM/GPU", "LAMMPS VRAM/GPU"]
ax.legend(handles, labels, loc="upper left", framealpha=0.92, ncol=2)
fig.tight_layout()
fig.savefig(f"{HERE}/gpu_multi.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{HERE}/gpu_multi.pdf", bbox_inches="tight")
plt.close(fig)

print("wrote gpu_single.png and gpu_multi.png  (per-card =", per_card, ")")

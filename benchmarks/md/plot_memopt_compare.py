#!/usr/bin/env python3
"""Before/after comparison of the per-atom memory optimizations (AV only).
Two pinned runs (baseline vs optimized); 2x2 grid:
  row 0 = memory (left axis): GPU VRAM [GB] / CPU RSS [GB]
  row 1 = throughput [Matom·step/s]
  cols  = LJ, PE.  weak+strong both drawn; baseline grey, optimized blue.
Mean ± sample-std over reps. Annotates the mean % change at each point.

  python plot_memopt_compare.py --old c73d896_2026-06-08 --new 968d036_2026-06-08 \
         --backend gpu --out memopt_gpu.png
"""
import argparse, csv, glob, os, collections, statistics
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BD = os.path.dirname(os.path.abspath(__file__))

def load(system, run_id, backend):
    path = os.path.join(BD, system, "runs", run_id, f"results_{backend}.csv")
    acc = collections.defaultdict(lambda: collections.defaultdict(list))
    if os.path.exists(path):
        for r in csv.DictReader(open(path)):
            if r["engine"] != "AV":
                continue
            k = (r["mode"], int(r["ng"]))
            acc[k]["thr"].append(float(r["matom_step_s"]))
            acc[k]["mem"].append(float(r["vram_gb"]))
    return acc

def agg(xs):
    m = statistics.mean(xs) if xs else float("nan")
    s = statistics.stdev(xs) if len(xs) > 1 else 0.0
    return m, s

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--backend", choices=["gpu", "cpu"], default="gpu")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    xlab = "GPUs" if a.backend == "gpu" else "MPI ranks (cores)"
    mlab = "peak VRAM [GB]" if a.backend == "gpu" else "peak RSS [GB]"
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    OLD = ("#999999", f"baseline ({a.old.split('_')[0]})")
    NEW = ("#1f77b4", f"optimized ({a.new.split('_')[0]})")
    for col, system in enumerate(["LJ", "PE"]):
        do = load(system.lower(), a.old, a.backend)
        dn = load(system.lower(), a.new, a.backend)
        for row, metric in enumerate(["mem", "thr"]):
            ax = axes[row, col]
            for run, (color, label) in [(do, OLD), (dn, NEW)]:
                for mode, ls, mk in [("weak", "-", "o"), ("strong", "--", "s")]:
                    ks = sorted(k for k in run if k[0] == mode)
                    if not ks:
                        continue
                    xs = [k[1] for k in ks]
                    ys = [agg(run[k][metric])[0] for k in ks]
                    es = [agg(run[k][metric])[1] for k in ks]
                    ax.errorbar(xs, ys, yerr=es, fmt=mk + ls, color=color,
                                lw=2, ms=6, capsize=3,
                                label=f"{label} · {mode}")
            # annotate % change (optimized vs baseline) per mode/point
            for mode in ["weak", "strong"]:
                for k in sorted(x for x in dn if x[0] == mode):
                    if k not in do:
                        continue
                    o = agg(do[k][metric])[0]; n = agg(dn[k][metric])[0]
                    if not o or o != o:
                        continue
                    pct = (n - o) / o * 100
                    ax.annotate(f"{pct:+.0f}%", (k[1], n),
                                textcoords="offset points", xytext=(4, 5),
                                fontsize=7, color="#1f77b4")
            ax.set_title(f"{system} · {'memory' if row==0 else 'throughput'}",
                         fontsize=12, fontweight="bold")
            ax.set_xlabel(xlab)
            ax.set_ylabel(mlab if row == 0 else "throughput [Matom·step/s]")
            ax.set_ylim(bottom=0); ax.grid(True, alpha=0.3)
            allx = sorted({k[1] for k in dn})
            if allx:
                ax.set_xticks(allx)
            if row == 0 and col == 0:
                ax.legend(fontsize=7.5, loc="upper left")
    fig.suptitle(
        f"AV per-atom memory optimization — {a.backend.upper()} "
        f"({a.old.split('_')[0]} → {a.new.split('_')[0]})\n"
        "bins→bin_pos.w + vl_max 116 (GPU) · lazy bonded storage (CPU); "
        "mean ± std, 3 reps", fontsize=12.5, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(a.out, dpi=130, bbox_inches="tight")
    print("wrote", a.out)

if __name__ == "__main__":
    main()

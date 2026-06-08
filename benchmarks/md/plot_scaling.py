#!/usr/bin/env python3
"""2x2 scaling figure from per-system run dirs (bench_scaling*.sbatch).
Rows = system (LJ, PE); cols = weak, strong.
  bottom x-axis = device count;  TOP x-axis = total atoms at that point
  left y-axis   = throughput [Matom-step/s] (AV vs LAMMPS, solid, value-labelled)
  right y-axis  = peak memory [GB] (dashed; VRAM for gpu, RSS for cpu)

  python plot_scaling.py --backend gpu     # reads <sys>/runs/<latest>/results_gpu.csv
  python plot_scaling.py --backend cpu
  python plot_scaling.py --backend gpu --id <sha>_<date> --out f.png
"""
import argparse, csv, glob, os, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BD = os.path.dirname(os.path.abspath(__file__))
SYS_DESC = {"LJ": "Lennard-Jones (argon)", "PE": "polyethylene melt (UA)"}

def latest_run(system, backend):
    g = sorted(d for d in glob.glob(os.path.join(BD, system, "runs", "*"))
               if os.path.exists(os.path.join(d, f"results_{backend}.csv")))
    return os.path.basename(g[-1]) if g else None

def load(system, run_id, backend):
    path = os.path.join(BD, system, "runs", run_id, f"results_{backend}.csv")
    d = collections.defaultdict(dict)  # (mode,engine)->{n:(thrpt,mem,atoms)}
    if run_id and os.path.exists(path):
        for r in csv.DictReader(open(path)):
            try:
                d[(r["mode"], r["engine"])][int(r["ng"])] = (
                    float(r["matom_step_s"]), float(r["vram_gb"]), int(r["atoms"]))
            except (ValueError, KeyError):
                pass
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=["gpu", "cpu"], default="gpu")
    ap.add_argument("--id", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    xlab = "GPUs" if a.backend == "gpu" else "MPI ranks (cores)"
    mlab = "peak VRAM [GB]" if a.backend == "gpu" else "peak RSS [GB]"
    out = a.out or os.path.join(BD, f"scaling_2x2_{a.backend}.png")
    ENG = {"AV": ("#1f77b4", "AV"), "LAMMPS": ("#d62728", "LAMMPS-KOKKOS")}
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ids = {}
    for i, system in enumerate(["LJ", "PE"]):
        rid = a.id or latest_run(system.lower(), a.backend)
        ids[system] = rid
        data = load(system.lower(), rid, a.backend)
        for j, mode in enumerate(["weak", "strong"]):
            ax = axes[i, j]; axr = ax.twinx()
            ref = data.get((mode, "AV")) or data.get((mode, "LAMMPS")) or {}
            for eng, (color, label) in ENG.items():
                dd = data.get((mode, eng), {})
                if not dd:
                    continue
                xs = sorted(dd); ys = [dd[x][0] for x in xs]
                ax.plot(xs, ys, "o-", color=color, lw=2.2, ms=7, label=label)
                for x, y in zip(xs, ys):
                    ax.annotate(f"{y:.0f}", (x, y), textcoords="offset points",
                                xytext=(0, 7), ha="center", fontsize=7.5, color=color)
                axr.plot(xs, [dd[x][1] for x in xs], "s--", color=color, lw=1.2,
                         ms=4, alpha=0.5)
            ax.set_title(f"{system} · {mode} scaling", fontsize=12, fontweight="bold")
            ax.set_xlabel(xlab); ax.set_ylabel("throughput [Matom·step/s]")
            ax.set_ylim(bottom=0); ax.grid(True, alpha=0.3)
            axr.set_ylabel(mlab, color="gray"); axr.set_ylim(bottom=0)
            axr.tick_params(axis="y", colors="gray")
            if ref:
                xs = sorted(ref); ax.set_xticks(xs)
                sec = ax.secondary_xaxis("top")
                sec.set_xticks(xs)
                sec.set_xticklabels([f"{ref[n][2]/1e6:.0f}M" for n in xs])
                sec.set_xlabel("total atoms", fontsize=9); sec.tick_params(labelsize=7.5)
            if j == 0:
                ax.text(-0.20, 0.5, f"{system}\n{SYS_DESC[system]}", transform=ax.transAxes,
                        rotation=90, va="center", ha="center", fontsize=11, fontweight="bold")
            if i == 0 and j == 0:
                ax.plot([], [], "ko-", label="throughput (L)")
                ax.plot([], [], "ks--", alpha=0.5, label="memory (R)")
                ax.legend(fontsize=8, loc="upper left")
    sub = "  |  ".join(f"{s}: {ids[s]}" for s in ids if ids[s])
    fig.suptitle(f"AV vs LAMMPS-KOKKOS — {a.backend.upper()} scaling (GH200)\n{sub}",
                 fontsize=13, y=1.0)
    fig.tight_layout(rect=[0.03, 0, 1, 1]); fig.savefig(out, dpi=130, bbox_inches="tight")
    print("wrote", out, "from", ids)

if __name__ == "__main__":
    main()

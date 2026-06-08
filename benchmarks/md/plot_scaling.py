#!/usr/bin/env python3
"""Plot the 2x2 GPU scaling figure from per-system run dirs written by
bench_scaling.sbatch. Rows = system (LJ, PE); cols = weak, strong. Left y-axis =
throughput [Matom-step/s] (AV vs LAMMPS, solid); right y-axis = peak VRAM [GB]
(dashed). x = GPU count.

  python plot_scaling.py                  # latest run per system -> scaling_2x2.png
  python plot_scaling.py --id <sha>_<date>  # pin a specific run id (both systems)
  python plot_scaling.py --out foo.png

The figure is NOT committed (only results.csv/meta.json are). Re-generate anytime.
"""
import argparse, csv, glob, os, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BD = os.path.dirname(os.path.abspath(__file__))

def latest_run(system):
    dirs = sorted(glob.glob(os.path.join(BD, system, "runs", "*")))
    return os.path.basename(dirs[-1]) if dirs else None

def load(system, run_id):
    path = os.path.join(BD, system, "runs", run_id, "results.csv")
    if not os.path.exists(path):
        return {}
    d = collections.defaultdict(dict)  # (mode,engine) -> {ng:(mas,vram)}
    for r in csv.DictReader(open(path)):
        try:
            d[(r["mode"], r["engine"])][int(r["ng"])] = (float(r["matom_step_s"]), float(r["vram_gb"]))
        except (ValueError, KeyError):
            pass
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", default=None, help="run id <sha>_<date> (default: latest per system)")
    ap.add_argument("--out", default=os.path.join(BD, "scaling_2x2.png"))
    a = ap.parse_args()
    ENG = {"AV": ("#1f77b4", "AV"), "LAMMPS": ("#d62728", "LAMMPS-KOKKOS")}
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    ids = {}
    for i, system in enumerate(["LJ", "PE"]):
        sysl = system.lower()
        rid = a.id or latest_run(sysl)
        ids[system] = rid
        data = load(sysl, rid) if rid else {}
        for j, mode in enumerate(["weak", "strong"]):
            ax = axes[i, j]; axr = ax.twinx()
            for eng, (color, label) in ENG.items():
                dd = data.get((mode, eng), {})
                if not dd:
                    continue
                xs = sorted(dd)
                ax.plot(xs, [dd[x][0] for x in xs], "o-", color=color, lw=2, ms=7, label=f"{label} thrpt")
                axr.plot(xs, [dd[x][1] for x in xs], "s--", color=color, lw=1.3, ms=5, alpha=0.55, label=f"{label} VRAM")
            ax.set_title(f"{system} · {mode} scaling", fontsize=12, fontweight="bold")
            ax.set_xlabel("GPUs"); ax.set_xticks([1, 2, 3, 4]); ax.set_ylim(bottom=0)
            ax.set_ylabel("throughput [Matom·step/s]")
            axr.set_ylabel("peak VRAM [GB]", color="gray"); axr.set_ylim(bottom=0)
            axr.tick_params(axis="y", colors="gray"); ax.grid(True, alpha=0.3)
            if i == 0 and j == 0:
                h1, l1 = ax.get_legend_handles_labels(); h2, l2 = axr.get_legend_handles_labels()
                ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left")
    sub = " | ".join(f"{s}:{ids[s]}" for s in ids if ids[s])
    fig.suptitle(f"AV vs LAMMPS-KOKKOS — GPU scaling (GH200)   [{sub}]", fontsize=12, y=1.0)
    fig.tight_layout(); fig.savefig(a.out, dpi=130, bbox_inches="tight")
    print("wrote", a.out, "from", ids)

if __name__ == "__main__":
    main()

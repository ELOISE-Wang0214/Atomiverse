#!/usr/bin/env python3
"""2x2 scaling figure from per-system run dirs (bench_scaling*.sbatch).
Rows = system (LJ, PE); cols = weak, strong.

PLOT CONVENTION (fixed — do not change; mirrored in bench_scaling.README.md):
  bottom x-axis = device count;  TOP x-axis = total atoms at that point
  left y-axis   = throughput [Matom-step/s] (AV vs LAMMPS, solid, value-labelled)
  right y-axis  = peak memory [GB] (dashed; VRAM for gpu, RSS for cpu)

  python plot_scaling.py --backend gpu     # reads <sys>/performance/runs/<latest>/results_gpu.csv
  python plot_scaling.py --backend cpu
  python plot_scaling.py --backend gpu --id <sha>_<date> --out f.png
"""
import argparse, csv, glob, os, collections, statistics
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# This script lives in benchmarks/md/common/; systems live one level up
# at benchmarks/md/<sys>/performance/runs/.
BD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SYS_DESC = {"LJ": "Lennard-Jones (argon)", "PE": "polyethylene melt (UA)"}

def latest_run(system, backend):
    g = sorted(d for d in glob.glob(os.path.join(BD, system, "performance", "runs", "*"))
               if os.path.exists(os.path.join(d, f"results_{backend}.csv")))
    return os.path.basename(g[-1]) if g else None

def fmt_atoms(a):
    """Adaptive atom-count label: K below 1 M, else M (CPU runs are far smaller)."""
    return f"{a/1e6:.0f}M" if a >= 1e6 else f"{a/1e3:.0f}K"

def _agg(vals):
    """(mean, sample-std) of a list; std=0.0 for a single rep."""
    m = statistics.mean(vals)
    s = statistics.stdev(vals) if len(vals) > 1 else 0.0
    return m, s

def load(system, run_id, backend):
    """Aggregate repeats: rows sharing (mode,engine,ng) are pooled into
    mean ± sample-std. Schema is rep-aware but back-compatible — a CSV with no
    `rep` column (one row per key) yields std 0."""
    path = os.path.join(BD, system, "performance", "runs", run_id, f"results_{backend}.csv")
    acc = collections.defaultdict(lambda: collections.defaultdict(
        lambda: {"thr": [], "mem": [], "atoms": None}))  # (mode,eng)->n->lists
    if run_id and os.path.exists(path):
        for r in csv.DictReader(open(path)):
            try:
                cell = acc[(r["mode"], r["engine"])][int(r["ng"])]
                cell["thr"].append(float(r["matom_step_s"]))
                cell["mem"].append(float(r["vram_gb"]))
                cell["atoms"] = int(r["atoms"])
            except (ValueError, KeyError):
                pass
    d = collections.defaultdict(dict)  # (mode,eng)->{n:(thr,thr_std,mem,mem_std,atoms,reps)}
    for key, byn in acc.items():
        for n, cell in byn.items():
            tm, ts = _agg(cell["thr"]); mm, ms = _agg(cell["mem"])
            d[key][n] = (tm, ts, mm, ms, cell["atoms"], len(cell["thr"]))
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=["gpu", "cpu"], default="gpu")
    ap.add_argument("--id", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--md", default=None,
                    help="also write a mean±std markdown table to this path")
    a = ap.parse_args()
    xlab = "GPUs" if a.backend == "gpu" else "MPI ranks (cores)"
    mlab = "peak VRAM [GB]" if a.backend == "gpu" else "peak RSS [GB]"
    out = a.out or os.path.join(BD, f"scaling_2x2_{a.backend}.png")
    ENG = {"AV": ("#1f77b4", "AV"), "LAMMPS": ("#d62728", "LAMMPS-KOKKOS")}
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ids = {}; alldata = {}
    for i, system in enumerate(["LJ", "PE"]):
        rid = a.id or latest_run(system.lower(), a.backend)
        ids[system] = rid
        data = load(system.lower(), rid, a.backend); alldata[system] = data
        for j, mode in enumerate(["weak", "strong"]):
            ax = axes[i, j]; axr = ax.twinx()
            ref = data.get((mode, "AV")) or data.get((mode, "LAMMPS")) or {}
            for eng, (color, label) in ENG.items():
                dd = data.get((mode, eng), {})
                if not dd:
                    continue
                xs = sorted(dd)
                ys = [dd[x][0] for x in xs]; yerr = [dd[x][1] for x in xs]
                ax.errorbar(xs, ys, yerr=yerr, fmt="o-", color=color, lw=2.2, ms=7,
                            capsize=4, elinewidth=1.3, label=label)
                for x, y in zip(xs, ys):
                    ax.annotate(f"{y:.0f}", (x, y), textcoords="offset points",
                                xytext=(5, 6), ha="left", fontsize=7.5, color=color)
                axr.errorbar(xs, [dd[x][2] for x in xs], yerr=[dd[x][3] for x in xs],
                             fmt="s--", color=color, lw=1.2, ms=4, alpha=0.5, capsize=3)
            ax.set_title(f"{system} · {mode} scaling", fontsize=12, fontweight="bold")
            ax.set_xlabel(xlab); ax.set_ylabel("throughput [Matom·step/s]")
            # headroom so value labels + top error caps stay inside the panel
            top = max((dd[x][0] + dd[x][1] for eng in ENG
                       for dd in [data.get((mode, eng), {})] for x in dd), default=0)
            ax.set_ylim(bottom=0, top=top * 1.22 if top else None)
            ax.grid(True, alpha=0.3)
            axr.set_ylabel(mlab, color="gray"); axr.set_ylim(bottom=0)
            axr.tick_params(axis="y", colors="gray")
            if ref:
                xs = sorted(ref)
                # Wide rank sweeps (CPU 1..1024) span >1 decade — a log x-axis
                # keeps the low end legible; GPU 1..4 stays linear.
                if a.backend == "cpu" and xs and xs[-1] / max(xs[0], 1) >= 16:
                    import matplotlib.ticker as mt
                    ax.set_xscale("log", base=2)
                    ax.xaxis.set_major_formatter(mt.ScalarFormatter())
                    ax.xaxis.set_minor_locator(mt.NullLocator())
                ax.set_xticks(xs); ax.set_xticklabels([str(x) for x in xs])
                sec = ax.secondary_xaxis("top")
                sec.set_xticks(xs)
                sec.set_xticklabels([fmt_atoms(ref[n][4]) for n in xs])
                sec.set_xlabel("total atoms", fontsize=9); sec.tick_params(labelsize=7.5)
            if j == 0:
                ax.text(-0.20, 0.5, f"{system}\n{SYS_DESC[system]}", transform=ax.transAxes,
                        rotation=90, va="center", ha="center", fontsize=11, fontweight="bold")
            if i == 0 and j == 0:
                ax.plot([], [], "ko-", label="throughput (L)")
                ax.plot([], [], "ks--", alpha=0.5, label="memory (R)")
                ax.legend(fontsize=8, loc="upper left")
    sub = "  |  ".join(f"{s}: {ids[s]}" for s in ids if ids[s])
    hw = "GH200" if a.backend == "gpu" else "MPI ranks, up to 4 nodes"
    fig.suptitle(f"AV vs LAMMPS-KOKKOS — {a.backend.upper()} scaling ({hw})\n{sub}",
                 fontsize=13, y=1.0)
    fig.tight_layout(rect=[0.03, 0, 1, 1]); fig.savefig(out, dpi=130, bbox_inches="tight")
    print("wrote", out, "from", ids)
    if a.md:
        write_markdown(a.md, alldata, ids, a.backend, xlab)
        print("wrote", a.md)


def write_markdown(path, alldata, ids, backend, xlab):
    """Mean ± sample-std table over reps, one block per system, AV vs LAMMPS."""
    unit = xlab.split()[0]
    lines = [f"### {backend.upper()} scaling — throughput (Matom·step/s), mean ± std\n"]
    for system in ["LJ", "PE"]:
        data = alldata.get(system, {})
        if not data:
            continue
        reps = max((c[5] for byn in data.values() for c in byn.values()), default=0)
        lines.append(f"**{system} — {SYS_DESC[system]}**  (run `{ids.get(system)}`, "
                     f"n={reps} reps)\n")
        for mode in ["weak", "strong"]:
            av = data.get((mode, "AV"), {}); lm = data.get((mode, "LAMMPS"), {})
            ns = sorted(set(av) | set(lm))
            if not ns:
                continue
            lines.append(f"_{mode} scaling_\n")
            lines.append(f"| {unit} | atoms | AV | LAMMPS-KOKKOS | AV/LAMMPS |")
            lines.append("|---:|---:|---:|---:|---:|")
            for n in ns:
                a_ = av.get(n); l_ = lm.get(n)
                atoms = (a_ or l_)[4]
                ac = f"{a_[0]:.0f} ± {a_[1]:.0f}" if a_ else "—"
                lc = f"{l_[0]:.0f} ± {l_[1]:.0f}" if l_ else "—"
                ratio = f"{a_[0]/l_[0]:.2f}×" if a_ and l_ and l_[0] else "—"
                lines.append(f"| {n} | {fmt_atoms(atoms)} | {ac} | {lc} | {ratio} |")
            lines.append("")
    open(path, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

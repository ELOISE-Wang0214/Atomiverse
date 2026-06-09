# Multi-GPU scaling benchmark — AV vs LAMMPS-KOKKOS

Permanent, auditable, comparable. **One** scan script + **one** plot script;
data is stored per system and keyed to the exact code that produced it. Never
copy these scripts into a dated folder — re-run them.

## Layout (organized by system, not by date)
```
benchmarks/md/
├── bench_scaling.sbatch     # the ONE scan (run per system)
├── plot_scaling.py          # the ONE plotter
├── lib_scaling.sh           # shared helpers (peak-VRAM, cbrt, LAMMPS loop-time parse)
├── lj/  in.lj_kk  config.sh  runs/<sha>_<date>/{results.csv, meta.json}
└── pe/  in.pe_kk  pe_base.data  gen_pe_chain_mol.py  config.sh  runs/<sha>_<date>/{...}
```
Committed: `results.csv` + `meta.json` per run (data only — an auditable history
of numbers vs code). Gitignored: `logs/`, `*.png`, `.scanout/` (figures are
regenerated, never committed).

## Run (each run builds the current committed HEAD first)
Each cell is measured `REPS` times (default **3**) so the plot can show error bars.
```bash
sbatch benchmarks/md/bench_scaling.sbatch lj      # -> lj/runs/<sha>_<date>/  (3 reps/cell)
sbatch benchmarks/md/bench_scaling.sbatch pe      # -> pe/runs/<sha>_<date>/
REPS=5 sbatch benchmarks/md/bench_scaling.sbatch lj   # more reps -> tighter error bars
# plot mean ± std error bars (latest run per system) + emit a markdown table:
python benchmarks/md/plot_scaling.py --out scaling_2x2_gpu_errorbar.png --md table.md
python benchmarks/md/plot_scaling.py --id <sha>_<date>   # pin a specific run
```

**CPU (1–1024 cores):** the simplified self-contained submit script —
```bash
cd <repo>                                                    # submit from inside the repo
LMP=/path/to/lmp sbatch benchmarks/md/bench_scaling_cpu.sbatch          # LJ (omit LMP ⇒ AV only)
SYS=pe LMP=/path/to/lmp sbatch benchmarks/md/bench_scaling_cpu.sbatch   # PE
# overrides: CORES="1 8 32 ..."  STEPS=  WEAK_PER_RANK=  STRONG_TOTAL=
# account/partition:  sbatch -A <proj> -p <part> ...
```
It writes the same `runs/<sha>_<date>/results_cpu.csv` layout and **auto-plots**
at the end (then prints the `plot_scaling.py` command to re-plot on the login node).

## Plot convention (every scaling figure — `plot_scaling.py`)
**Always** rendered the same way, do not change:
- **left y-axis** = throughput [Matom·step/s] (solid, value-labelled), AV vs LAMMPS.
- **right y-axis** = peak memory (dashed) — RSS [GB] for CPU, VRAM [GB] for GPU.
- **bottom x-axis** = device count (MPI ranks / GPUs); **top x-axis** = total atoms
  at that point.

`plot_scaling.py` implements this (`twinx()` for the right memory axis,
`secondary_xaxis("top")` for the atom-count axis); the CPU sbatch invokes it on
every run. Any new plotter must follow the same convention.

## Why this is auditable + comparable (可审计 / 可环比)
- **One code set:** `bench_scaling.sbatch`/`plot_scaling.py`/`lib_scaling.sh` are
  committed once and reused; no per-date duplicates.
- **Builds committed HEAD** and writes the git SHA into `meta.json`; the run dir
  is `runs/<sha>_<date>/` → every number traces to exact code.
- **`config.sh` per system** is the only place sizes/params live → no drift.
- **`results.csv` schema is fixed:** `engine,system,mode,ng,atoms,matom_step_s,vram_gb,rep`
  (`rep` = 1…`REPS`; rows sharing `engine,system,mode,ng` are repeats pooled into mean ± std).
- Compare two runs (before/after an optimization) by diffing their `results.csv`,
  or `plot_scaling.py --id` each.

## Methodology (held identical across runs)
- **GH200, 1–4 GPUs, single node.** weak = fixed per-GPU N* (LJ 8M, PE 4M);
  strong = fixed total (LJ 32M, PE 96³≈16M). Per-GPU N* from the Phase-0
  single-card saturation sweep.
- **LJ:** AV `bench_multi_gpu_scaling` (sc ρ*=0.712, rc=2σ, reneigh 20) vs LAMMPS
  `lj/in.lj_kk` (`lattice sc 0.712`, `lj/cut 2.0`).
- **PE:** AV `polyethylene_tension --pe-melt-grid` (lj/cut 10.5 + harmonic
  bond/angle + multi/harmonic dihedral, special_bonds 0 0 0, reneigh 10) vs
  LAMMPS `pe/in.pe_kk` (`read_data pe_base.data` + `replicate nx ny nz`, **same
  grid as AV** → identical structures), MSU Tschopp 2010 Dreiding FF.
- AV float vs LAMMPS kk **mixed** precision; `srun --mpi=pmi2`; 4-GPU
  `gpu/aware on` + `MPIR_CVAR_CH4_IPC_GPU_HANDLE_CACHE=generic`. VRAM = peak
  summed `nvidia-smi memory.used` over the GPUs used.
- LAMMPS throughput uses the **production** `Loop time` (the deck does an untimed
  warmup `run` first) — `lib_scaling.sh::lmp_throughput` parses the last one.
- **CPU: 1–1024 MPI ranks** (1 8 32 128 256 512 1024; 256 cores/node, so 512/1024
  span 2/4 nodes). **RSS = sum of per-rank peak VmHWM** for BOTH engines —
  AV self-reports it (`peak_rss_total_kb`, MPI-summed `/proc` VmHWM); LAMMPS is
  wrapped in `/usr/bin/time -v` (`lib_scaling.sh::rss_hwm_gb`). The old `pv_cpu`
  node-`free`-delta read LAMMPS as ~0 and inflated AV/LAMMPS to a bogus ~12×.

## Memory-optimization pass (2026-06-08)

Five per-atom / per-rank memory reductions plus a neighbor-build speedup, each
verified parity-clean (GPU `test_cuda_multi_gpu` 18/18 etc.; CPU MPI 18/18):

- **GPU** (`968d036`): fold `atom_bins[]` into `bin_pos.w` (one coalesced load
  carries position + index) and tighten the Verlet start `vl_max` 128→116 →
  **LJ VRAM −6…−9 %** with a **+2…+5 % throughput** bonus (fewer build loads);
  the same build-loop saving lifts **PE GPU throughput +1.6…+8 %** (the ng4
  reneighbor lever). Morton/cluster arrays were already gated off the MD path.
- **CPU** (`968d036`): lazy `atom_style`-molecular storage — the ~0.5 KB/atom
  wide bond/angle/dihedral/special arrays allocate only when topology exists, so
  pure-LJ pays nothing. Drops **LJ RSS substantially** and lets the 8 M-atom
  `strong @1–2 ranks` cases run (they used to OOM).
- **CPU** (`09f79a8`): per-step atom migration rewritten from `MPI_Alltoallv`
  (dense all-to-all) to **LAMMPS-style staged nearest-neighbor exchange**
  (x→y→z, multi-hop forwarding). O(1) point-to-point + O(log P) reduce fabric
  state instead of O(world). Correctness-preserving (identical end state; energy
  bit-parity). _Note: this did **not** move the cross-node RSS — at weak scaling
  migration ships only a few boundary atoms, so its buffers were never the RSS
  peak. The residual cross-node growth (below) is MPI-runtime per-rank state, not
  AV data._

## GPU results — AV vs LAMMPS-KOKKOS (GH200, run `968d036_2026-06-08`)

Solid = throughput (left axis); dashed = peak VRAM (right). Mean ± sample-std,
3 reps. AV leads at 1–3 GPU (kernel efficiency); LAMMPS-KOKKOS scales steeper
and edges ahead at 4 GPU (AV's multi-GPU comm/decomposition overhead).

![AV vs LAMMPS-KOKKOS GPU scaling](scaling_2x2_gpu_errorbar.png)


**LJ — Lennard-Jones (argon)**  (run `968d036_2026-06-08`, n=3 reps)

_weak scaling_

| GPUs | atoms | AV | LAMMPS-KOKKOS | AV/LAMMPS |
|---:|---:|---:|---:|---:|
| 1 | 8M | 2892 ± 5 | 1752 ± 145 | 1.65× |
| 2 | 16M | 4814 ± 3 | 2877 ± 31 | 1.67× |
| 3 | 24M | 6358 ± 12 | 3993 ± 71 | 1.59× |
| 4 | 32M | 5333 ± 10 | 5499 ± 287 | 0.97× |

_strong scaling_

| GPUs | atoms | AV | LAMMPS-KOKKOS | AV/LAMMPS |
|---:|---:|---:|---:|---:|
| 1 | 32M | 2972 ± 7 | 1779 ± 1 | 1.67× |
| 2 | 32M | 3757 ± 8 | 3325 ± 23 | 1.13× |
| 3 | 32M | 6455 ± 7 | 4167 ± 32 | 1.55× |
| 4 | 32M | 5327 ± 12 | 5830 ± 15 | 0.91× |

**PE — polyethylene melt (UA)**  (run `968d036_2026-06-08`, n=3 reps)

_weak scaling_

| GPUs | atoms | AV | LAMMPS-KOKKOS | AV/LAMMPS |
|---:|---:|---:|---:|---:|
| 1 | 4M | 426 ± 9 | 216 ± 5 | 1.97× |
| 2 | 8M | 429 ± 8 | 441 ± 15 | 0.97× |
| 3 | 12M | 562 ± 0 | 529 ± 47 | 1.06× |
| 4 | 16M | 703 ± 12 | 936 ± 29 | 0.75× |

_strong scaling_

| GPUs | atoms | AV | LAMMPS-KOKKOS | AV/LAMMPS |
|---:|---:|---:|---:|---:|
| 1 | 16M | 440 ± 1 | 231 ± 21 | 1.91× |
| 2 | 16M | 462 ± 2 | 472 ± 21 | 0.98× |
| 3 | 16M | 587 ± 2 | 685 ± 10 | 0.86× |
| 4 | 16M | 695 ± 7 | 937 ± 25 | 0.74× |


## CPU results — AV vs LAMMPS (1–1024 MPI ranks, run `09f79a8_2026-06-08`)

Fair per-rank-VmHWM-sum RSS for both engines (see Methodology). AV CPU
throughput leads LAMMPS across almost the whole range (LJ 1.0–2.6×, PE
1.1–1.4×). **Memory:** at 1 rank AV ≈ LAMMPS (73 vs 70 MB/rank); single-node
(≤256 ranks) AV is **1.0–1.8×** LAMMPS — _not_ the old "12×", which was the
`pv_cpu` measurement artifact. Cross-node, AV's **per-rank RSS grows with rank
count** (LJ-weak → ~3× LAMMPS at 1024) while LAMMPS stays flat (~150 MB/rank).
This was traced to MPI-runtime per-rank state (libfabric/cxi), **not** atom
migration (the staged rewrite above left it unchanged) and **not** any
O(world_size) allocation in AV's code — both engines link the same
MPICH 4.3.2 + libfabric. Pinpointing it needs runtime memory profiling; it is
within the node memory budget (487 MB/rank of ~2.9 GB available at 1024).

![AV vs LAMMPS CPU scaling](scaling_2x2_cpu_errorbar.png)


**LJ — Lennard-Jones (argon)**  (run `09f79a8_2026-06-08`, n=3 reps)

_weak scaling_

| MPI | atoms | AV | LAMMPS-KOKKOS | AV/LAMMPS |
|---:|---:|---:|---:|---:|
| 1 | 195K | 14 ± 0 | 5 ± 0 | 2.66× |
| 8 | 2M | 74 ± 1 | 42 ± 0 | 1.76× |
| 32 | 6M | 271 ± 0 | 154 ± 0 | 1.76× |
| 128 | 26M | 698 ± 3 | 496 ± 20 | 1.41× |
| 256 | 51M | 954 ± 2 | 865 ± 5 | 1.10× |
| 512 | 103M | 1731 ± 10 | 1719 ± 3 | 1.01× |
| 1024 | 204M | 2683 ± 6 | 2764 ± 2 | 0.97× |

_strong scaling_

| MPI | atoms | AV | LAMMPS-KOKKOS | AV/LAMMPS |
|---:|---:|---:|---:|---:|
| 1 | 8M | 14 ± 0 | 5 ± 0 | 2.55× |
| 8 | 8M | 96 ± 0 | 42 ± 0 | 2.30× |
| 32 | 8M | 331 ± 3 | 152 ± 1 | 2.18× |
| 128 | 8M | 1017 ± 0 | 517 ± 2 | 1.97× |
| 256 | 8M | 1854 ± 38 | 994 ± 21 | 1.87× |
| 512 | 8M | 3199 ± 136 | 2041 ± 53 | 1.57× |
| 1024 | 8M | 2873 ± 30 | 2639 ± 34 | 1.09× |

**PE — polyethylene melt (UA)**  (run `09f79a8_2026-06-08`, n=3 reps)

_weak scaling_

| MPI | atoms | AV | LAMMPS-KOKKOS | AV/LAMMPS |
|---:|---:|---:|---:|---:|
| 1 | 74K | 1 ± 0 | 1 ± 0 | 1.40× |
| 8 | 590K | 11 ± 0 | 8 ± 0 | 1.42× |
| 32 | 2M | 41 ± 0 | 30 ± 0 | 1.39× |
| 128 | 9M | 140 ± 0 | 114 ± 1 | 1.23× |
| 256 | 19M | 241 ± 0 | 217 ± 1 | 1.11× |
| 512 | 38M | 509 ± 2 | 416 ± 3 | 1.23× |
| 1024 | 75M | 842 ± 3 | 694 ± 1 | 1.21× |

_strong scaling_

| MPI | atoms | AV | LAMMPS-KOKKOS | AV/LAMMPS |
|---:|---:|---:|---:|---:|
| 1 | 988K | 1 ± 0 | 1 ± 0 | 1.40× |
| 8 | 988K | 11 ± 0 | 8 ± 0 | 1.44× |
| 32 | 988K | 40 ± 0 | 30 ± 0 | 1.35× |
| 128 | 988K | 137 ± 1 | 110 ± 1 | 1.25× |
| 256 | 988K | 237 ± 1 | 200 ± 1 | 1.19× |
| 512 | 988K | 395 ± 2 | 321 ± 3 | 1.23× |
| 1024 | 988K | 352 ± 15 | 491 ± 1 | 0.72× |


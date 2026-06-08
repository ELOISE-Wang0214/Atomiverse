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
```bash
sbatch benchmarks/md/bench_scaling.sbatch lj      # -> lj/runs/<sha>_<date>/
sbatch benchmarks/md/bench_scaling.sbatch pe      # -> pe/runs/<sha>_<date>/
python benchmarks/md/plot_scaling.py              # latest run per system -> scaling_2x2.png
python benchmarks/md/plot_scaling.py --id <sha>_<date>   # pin a specific run
```

## Why this is auditable + comparable (可审计 / 可环比)
- **One code set:** `bench_scaling.sbatch`/`plot_scaling.py`/`lib_scaling.sh` are
  committed once and reused; no per-date duplicates.
- **Builds committed HEAD** and writes the git SHA into `meta.json`; the run dir
  is `runs/<sha>_<date>/` → every number traces to exact code.
- **`config.sh` per system** is the only place sizes/params live → no drift.
- **`results.csv` schema is fixed:** `engine,system,mode,ng,atoms,matom_step_s,vram_gb`.
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

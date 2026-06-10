# MD benchmarks — AV vs LAMMPS

Every benchmark system has exactly two parts:

- **`<sys>/validation/`** — scientific correctness vs an external reference
  (LAMMPS). Plain executables on the `common/bench_check.h` PASS/FAIL
  harness; non-zero exit on any failing case. Run manually, not by CI.
- **`<sys>/performance/`** — speed / scaling vs LAMMPS. CSV-emitting
  executables + the system's LAMMPS deck and sizing config; results are
  committed per run under `performance/runs/<sha>_<date>/`.

```
benchmarks/md/
├── README.md
├── common/                      shared infra
│   ├── bench_check.h            PASS/FAIL harness for validation benches
│   ├── lib_scaling.sh           helpers (peak RSS/VRAM, LAMMPS log parse)
│   ├── plot_scaling.py          the ONE plotter (CPU + GPU)
│   ├── bench_scaling.sbatch     GPU scan driver (run per system)
│   └── bench_scaling_cpu.sbatch CPU scan driver, 1–1024 MPI ranks
├── lj/                          argon Lennard-Jones
│   ├── validation/              7 LAMMPS API-alignment / ensemble-parity benches
│   └── performance/             scaling.cpp (CPU MPI), bench_md_step.cu,
│                                bench_multi_gpu_scaling.cu, bench_torch.cu,
│                                bench_minimize.cu + in.lj_min_kk (CUDA
│                                minimizer CG/LBFGS vs LAMMPS minimize on
│                                identical coordinates via data-file handoff),
│                                in.lj_kk, config.sh, runs/
└── pe/                          polyethylene melt (Dreiding-UA)
    ├── validation/              bench_pe_bonded_lammps_xcheck (C20H42 vs
    │                            LAMMPS, 1e-6 eV) + reference/ provenance
    └── performance/             polyethylene_tension.cpp (MSU 5-stage
                                 protocol), in.pe_kk, pe_base.data,
                                 gen_pe_chain_mol.py, config.sh, runs/
```

Committed: `results*.csv` per run (an auditable history of numbers vs
code — the run dir is `<sha>_<date>`, every number traces to exact
code). Gitignored: `logs/`, `*.png`, `.scanout/` (figures are
regenerated, never committed).

## Validation (scientific correctness)

Build with `-DATV_BUILD_BENCHMARKS=ON`, then run the binaries directly;
each prints PASS/FAIL per case and exits non-zero on any failure:

```bash
./build-cpu/benchmarks/bench_setup_lammps_align        # lj: setup parity
./build-cpu/benchmarks/bench_lammps_npt_parity         # lj: NPT ensemble
./build-cpu/benchmarks/bench_lammps_langevin_parity    # lj: NVT ensemble
./build-cpu/benchmarks/bench_pe_bonded_lammps_xcheck   # pe: bonded E vs LAMMPS
# ... (deform / nve_limit / thermostat / pressure alignment)
```

The same physics is additionally gated in CI by the unit-test suite
(e.g. `test_pe_energy_xcheck`, `test_pe_tension_protocol`); these
standalone benches exist for manual third-party comparison without
gtest.

## Performance (speed / scaling)

```bash
# CPU, 1–1024 MPI ranks (256 cores/node; 512/1024 span 2/4 nodes):
cd <repo>     # submit from inside the repo (worktrees work too)
LMP=/path/to/lmp sbatch benchmarks/md/common/bench_scaling_cpu.sbatch          # LJ
SYS=pe LMP=/path/to/lmp sbatch benchmarks/md/common/bench_scaling_cpu.sbatch   # PE
# overrides: CORES="1 8 32 ..." STEPS= WEAK_PER_RANK= STRONG_TOTAL=

# GPU (GH200, 1–4 GPUs):
sbatch benchmarks/md/common/bench_scaling.sbatch lj
sbatch benchmarks/md/common/bench_scaling.sbatch pe

# plot (the sbatch scripts auto-plot; to re-plot on the login node):
python3 benchmarks/md/common/plot_scaling.py --backend cpu --out scaling_cpu.png
python3 benchmarks/md/common/plot_scaling.py --backend cpu --id <sha>_<date> --md table.md
```

### Plot convention (every scaling figure — do not change)

- **left y-axis** = throughput [Matom·step/s] (solid, value-labelled), AV vs LAMMPS.
- **right y-axis** = peak memory (dashed) — RSS [GB] for CPU, VRAM [GB] for GPU.
- **bottom x-axis** = device count (MPI ranks / GPUs); **top x-axis** = total atoms.

### Methodology (held identical across runs)

- **LJ:** AV `scaling` (sc ρ\*=0.712, rc=2σ) vs LAMMPS `lj/performance/in.lj_kk`
  (`lattice sc 0.712`, `lj/cut 2.0`). CPU weak = 200k atoms/rank, strong = 8M total.
- **PE:** AV `polyethylene_tension --pe-melt-grid` (lj/cut 10.5 + harmonic
  bond/angle + multi/harmonic dihedral, special_bonds 0 0 0) vs LAMMPS
  `pe/performance/in.pe_kk` (**same grid as AV** → identical structures),
  MSU Tschopp 2010 Dreiding-UA FF.
- **CPU memory:** RSS = sum of per-rank peak VmHWM for BOTH engines — AV
  self-reports it (`peak_rss_total_kb`, MPI-summed /proc VmHWM); LAMMPS is
  wrapped in `/usr/bin/time -v` (`lib_scaling.sh::rss_hwm_gb`).
- **GPU memory:** VRAM = peak summed `nvidia-smi memory.used` over the GPUs used.
- LAMMPS throughput uses the **production** `Loop time` (untimed warmup first).

## Latest CPU results — run `2bca215_2026-06-09` (1–1024 ranks)

After the 2bca215 memory/scaling pass (subdomain-local cell grid, chunked
capacity, conditional tag maps, bonded tuple cache):

**LJ throughput [Matom·step/s]** — AV leads at every point:

| ranks | weak AV | weak LMP | | strong AV | strong LMP |
|---:|---:|---:|---|---:|---:|
| 8    | 78   | 42   | | 96   | 42   |
| 128  | 775  | 479  | | 1062 | 517  |
| 512  | 2285 | 1717 | | 4360 | 2082 |
| 1024 | 4059 | 2762 | | 4678 | 2691 |

**LJ memory (sum VmHWM)** — per-rank RSS is now flat (~192 MB/rank weak;
was growing 187→487 MB/rank at 8→1024 ranks before 2bca215):

| | AV | LAMMPS | ratio |
|---|---:|---:|---:|
| weak 1024 (204M atoms) | 196 GB | 158 GB | 1.24× (was 3.15×) |
| strong 1 (8M atoms)    | 1.74 GB | 2.07 GB | **0.84×** |

**PE:** weak throughput leads LAMMPS 1.2–1.5× through 1024 ranks
(913 vs 693 @1024); strong-512 memory 55.3 vs 64.6 GB (below LAMMPS).
Known gap: PE strong-1024 (≈960 atoms/rank, latency-bound) trails
LAMMPS (350 vs 471).

## Latest GPU results — run `968d036_2026-06-08` (GH200, 1–4 GPUs)

AV leads at 1–3 GPU (kernel efficiency); LAMMPS-KOKKOS scales steeper
and edges ahead at 4 GPU (AV multi-GPU comm/decomposition overhead).
LJ weak AV/LAMMPS: 1.65× / 1.67× / 1.59× / 0.97× at 1/2/3/4 GPU;
PE weak: 1.97× / 0.97× / 1.06× / 0.75×.

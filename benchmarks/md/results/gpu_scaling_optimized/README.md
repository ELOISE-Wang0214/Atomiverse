# CUDA LJ GPU scaling — optimized Atomiverse vs LAMMPS KOKKOS (2026-06-06)

Single- and multi-GPU Lennard-Jones throughput of the **optimized** Atomiverse
CUDA backend against **LAMMPS KOKKOS/CUDA**, on one NVIDIA GH200 node, on one
shared LJ system so the numbers are directly comparable.

Single-card (throughput + VRAM vs atom count) and multi-card (weak scaling) are
separate figures; each overlays throughput (solid, left axis) and peak VRAM
(dashed, right axis) for both codes.

![Single GH200](gpu_single.png)
![Multi GH200 weak scaling](gpu_multi.png)

## The system

Simple-cubic LJ fluid, identical physics in both codes:

| quantity | value |
|----------|-------|
| reduced density ρ\* | 0.712 (lattice spacing 1.12 σ) |
| cutoff | 2 σ |
| Verlet skin | 0.6 σ (= AV's 0.3·cutoff auto-skin) |
| reneighbor | every 20 steps (`neigh_modify every 20`) |
| precision | mixed (FP32 geometry/force, FP64 energy/virial) |
| timed steps | 100 (single), 500 (multi) after warmup |

AV is built by `bench_multi_gpu_scaling.cu` (cluster lattice + jitter); LAMMPS by
`benchmarks/md/reference/lj_scaling/in.lj_kk` (`units lj`, same ρ\*, rc, skin).

## The optimization (this is the headline)

The baseline AV CUDA LJ force ran at ~1000 Matom-step/s — **~2× behind** Kokkos.
An `nsys` + `ncu` profile showed the LJ force kernel was **compute-bound with
only ~7 of 32 warp lanes active**: the cluster-CSR force walk used 8 threads per
cluster, each redundantly iterating the cluster's *union* neighbor list with
divergent control flow. Memory was nearly idle (~19 %), so the README's earlier
"uncoalesced gather" diagnosis was wrong.

**Fix:** replace the cluster-CSR force traversal with a **tight per-atom Verlet
list**, built once per rebuild from the cluster CSR (`build_verlet_list`) and
iterated by a **one-thread-per-atom** force kernel (`lj_force_verlet`). This
removes the 8×-per-cluster redundancy and the divergence, and the column
capacity is **auto-sized from density** so it also works for the much denser
PE-melt shell. All CUDA correctness tests pass, including the nlist-vs-brute-force
force-parity check.

Result: **single-GPU LJ went from ~0.5× to ~1.15× LAMMPS Kokkos.**

## Single-GPU result (panel a)

| atoms | AV Matom-step/s | Kokkos | AV / Kokkos |
|------:|----------------:|-------:|:-----------:|
|   1 M | 1736 | 1251 | **1.39×** |
|   2 M | 1998 | 1662 | **1.20×** |
|   4 M | 1997 | 1717 | **1.16×** |
|   8 M | 2097 | 1775 | **1.18×** |
|  16 M | 2039 | 1791 | **1.14×** |
|  32 M | 2137 | 1807 | **1.18×** |

Atomiverse **beats LAMMPS Kokkos at every size ≥ 1 M** (≈ +14–39 %). Below
~262 k atoms Kokkos still wins (AV kernel-launch/sync overhead dominates a tiny
GPU).

**VRAM tradeoff (honest):** AV uses more memory — the per-atom Verlet list (a
density-sized column buffer) plus some still-allocated cluster/flat buffers. At
32 M: AV ≈ 41 GB vs Kokkos ≈ 14 GB. The list footprint was already cut ~40 %
(tight single-GPU stride + a 128-floor cap); the remaining gap (still-allocated
mask/flat buffers that the Verlet force path no longer reads) is a clear
follow-up. Throughput is unaffected by the cap.

## Multi-GPU result (panel b) — WEAK scaling

This panel is **weak scaling**: the per-card load is fixed at the single-card
Pareto size (≈ 48 M atoms/card, near the VRAM-limited max where single-card
throughput is highest) and N cards run N × that. Ideal = total throughput grows
linearly with GPU count.

First, two real multi-GPU **bugs** were found and fixed (a delta measurement
excluding one-time setup proved it wasn't a timing artifact):
1. `MultiGpu::migrate()` ran its full classify + gather + cross-GPU P2P + publish
   **every reneighbor even when no atom changed subdomain** → fixed with an
   early-out (skip all data movement in the common stable case).
2. A density-based Verlet-cap heuristic under-sized subdomains
   (`n_subdomain / V_global`) → 4-GPU overflow crash → fixed with a fixed-128
   floor + **self-healing regrow**. All CUDA correctness tests pass.

A profile then showed the real per-step bottleneck was **host-side**: ~1000
`cudaMalloc`/`cudaFree` calls/run (~600 ms of a 905 ms 30-step 4-GPU wall) from
per-reneighbor halo/migrate scratch — **not** the P2P copies (12 ms total) and
**not** the NVLink (`nvidia-smi topo`: full mesh, every pair `NV6` = 6 NVLinks,
~300 GB/s, P2P enabled). **Fix:** stream-ordered **pooled** allocation
(`cudaMallocAsync`/`cudaFreeAsync` + a retain-threshold mem-pool, drained at
teardown) so per-reneighbor allocations reuse the pool; + migrate's counting
sort is skipped on the no-migration path. Suite re-run 5× — all green.

Weak scaling (500 steps, 48 M/card) after these fixes:

| GPUs | total Matom-step/s | per-card | parallel eff. |
|----:|-------------------:|---------:|:-------------:|
| 1 | 2343 | 2343 | 100% |
| 2 | **2909** | 1455 | **62%** |
| 4 | 1759 | 440 | 19% |

So eliminating the malloc storm lifted **2-GPU weak scaling to 1.24×** (was
anti-scaling). 4-GPU weak at this large per-card is still limited: the malloc
storm is gone, but the halo is still **synchronous** (blocking `cudaMemcpyPeer`
+ a full cross-device barrier every step) and the 2×2×1 decomposition has more,
larger face exchanges on the critical path.

Two corroborating points:
- At **smaller** per-card (strong scaling, 64 M total → 16 M/card at 4 GPU) the
  engine scales **positively** (2373 → 2964 Matom-step/s, 1.25×) — the halo is a
  smaller fraction of the per-card work.
- LAMMPS KOKKOS (1 rank/GPU, `gpu/aware`) also does not weak-scale here
  (per-card 1634 → 376), so both codes are halo-bound at this atoms/GPU ratio.

**Remaining lever (software):** an async stream/event halo with interior/boundary
comm-compute overlap. The NVLink mesh makes near-linear scaling achievable — this
is the documented next step for 4-GPU weak scaling.

## Files

| file | contents |
|------|----------|
| `av_single.csv` / `lmp_single.csv` | single-GPU sweep (32k–32M / 64M) |
| `av_multi.csv` / `lmp_multi.csv` | multi-GPU strong scaling (16/32/64 M × 1/2/4 GPU, 500 steps) |
| `plot_gpu.py` | the figure script (reads the CSVs co-located here) |
| `in.lj_kk` (under `reference/lj_scaling/`) | the LAMMPS deck |

## Reproduce

```bash
# on a GH200 node, from the worktree
sbatch -p gpu --gpus 4 -t 01:30:00 -c 32 --mem=180G \
  --export=ALL,WORK_DIR=$PWD,RESULTS_DIR=$PWD/benchmarks/md/results/gpu_scaling_optimized \
  tools/gpu_final_all.sh        # build + correctness + single + multi sweeps
python benchmarks/md/results/gpu_scaling_optimized/plot_gpu.py
```

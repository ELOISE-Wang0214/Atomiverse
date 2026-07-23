# CUDA OpenMM MD — GH200 timing note

**Date:** 2026-07-23  
**Chain:** `cuda-openmm-md` (01–05)

## What changed

CUDA classical MD neighbor + nonbonded force is an **OpenMM cluster/tile engine**:

- `kClusterSize = 32`
- FULL cluster CSR + multi-word interaction masks (no `memset 0xFF`)
- Tile force kernels with per-warp shared-memory neighbor loads (`lj_force_tile`, `coul_long_tile`)
- multi-GPU interior classification from cluster CSR (no `vl_*`)
- Deleted `nlist_direct.cu` / `launch_build_verlet_direct` default path

## Numbers

**Login node cannot build CUDA** — fill in after GH200 `ctest` + optional
`tools/cuda_nonbonded_style_bench` (now cluster-only):

| Case | N | steps | rebuild cadence | ms/step (or step/s) | notes |
|------|---|-------|-----------------|---------------------|-------|
| TBD  |   |       |                 |                     | GH200 |

Historical context: commit `0431eb13` claimed ~16× by switching **from** a
half-finished CS=8 cluster shell **to** Verlet — that measured an incomplete
OpenMM path, not engine architecture. Re-measure after this rewrite before
any speed claims.

## Honesty rules

- Do **not** claim “faster than OpenMM official” without a real binary A/B.
- Do **not** cite the old style A/B bench as architecture truth.

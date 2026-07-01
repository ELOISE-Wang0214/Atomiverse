# pw/si/convergence — plane-wave energy vs ecut & k-mesh

Convergence study for the manual "plane-wave DFT for solids" chapter
(manual-system §3 chapter 8, which calls for a convergence study alongside the
Si EOS validation). A third study type alongside the case's `validation/`
(EOS a0/B0 vs CP2K) and `performance/` (SCF time-to-solution) parts.

## What it checks

`bench_pw_convergence` runs two independent 1-D sweeps on a **fixed 2-atom fcc
Si cell** (a = 5.43 Å, the ONCVPSP LDA UPF shipped at
`../validation/reference/Si_lda_nc.upf`) and asserts the
**physical-approximation** result class (manual-system §5) — monotone
convergence toward the converged limit, *not* a single reference energy
(absolute PW total energies are not comparable across codes/pseudopotentials;
that is exactly why the `validation/` bench compares the pseudo-insensitive
a0/B0 instead):

1. **ecut sweep** at fixed k-mesh — **variational**. Enlarging the plane-wave
   basis (higher `ecut`) can only lower the total energy; the per-step lowering
   shrinks and plateaus. The gate: `E(ecut_{k+1}) ≤ E(ecut_k) + slack`
   (`--slack`, default 1e-4 Ry) and the deltas shrink.
2. **k-mesh sweep** at fixed ecut — Brillouin-zone sampling. The energy is
   **not** variational in the k-mesh, so the gate is only that the successive
   change *shrinks* (`|dE|` monotone-decreasing) — the signature of a
   converging k-integration.

Exit code 0 = both sweeps converge; 1 = a point raised the energy past the
variational slack, a k-mesh step failed to shrink, or an SCF diverged.

## Expected shape (orientation, not a gate)

- **E vs ecut**: the ONCV Si pseudo needs ~60 Ry for sub-mRy total-energy
  convergence (see `../README.md`); below that the energy is still dropping.
  The curve is monotone-decreasing and flattening — that flattening is what the
  benchmark asserts, at whatever the current solver produces.
- **E vs k-mesh**: 2×2×2 leaves several-mRy BZ-sampling error; 4×4×4 and finer
  converge it. `|dE|` between successive meshes shrinks toward zero.

## Running

```bash
cmake -B build -DATV_BUILD_PW=ON -DATV_BUILD_RHF=ON -DATV_BUILD_BENCHMARKS=ON ...
cmake --build build --target bench_pw_convergence

# CPU (login node) — SLOW (single-threaded PW SCF): keep the lists small.
./build/benchmarks/bench_pw_convergence --ecut-list 15,20,25 --kmesh-list 2,3

# GPU node — the Davidson hot path runs on device; production-scale sweep.
./build/benchmarks/bench_pw_convergence --backend cuda \
    --ecut-list 20,30,40,60 --kmesh-fixed 4 \
    --kmesh-list 2,4,6 --ecut-fixed 40
```

CPU plane-wave SCF is single-threaded and slow (minutes per point). The
defaults are a small, fast **smoke** sweep; production convergence uses the
larger `--backend cuda` lists above. The `cuda::pw` and `cpu::pw` backends
produce bit-identical energies (`pw/README.md`), so the shape is backend-independent.

The UPF path defaults to the committed `validation/reference/Si_lda_nc.upf`
(via `ATV_PW_DATA_DIR`, set by `benchmarks/CMakeLists.txt`); override with
`--upf PATH`.

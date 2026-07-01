# md/lj/consistency — cpu↔cuda backend-consistency parity

Realizes the **manual-system §6 backend-consistency policy** as a committed
benchmark case: *CPU is the reference; CUDA is compared, not assumed identical.*
This is the "cpu↔cuda parity harness" listed as a required, high-priority
industrial deliverable (manual-system §10, discovery §8 gap #2).

It is a **third part** of the `md/lj` case, alongside `validation/` (science vs
LAMMPS) and `performance/` (speed/scaling). Backend parity compares two
Atomiverse backends against each other, so — like every cross-engine check —
it lives in `benchmarks/`, never in the `ctest` gate.

```
consistency/
├── consistency_system.h   shared FCC-argon system + golden-file I/O
├── cpu_reference.cpp       CPU (reference) → writes energy+forces golden
├── cuda_compare.cu         GPU → replays, reads golden, asserts tolerances
└── README.md               this file
```

## The system

A 3×3×3 FCC argon crystal (a = 5.26 Å → **108 atoms**) with a deterministic
sinusoidal jitter (0.05 Å) that breaks the lattice symmetry — a perfect FCC
lattice has a vanishing net per-atom force by symmetry and is useless as a
parity target. LJ parameters are the canonical argon set shared with the
`md/lj` validation benches (ε = 0.0104 eV, σ = 3.4 Å, rc = 6.8 Å). Both
backends build this system from the single `make_consistency_system()` in the
shared header, so they evaluate the *identical* configuration.

## Why this is NOT bitwise identical

The CUDA backend stores per-atom state at `real_t` (**float** —
`cuda/frame/types.cuh`) and does the pair math at `work_t`; the CPU reference
is all-`double`. Float-stored positions carry ~1e-7 relative error, which the
`r^-12` LJ core amplifies. Energy accumulates in `AccReal = double` on both
backends, but the CUDA energy sums float-stored contributions, so it too
carries the float floor. This is the **`backend-dependent` result class**
(manual-system §5): a per-quantity tolerance to the storage-precision floor,
documented, never a bit-for-bit or `EXPECT_EQ` comparison.

## Tolerance policy (per quantity)

| quantity | rule | default | rationale |
|---|---|---|---|
| energy | `\|ΔE\| ≤ etol_rel · \|E_cpu\|` | `etol_rel = 2e-4` | float-stored positions → LJ energy; double accumulator over 108 atoms |
| force (per component) | `\|ΔFc\| ≤ ftol_rel · \|Fc_cpu\| + ftol_abs` | `ftol_rel = 2e-3`, `ftol_abs = 1e-3 · F_rms` | relative floor for the r^-12 amplification; the absolute term bounds near-zero components against the force scale, not against 0 |

The defaults are the **initial** float-storage floor. On the first real GPU
run, read off the reported `max rel` / `rms |dF|` and **calibrate**: set the
tolerances just above the observed mixed-precision residue. Thereafter only
ever **tighten** them — loosening a parity tolerance to make a red run green
hides exactly the regression this case exists to catch.

Comparing at `--steps 0` (post-setup single force field) is the tightest,
least-chaotic point. `--steps N > 0` lets mixed-precision divergence accumulate
along the trajectory (LJ dynamics is chaotic); use it only to *probe* drift and
loosen tolerances accordingly — it is not the parity gate.

## Running (GPU node only)

The comparison half is `ATV_BUILD_CUDA`-gated and runs only on the GH200 node;
the login nodes have no CUDA toolkit. Build both halves with
`-DATV_BUILD_CUDA=ON -DATV_BUILD_BENCHMARKS=ON`, then:

```bash
# 1. CPU reference writes the golden (all-double reference values).
./build-cuda/benchmarks/bench_lj_consistency_cpu  --golden lj_golden.txt

# 2. GPU replay reads the golden and asserts the mixed-precision tolerance.
./build-cuda/benchmarks/bench_lj_consistency_cuda --golden lj_golden.txt
#    exit 0 = within tolerance, 1 = FAIL (prints the worst |dE| / |dF|)
```

Both halves must use the **same** `--steps` / `--dt` so the phase points match.

### Golden file is regenerated, not committed

The golden is produced by our *own* reference backend (CPU, deterministic), so
it is regenerated on every run rather than checked in. A committed golden would
either drift silently from the code or, worse, tempt hand-editing — and a
fabricated parity reference defeats the check. (This differs from `validation/`
reference data, which comes from an *external* engine and therefore is
committed with provenance in `SOURCES.md`.)

## CI wiring

- `bench_lj_consistency_cpu` — CPU-only, builds anywhere; can smoke-run in the
  CPU matrix to prove the golden writer stays green.
- `bench_lj_consistency_cuda` — **nightly / self-hosted GPU job only.** Add the
  two-step run above to `.github/workflows/benchmark.yml` (the self-hosted GPU
  workflow), never to `ci.yml` (CPU-only, no CUDA toolkit).

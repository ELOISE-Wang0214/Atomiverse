# Performance-regression gate — design

**Status:** design / not yet wired. Addresses discovery §8 gap #3 — *"No
automated perf-regression gate (perf CSVs are committed history)."* Today
`*/performance/runs/<sha>_<date>/results_*.csv` is an auditable **archive** of
numbers vs code, but nothing compares a new run against the archive: a
throughput regression only shows up if a human eyeballs `plot_scaling.py`. This
document specifies the gate that closes that loop.

## 1. What we already commit (the raw material)

Every performance run is committed under its case's `performance/runs/`:

```
benchmarks/md/lj/performance/runs/<sha>_<date>/
  results_cpu.csv        engine,system,mode,ng,atoms,matom_step_s,vram_gb
  results_gpu.csv        (same schema; vram_gb = VRAM on GPU, RSS on CPU)
  meta_cpu.json / meta.json
      { system, av_git_sha, date, host, backend, cores, steps, lammps }
benchmarks/pw/si/performance/            (CSV emitted by bench_pw_scf_perf)
  backend,n_atoms,a_angstrom,ecut_ry,k_mesh,threads,iterations,converged,
  energy_ry,wall_s,s_per_iter
```

google-benchmark targets (`bench_md_step`, `bench_native`) instead emit JSON
(`--benchmark_format=json`) with `real_time` / `items_per_second` / counters.

Two metric polarities, handled oppositely by the gate:

| metric | source | better = | regression = |
|---|---|---|---|
| `matom_step_s` (throughput) | md `results_*.csv`, `Matom-step/s` counter | **higher** | current < baseline |
| `s_per_iter`, `wall_s` (latency) | pw perf CSV | **lower** | current > baseline |
| `real_time`, `cpu_time` | google-benchmark JSON | **lower** | current > baseline |
| `items_per_second` | google-benchmark JSON | **higher** | current < baseline |
| `vram_gb` (memory) | md `results_*.csv` | **lower** | current > baseline (looser band) |

## 2. Comparison key — never compare across hardware

A perf number is only meaningful against another number **on the same machine
with the same physical work**. The join key is:

```
(host, backend, engine, system, mode, ng)          # md CSVs
(backend, n_atoms, ecut_ry, k_mesh, threads)        # pw CSV
(benchmark_name, arg)                               # google-benchmark JSON
```

- `host` + `backend` come from `meta.json`; a run on `n215` never gates a run
  on the GH200. The gate **skips** (does not fail) any key without a
  same-host baseline.
- Only `engine == AV` rows are gated. LAMMPS/CP2K rows are an external
  reference plotted alongside — their speed is not ours to regress.

## 3. Baseline selection

Two supported modes; **pinned** is the default (explicit, auditable):

### Pinned baseline (recommended)
A committed `benchmarks/<domain>/performance/baseline.json` names the accepted
baseline run id per `(host, backend)`:

```json
{
  "n215/cpu":       { "run": "2bca215_2026-06-09", "blessed": "2026-06-09",
                      "note": "post memory/scaling pass" },
  "arrhenius/cuda": { "run": "968d036_2026-06-08", "blessed": "2026-06-08" }
}
```

The gate loads the baseline run's CSV, joins on the key, and compares. A
baseline is only moved by the **bless procedure** (§6) — so a regression can
never be laundered into the baseline silently.

### Rolling baseline (optional)
Baseline = per-key median of the last `K` committed runs on the same host
(default `K = 5`), excluding the current one. Cheaper to maintain (no
`baseline.json`), but a slow drift can walk the baseline downward unnoticed;
use only for noisy secondary metrics (e.g. `vram_gb`).

## 4. Threshold band + noise model

Perf is noisy; the band must exceed run-to-run variance so the gate fires on
real regressions, not jitter.

- **Per-run noise reduction (measurement side).** Each config is already a
  median-of-repeats with an untimed warmup (google-benchmark repetitions; the
  scaling drivers report the production `Loop time` after a warmup — see
  `md/README.md`). The gate consumes those already-reduced medians; it does not
  re-time.
- **Two-band verdict per key** (`r` = current/baseline, oriented so **>1 is
  always worse** — invert throughput ratios):

  | band | default | action |
  |---|---|---|
  | `r ≤ warn` | `warn = 1.03` (3%) | PASS |
  | `warn < r ≤ fail` | `fail = 1.08` (8%) | WARN (report, exit 0) |
  | `r > fail` | | FAIL (exit non-zero) |

  Memory (`vram_gb`) uses a looser `fail = 1.15` (allocator/page noise).

- **Aggregate guard against single-point noise.** A sweep has many `ng` points;
  one noisy point should not fail the job. FAIL the run only when **either**:
  (a) the **geometric-mean** ratio across the sweep exceeds `fail`, **or**
  (b) **≥ 2 keys** individually exceed `fail`. A lone point over `fail` is a
  WARN. (Tune the "≥ 2" to sweep length.)

- **Calibrate the bands from history.** Compute the per-key ratio between
  consecutive same-host runs already committed (e.g. `e375593` → `2bca215` on
  `n215/cpu`) to estimate real variance, and set `warn`/`fail` a comfortable
  margin above it. Bands are config, not code.

## 5. The gate (algorithm)

Proposed tool: `benchmarks/md/common/regression_gate.py` (sibling to
`plot_scaling.py`; reuse its CSV loader). Domain-agnostic — pass the CSV schema
and metric.

```
regression_gate.py --current runs/<sha>_<date>/results_cpu.csv
                   --baseline auto|<run_id>|baseline.json
                   --metric matom_step_s --higher-is-better
                   --warn 0.03 --fail 0.08 [--memory-metric vram_gb]

1. load current CSV + meta.json (host, backend)
2. resolve baseline CSV for (host, backend)          # pin or rolling; skip if none
3. filter engine == AV on both
4. inner-join on the key (§2); report keys present in one side only
5. per key: ratio = worse-is-larger(current, baseline)
6. classify each key PASS / WARN / FAIL against warn/fail bands
7. aggregate: geomean ratio + count of FAIL keys → run verdict (§4)
8. print a table (key, baseline, current, ratio, verdict) + PASS/WARN/FAIL
9. exit non-zero iff run verdict == FAIL
```

For the google-benchmark JSON targets (`bench_md_step`, `bench_native`), do
**not** re-implement this — use upstream `benchmark/tools/compare.py`:

```
compare.py benchmarks <baseline.json> <current.json> \
    --benchmark_repetitions=5     # U-test on medians
```

Commit the accepted baseline JSONs under
`benchmarks/<domain>/performance/perf_baselines/<host>_<backend>.json` and gate
on `real_time` / `items_per_second` with the same band, blessed via §6.

## 6. Bless procedure (moving a baseline)

A legitimate perf change (optimization, or an accepted cost from a correctness
fix) must update the baseline **deliberately and reviewably**:

1. Land the code change; run the perf suite → new `runs/<sha>_<date>/`.
2. Confirm the delta is intended (the commit that caused it is the perf commit).
3. Edit `baseline.json` / `perf_baselines/*.json` to the new run id, with a
   `note` and `blessed` date, in the **same PR**.
4. Reviewer sees both the code change and the baseline move together — the
   audit trail is the git history of `baseline.json`.

Never move a baseline in a PR that did not intend the perf change.

## 7. Where it runs

- **Self-hosted / nightly only.** Perf gates belong on dedicated hardware, not
  shared CI runners (the ci.yml CPU matrix is correctness-only and its runners
  are too noisy for a sub-10% gate). Wire `regression_gate.py` into the
  self-hosted `benchmark.yml` job after the run, and into the CPU scaling
  sbatch drivers (`md/common/bench_scaling*.sbatch`) after they emit the CSV.
- **On push to `master`** for the GPU job (matches `benchmark.yml`'s trigger),
  plus a manual/nightly CPU scaling run.
- **Non-blocking first.** Ship the gate in WARN-only mode (never exit non-zero)
  for a few runs to validate the bands against real variance, then flip the
  FAIL band on.

## 8. Scope / non-goals

- Not a per-PR CPU-CI check — too noisy, wrong hardware.
- Does not re-time or re-run benchmarks; consumes the committed medians.
- Does not gate external-reference (LAMMPS/CP2K) rows.
- Absolute-number drift (energies) is a *validation* concern (result classes,
  manual-system §5), not this gate — this gate is wall-clock/throughput/memory
  only.

# Atomiverse 1-rank PE-tension run (LAMMPS-comparable)

Captured for `parallel-md-finish-04-scaling-report`'s Table 2
single-core row.

## Command

```bash
OMP_NUM_THREADS=1 /usr/bin/time -l \
  ./build/benchmarks/polyethylene_tension \
    --mode full --box-len 90 \
    --out-dir benchmarks/md/reference/polyethylene/runs/atomiverse_n1 \
    --quiet
```

`OMP_NUM_THREADS=1` chosen for apples-to-apples comparison with `lmp_serial`
(LAMMPS uses 1 thread per rank by default). `--box-len 90` matches the
self-contained LAMMPS reference in `lammps.in`.

## Results

| Item | Value |
|---|---|
| Exit code | 0 |
| Wall time (under `time -l`) | 3.64 s |
| Wall time (under `time -p`, less overhead) | 2.85 s — number quoted in `scaling_report.md` Table 2 |
| Peak RSS | 21,708,800 B (20.71 MiB) |
| LAMMPS n=1 comparison wall (`lmp_serial`) | 2.46 s |
| Atomiverse / LAMMPS ratio | **1.16×** (≤ 1.20× target ✓) |

## Output artifacts (same schema as `runs/atomiverse_n8/`)

| File | Lines |
|---|---|
| `thermo.csv` | 281 |
| `stress_strain.csv` | 101 |
| `time.log` | 25 (`/usr/bin/time -l` raw output) |

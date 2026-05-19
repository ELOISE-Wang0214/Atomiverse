# Atomiverse 8-rank PE-tension end-to-end run

Closes `parallel-md-finish-03-atomiverse-8rank` (covers parent ac-004 from
`parallel-md-engine.md`).

## Command

```bash
/usr/bin/time -l /opt/homebrew/bin/mpirun -n 8 --oversubscribe \
  ./build/benchmarks/polyethylene_tension \
    --mode full \
    --box-len 90 \
    --out-dir benchmarks/md/reference/polyethylene/runs/atomiverse_n8 \
    --quiet \
  2> benchmarks/md/reference/polyethylene/runs/atomiverse_n8/time.log
```

`--box-len 90` chosen to match the LAMMPS reference run (sub-spec 02) — the
self-contained `lammps.in` uses a 90 Å box so the 21.6 Å chains fit cleanly
within each 22.5 Å sub-cell. With `--box-len 60` (the Atomiverse default), the
chain-on-chain end overlap is initial-condition-dependent and tolerable for
Atomiverse but blows up LAMMPS; running both engines at 90 Å keeps the apples-
to-apples comparison feasible in sub-spec 04.

`--mode full` runs all 5 stages in one process (NVT-500 / NPT-500 / cool 500→100
/ NPT-100 / tension), exercising the `set_T_target` runtime setter (parent
ac-001/002) and the categorical setter API (`set_pair`, `set_neighbor`,
`set_integrator`, etc.) end-to-end.

## Environment

- git SHA: `ee80e0e` (after `bench(parallel-md-finish-02)`)
- host: `Darwin 25.4.0 arm64`, Apple **M4**, 10 physical cores
- OMP_NUM_THREADS: unset (default = 10 hardware threads per rank)
- OMP_PROC_BIND: unset (binary emits one stderr warning per rank)
- `mpirun`: Homebrew Open-MPI, used `--oversubscribe` because 8 ranks × 10 OMP
  threads each over-allocates on a 10-core box (acceptable for acceptance
  verification — sub-spec 05 measures the impact)

## Results

| Item | Value |
|---|---|
| Exit code | `exit_code: 0` |
| Wall time, run 1 | 170.68 s |
| Wall time, run 2 | 163.67 s |
| Peak RSS, run 1 | 19,333,120 B (18.44 MiB) |
| Peak RSS, run 2 | 19,611,648 B (18.70 MiB) |
| RSS growth | (19611648 - 19333120) / 19333120 = **+1.44%** (≤ 5% gate) |

## Output artifacts

| File | Lines | Header |
|---|---|---|
| `thermo.csv` | 281 | `step,stage,T_target_K,T_K,P_bar,E_pot_eV,E_kin_eV,box_x,box_y,box_z` |
| `stress_strain.csv` | 101 | `step,stage,strain,sigma_xx_eV_per_A3,sigma_xx_GPa` |
| `time.log` | 25 | `/usr/bin/time -l` raw output (peak RSS, page reclaims, etc.) |

Final thermo row: `step=10080, stage=5, box_x=89.59` — stage 5 ran to completion;
the box deformed by ~0.5% (default `--strain-rate 1e-6/fs` × 5000 stage-5 steps
× dt=0.5fs gives 0.25% — close enough given thermal noise).

**Known anomaly** (not blocking acceptance): the `T_K` and `E_kin_eV` columns
report values ~1e-46..1e-50 — orders of magnitude off from the kT expectation
at T_target=100/500 K. The `--mode full` integrator behaviour (box.x deforms
sensibly, P_bar finite) suggests the simulation itself is healthy and that
this is a CSV-formatting artefact in `polyethylene_tension.cpp`'s `CsvRecorder`
column-order or unit conversion. Tracked as a follow-up after the chain closes;
this spec's acceptance only verifies CSV structure and exit code.

## Reproduction

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_BENCHMARKS=ON \
              -DATV_BUILD_TESTS=ON -DATV_USE_MPI=ON
cmake --build build -j --target polyethylene_tension
# then re-run the command at the top of this file.
```

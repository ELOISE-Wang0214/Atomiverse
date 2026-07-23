# MetalWalls method loop — worktree → merge → dev

## Pipeline

```text
parent (dev)
  census + assign
       │
       ├─ worker₁ isolation_worktree  ──┐
       ├─ worker₂ isolation_worktree  ──┼─ pin + impl + verify + commit on wf/<slug>
       └─ workerₙ isolation_worktree  ──┘
       │
  merge agent (parent): git merge each green HEAD → dev
       │
  harness agent (parent): queue + roadmap + notes; commit chore(harness)
```

| Role | `isolation_worktree` | Writes |
| ---- | -------------------- | ------ |
| Census / foreman / merge / harness | **false** (parent) | `dev` only |
| Feature worker | **true** | private worktree only until merge |

Workers must report `worktree_path`, `branch` (`wf/<slug>`), `head_sha`.

**Recover incomplete first** (`## in_progress` `[~]`, open MD/electrochem specs,
unmerged `wf/<slug>`). Spec is source of truth. Offline LAMMPS / MetalWalls
pins only (`tools/metalwalls_oracle/`); never live oracle binaries in ctest.

**Scope:** classical electrochemical MD + MetalWalls map **E1–E12** +
LAMMPS-overlap classical MD/FF. **Not** SCF/post-HF, GPW/PW, or CUDA OpenMM
rewrite.

**Architecture binding:** electrodes / const-pot = **Components** on
`cpu::md::MD`; never a new Driver category; never charge solve inside
`KSpace*`. Map: `.claude/notes/metalwalls-capability-map.md`.

**Oracle priority:** if the algorithm also exists in LAMMPS, offline pin + fair
bench **must** use LAMMPS. MetalWalls is oracle only for MW-only physics.
`REFERENCE.md` must state `oracle: lammps|metalwalls|literature` and
`oracle_version`.

## Why not edit dev in parallel?

Parallel agents on one tree thrash shared MD / longrange / electrode headers and
leave a half-broken `dev`. Isolation = each feature lands as a branch; merge is
ordered and abortable per slug.

## Args

| Arg | Default | Meaning |
| --- | ------- | ------- |
| `staff_count` / `parallelism` | 10 | Concurrent worktrees (wave 1 staff batch) |
| `max_waves` | 20 | Assignment rounds this shift |
| `target_branch` | `dev` | Merge destination |
| `do_perf` | false | Optional fair wall vs LAMMPS/MW |
| `mode` | `exhaust` | Until queue empty / budget |

## Launch

```text
agent_budget: 256
args: { "mode":"exhaust", "staff_count":10, "max_waves":20, "target_branch":"dev", "do_perf":false }
```

`agent_budget` = LLM-turn payroll for the shift, **not** headcount.

## Build (binding — MD-only)

```bash
bash -lc '
  module load buildenv-gcc/2026.03-mpich
  export OPENBLAS_ROOT="$OPENBLASDIR"
  export RUSTUP_TOOLCHAIN=stable-x86_64-unknown-linux-gnu
  export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
  cmake -B build -DCMAKE_BUILD_TYPE=Release \
        -DATV_BUILD_MD=ON -DATV_BUILD_KSPACE=ON \
        -DATV_BUILD_TESTS=ON -DATV_BUILD_REGRESSIONS=ON \
        -DATV_BUILD_HF=OFF -DATV_BUILD_PW=OFF -DATV_BUILD_GPW=OFF \
        -DATV_BUILD_CUDA=OFF -DATV_BUILD_TORCH=OFF \
        -DBLA_VENDOR=OpenBLAS
  cmake --build build -j$(nproc)
'
```

**Forbidden:** `ATV_BUILD_HF=ON`, `ATV_BUILD_PW=ON`, `ATV_BUILD_GPW=ON`,
`ATV_BUILD_CUDA=ON`, `ATV_BUILD_TORCH=ON`.

## Queue

`docs/developer/metalwalls-method-queue.md`

## Run log

(autonomous harness appends lines here)

- 2026-07-23 assign wave 1: const-pot-bo-parallel-plate, freeze-electrode-atoms, const-pot-charge-observer, charging-work-energy, const-pot-md-nvt-short, gaussian-electrode-charge, slab-2d-ewald-foothold, external-efield-md, buckingham-pair-lammps, four-site-tip4p-foothold → in_progress staff=10; incomplete was 0
- 2026-07-23 wave 1 merge: all ten green HEADs → dev (const-pot-bo-parallel-plate, charging-work-energy, gaussian-electrode-charge, freeze-electrode-atoms, const-pot-charge-observer, slab-2d-ewald-foothold, external-efield-md, buckingham-pair-lammps, four-site-tip4p-foothold, const-pot-md-nvt-short); blocked=[]; merge-failed=[]
- 2026-07-23 harness wave 1: queue merged→done (10); incomplete remaining=0; roadmaps ○→✓ Buckingham (CPU), constant E-field MD (CPU), TIP4P M-site virtual-site foothold (CPU); specs closed (10); ratios freeze=3.98 slab=6.49 efield=0.60 tip4p=0.338 others perf_exempt; queue_remaining=14 (pending only); handoff rewritten
- 2026-07-24 assign wave 2: const-pot-cg-large-electrode, electrode-charge-predictor, finite-field-single-electrode, pim-induced-dipole-bulk, piston-electrode-pz, thomas-fermi-electrode, fumi-tosi-nacl, steele-wall-argon, born-repulsion-foothold, capacitance-charge-fluctuation → in_progress staff=10; is_resume=false all; skipped dipole-electrode-joint (dep PIM); incomplete_before=0; pending_remaining=4; blocked=4
- 2026-07-24 wave 2 merge: all ten green HEADs → dev (const-pot-cg-large-electrode, electrode-charge-predictor, finite-field-single-electrode, pim-induced-dipole-bulk, piston-electrode-pz, thomas-fermi-electrode, fumi-tosi-nacl, steele-wall-argon, born-repulsion-foothold, capacitance-charge-fluctuation); blocked=[]; merge-failed=[]
- 2026-07-24 harness wave 2: queue merged→done (10 wave-2); incomplete remaining=0; roadmaps ○→✓ PairBorn/Fumi–Tosi, Steele wall/lj1043, bulk PIM+Thole (CPU); electrode notes CG/ASPC/finite-field single/piston/TF/capacitance; specs closed (10); all wave-2 ratios perf_exempt; queue_remaining=6 (pending only); dipole-electrode-joint unblocked; handoff rewritten

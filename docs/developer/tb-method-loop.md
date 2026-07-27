# TB / SE method loop — worktree → merge → dev

Abstract pattern: [`oracle-capability-loop.md`](oracle-capability-loop.md).  
Scaffold skill: `.agents/skills/create-oracle-loop/`.

## Continuous discovery (binding)

**Every wave starts with GapCensus:** re-diff OpenMOPAC / DFTB+ / xtb surfaces +
`docs/roadmap/semiempirical.md` ○ vs ATV on `target_branch`, refill pending,
update roadmap. Foothold ≠ capability closed. Stop `queue_empty` only when
`gaps_remaining=0` — **not** after one seed batch and **not** because a harness
wave marked a few slugs done.

**Family-locked oracles (binding):**

| Family | Oracle |
| ------ | ------ |
| NDDO (AM1/PM*/RM1/…) | **OpenMOPAC only** |
| DFTB0/2/3 / spin / PBC / TD-DFTB | **DFTB+ only** |
| GFN-xTB | **xtb only** (not DFTB3 successor; staff late) |
| Unified API shape (no E pin) | SCINE Sparrow (architecture) |
| Evaluate/MD embed patterns | CP2K (architecture only) |

**Forbidden:** OpenMOPAC as DFTB reference; DFTB+ as PM6/PM7 reference;
treating GFN-xTB as DFTB3++.

## Pipeline

```text
parent (dev)
  recover incomplete
       │
  for each wave:
    GapCensus → Assign → Worktrees → Merge? → Clearup → Harness? → Tutorial?
```

| Role | `isolation_worktree` | Writes |
| ---- | -------------------- | ------ |
| Census / foreman / merge / harness | **false** (parent) | `dev` only |
| Feature worker | **true** | private worktree only until merge |

Workers must report `worktree_path`, `branch` (`wf/<slug>`), `head_sha`.

**Recover incomplete first** (`## in_progress` `[~]`, open SE/TB specs,
unmerged `wf/<slug>`). Spec is source of truth. Offline pins only
(`tools/tb_oracle/`); never live oracle binaries in ctest.

**Scope IN:** `cpu/semiempirical` — DFTB0/2/3 ladder, NDDO PM depth, late
GFN-xTB; static in-tree SK/NDDO tables; Evaluate + Semiempirical Component.

**Scope OUT:** Gaussian SCF (PySCF loop); GPW/PW product (CP2K loop — arch-only
here); classical electrochem MD (MetalWalls); CUDA product; GFN-FF classical;
using OpenMOPAC to pin DFTB or DFTB+ to pin PM6.

**Architecture binding:** no public `*Config`; no `make_*` / `build_*` factories;
Frame is sole I/O; parameters are members / ctor args; static tables only.

## Assign ladder (binding)

DFTB0 → DFTB2 depth (analytic F, pins, 3ob) → DFTB3 → spin-DFTB → NDDO PM
ladder → **GFN-xTB last**.

## Args

| Arg | Default | Meaning |
| --- | ------- | ------- |
| `staff_count` / `parallelism` | 3–5 | Concurrent worktrees |
| `max_waves` | 20 | Assignment rounds this shift |
| `target_branch` | `dev` | Merge destination |
| `do_perf` | false | Optional fair wall vs oracle |
| `mode` | `exhaust` | Until gaps gone / budget |

## Build (binding — SE-only)

```bash
bash -lc '
  module load buildenv-gcc/2026.03-mpich
  export OPENBLAS_ROOT="$OPENBLASDIR"
  export RUSTUP_TOOLCHAIN=stable-x86_64-unknown-linux-gnu
  export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
  cmake -B build -DCMAKE_BUILD_TYPE=Release \
        -DATV_BUILD_TESTS=ON -DATV_BUILD_REGRESSIONS=ON \
        -DATV_BUILD_SEMIEMPIRICAL=ON -DATV_BUILD_EVALUATE=ON \
        -DATV_BUILD_MD=OFF -DATV_BUILD_MC=OFF -DATV_BUILD_MINIMIZER=OFF \
        -DATV_BUILD_VIB=OFF -DATV_BUILD_PATH=OFF -DATV_BUILD_POTENTIALS=OFF \
        -DATV_BUILD_NEIGHBOR=OFF -DATV_BUILD_HF=OFF -DATV_BUILD_PW=OFF \
        -DATV_BUILD_GPW=OFF -DATV_BUILD_KSPACE=OFF -DATV_BUILD_CUDA=OFF \
        -DATV_BUILD_TORCH=OFF -DBLA_VENDOR=OpenBLAS
  # SCOPED — never bare cmake --build build
  cmake --build build -j$(nproc) --target atomiverse <test_*> reg_<slug>
  ctest --test-dir build -R "^(test_...)$" --output-on-failure
  cmake --build build --target run_reg_<slug>
'
```

**Forbidden ON:** HF, PW, GPW, CUDA, TORCH, MD, MC, MINIMIZER, VIB, PATH,
KSPACE, POTENTIALS, NEIGHBOR.

**Do NOT build** `reg_am1_water_bomd` / `reg_am1_water_minimize` (need MD /
MINIMIZER). Domain suite (merge gate only):
`ctest -R '^test_cpu_(am1|rm1|dftb2|xtb_gfn0|nddo_)'`. Full unfiltered ctest
never.

## Queue / handoff

- Queue: `docs/developer/tb-method-queue.md`
- Handoff: `docs/developer/tb-method-handoff.md`
- Roadmap: `docs/roadmap/semiempirical.md`
- Offline oracle: `tools/tb_oracle/README.md`

## DFTB impl template (conceptual)

```text
ParameterSet → SK → H/S → eigen → P → SCC(Hubbard, γ) → repulsive → E/F/stress
```

## Temps

`/tmp/atv-*` `/tmp/tb_oracle_*` `/tmp/mopac_*` `/tmp/dftbplus_*` `/tmp/xtb_*`

## Run log

(autonomous harness appends lines here)

- 2026-07-27 GapCensus wave-1: re-diff roadmap ○ + OpenMOPAC/DFTB+/xtb surfaces vs `cpu/semiempirical` on `dev` @ `d80fd953`; pending refilled Phase A–C (23 capability slugs); blocked=9; incomplete=0; gaps_remaining=23+9 process/capability open; foothold ≠ closed
- 2026-07-27 assign wave 1: dftb2-3ob-organic, dftb2-analytic-forces, spin-dftb2-o2 → staffed (plus other ladder tips on `wf/*` retained incomplete)
- 2026-07-27 wave 1 merge: three green HEADs → dev (`dftb2-3ob-organic`, `dftb2-analytic-forces`, `spin-dftb2-o2`); fixups for test gates / CMake registration; merge-failed=[] for those three
- 2026-07-27 harness wave 1: queue three merged→done [x]; keep other incomplete [~] (dftb0-h2, dftb2-dftbplus-energy-pin, dftb3-water-mio, dftb2-dispersion-d3, dftb2-periodic-gamma); specs closed (3 deleted); roadmap ○→✓ DFTB analytic gradient (restricted) + spin-polarized DFTB (collinear); incomplete_remaining=5; queue_remaining=15 (pending only); blocked=9; handoff written; **queue not exhausted**

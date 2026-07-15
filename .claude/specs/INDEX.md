# Specs index

Breadth-first campaign (2026-07-12 → 2026-07-14): **complete**. Four probes +
one decision, all closed (history in git log; findings in `.claude/NOTES.md`;
the Frame verdict and its research report live durably at
`docs/developer/frame-block-redesign.md`).

Campaign **frame-block** — executing the signed Frame verdict (Frame =
mapping[string, Block], molrs-aligned schema, independent per-backend storage)
via the staged variant-C route. **Step 01 (registry + ColKey + mapping API +
typed meta) is done** (git log). Next: 02+ blockify directory by directory
(torch -> cuda -> cpu/md, minimizer, path -> cpu/potentials), deleting each
member once its directory is clear.

Campaign **mpi-launcher-embedded-python** (2026-07-15; 02/03/04 revised
2026-07-15 per user directives: AM1 single-point coverage, namespaced Python
API, ~/work/.x86_64 venv + molpy interop) — native MPI launcher: all ranks
start a native `atomiverse` executable, rank 0 alone embeds CPython and runs
`atv.py`; each `run()` is ONE task-granular in-memory broadcast (TaskSpec +
Frame) to worker ranks, no on-disk IR; the wheel ships both `_atomiverse.so`
and `bin/atomiverse`. Implementation order: 01 → 02 → 03 →
python-api-namespaces-am1 → 04 → 05. **Steps 01 (bcast_frame, done+deleted)
and 03 (protocol, done+deleted) are complete; 02 is code-complete** (git
log, 2026-07-15). python-api-namespaces-am1 and 04-embed are code-complete (the launcher
runs: `mpirun -np 2 atomiverse atv.py` ≡ `python atv.py`, LJ+AM1 parity
regression green). Next: `/mol:impl mpi-launcher-embedded-python-05-packaging`.

| slug | status | summary |
|---|---|---|
| [mpi-launcher-embedded-python-02-codec](mpi-launcher-embedded-python-02-codec.md) | code-complete | new `cpu/task/` cluster: versioned TaskSpec with TaskKind{MD, SinglePoint} (9 MD components + AM1EvaluatorParams), MPI-free pack/unpack, TaskMaterializer with task-kind-agnostic run(frame) |
| [python-api-namespaces-am1](python-api-namespaces-am1.md) | code-complete | deep canonical Python namespaces (atv.potentials.pair.LJCut) + 16-name top-level re-exports, flat layout deleted, AM1 single-point bindings, device= parameter (cpu-only v1), molpy interop in ~/work/.x86_64 venv |
| [mpi-launcher-embedded-python-04-embed](mpi-launcher-embedded-python-04-embed.md) | code-complete | `app/atomiverse/` executable (ATV_BUILD_LAUNCHER), rank-0 embedded CPython, PID-qualified LauncherGate (never throws), MD + AM1 single-point dispatch, LJ+AM1 parity regression |
| [mpi-launcher-embedded-python-05-packaging](mpi-launcher-embedded-python-05-packaging.md) | approved | wheel ships native `bin/atomiverse` via SKBUILD_SCRIPTS_DIR, plain cmake --install, runtime InterpreterLocator + ATV_PYTHON |

## Standing user directives (2026-07-12/14) — binding on every spec

1. **Frame is the only input/output currency.** No bare parallel scalar arrays.
2. **Frame = mapping[string, Block]** (signed 2026-07-14): molrs-aligned schema
   and keys, independent per-backend Column storage, bookkeeping in ONE typed
   `meta` region (string→MetaValue, extending molrs's string→string; Box rides
   simbox; revised 2026-07-14), dtype follows ATV precision tiers. Resolution
   stops at the launch boundary — in-loop key lookup is forbidden (13.7x) and
   machine-checked (debug resolution-count assertion + LJ ⚡ gate).
3. **N replicas come from an expansion strategy, not from the user.**

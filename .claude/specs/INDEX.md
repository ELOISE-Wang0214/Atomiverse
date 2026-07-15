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
log, 2026-07-15). **Campaign COMPLETE (2026-07-15)**: all six specs done and deleted (01-bcast,
02-codec, 03-protocol, python-api-namespaces-am1, 04-embed, 05-packaging;
history in git log). Delivered: `mpirun -np N atomiverse atv.py` ≡
`python atv.py` (LJ parity 7e-16, AM1 bitwise), namespaced Python API + AM1
single-point bindings, wheel ships `bin/atomiverse` with runtime interpreter
discovery. Plus one library fix en route: dual-image ghost double-counting
in the full-stencil neighbor path (fix(neighbor), f8c35b2).

| slug | status | summary |
|---|---|---|

## Standing user directives (2026-07-12/14) — binding on every spec

1. **Frame is the only input/output currency.** No bare parallel scalar arrays.
2. **Frame = mapping[string, Block]** (signed 2026-07-14): molrs-aligned schema
   and keys, independent per-backend Column storage, bookkeeping in ONE typed
   `meta` region (string→MetaValue, extending molrs's string→string; Box rides
   simbox; revised 2026-07-14), dtype follows ATV precision tiers. Resolution
   stops at the launch boundary — in-loop key lookup is forbidden (13.7x) and
   machine-checked (debug resolution-count assertion + LJ ⚡ gate).
3. **N replicas come from an expansion strategy, not from the user.**

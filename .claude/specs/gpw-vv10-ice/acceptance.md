# Acceptance — `gpw-vv10-ice`

## Gate result (wave 5 + wave 6 reconfirm)

| Gate | Result |
| ---- | ------ |
| Spec written | **yes** (this directory) |
| Offline CP2K pin | **no** — CP2K 2026.1 rejects `TYPE VV10` |
| Regression `reg_gpw-vv10-ice` | **not written** (no pin) |
| Implementation `GpwDft::vv10` | **not written** (would be false closure) |
| Dual gate green | **no** |
| Wave-6 re-probe | **ABORT** `invalid value for enumeration:VV10` (same as wave 5) |
| `acceptance_met` | **false** |

## Block reason (one line)

CP2K 2026.1 NON_LOCAL enum is only DRSLL|LMKLL|RVV10; stock VV10 has no
Quickstep oracle surface (probe ABORT `invalid value for enumeration:VV10`).

## Harness recommendation

Move queue row from `## in_progress` → `## blocked` as oracle-surface gap
(do not re-staff without unblock checklist). Keep roadmap `dft.md` §8 VV10 ○;
product VV-family nonlocal on GPW is rVV10 ✓ (`reg_gpw-rvv10-ice`).

## Unblock checklist

- [ ] CP2K (or approved alternate offline oracle) can produce a stock-VV10
      E_disp / E_c^{nl} pin for the ice water-dimer geometry
- [ ] Kernel tables or real-space Vydrov φ evaluation path chosen and
      embedded without runtime file I/O
- [ ] `GpwDft::vv10` exclusive flag + unit tests + `regressions/gpw/gpw-vv10-ice`
- [ ] Dual gate: pin tol + perf_exempt or ≤10× wall

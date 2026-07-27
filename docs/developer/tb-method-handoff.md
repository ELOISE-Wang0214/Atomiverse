# TB / SE method loop — handoff

Last assign: **2026-07-27 wave 3** on parent `dev` @ `223dcd4c`.
Last GapCensus: **2026-07-27 wave 3** (product surface still post harness wave-2
SE land; HEAD may include concurrent esmm/driver commits).
Last harness: **wave 2** merged `dftb2-dftbplus-energy-pin` +
`dftb2-periodic-gamma` → done.
**incomplete_remaining=8** (6 resume tips + 2 new staff). **pending=16**.
**blocked=9**. GapCensus still owns discovery — **queue not exhausted**.

## Workflow order (binding)

Recover incomplete first → Census → Assign (`incomplete ≫ pending`) →
Worktrees (**spec-first**) → Merge → Clearup → Harness → Tutorial. See
`docs/developer/tb-method-loop.md` and `docs/developer/oracle-capability-loop.md`.

## Assign wave 3 staff (2026-07-27)

| Slug | is_resume | Tip / note | oracle |
| ---- | --------- | ---------- | ------ |
| `dftb0-h2` | true | `wf/dftb0-h2` @ `bcff8c4a` ahead=1 behind=45 | dftb+ |
| `dftb3-water-mio` | true | `wf/dftb3-water-mio` @ `f26ba40b` ahead=1 behind=45 | dftb+ |
| `dftb2-dispersion-d3` | true | `wf/dftb2-dispersion-d3` @ `edabfa64` ahead=2 behind=45 | dftb+\|literature |
| `dftb2-analytic-forces-threeob` | true | `wf/dftb2-analytic-forces-threeob` @ `211dc386` ahead=1 behind=29 | dftb+\|literature |
| `dftb2-spin-analytic-forces` | true | `wf/dftb2-spin-analytic-forces` @ `1e8cac93` ahead=1 behind=29 | dftb+\|literature |
| `dftb2-analytic-stress` | true | `wf/dftb2-analytic-stress` @ `d2da3fe9` ahead=1 behind=29 | dftb+\|literature |
| `dftb2-periodic-forces` | false | **new** — lattice-image F under Γ-PBC | dftb+\|literature |
| `dftb2-kpoints` | false | **new** — k-mesh beyond Γ | dftb+ |

**Skipped:** `dftb3-analytic-forces` (dep on `dftb3-water-mio`); `td-dftb-h2`;
Phase B OpenMOPAC NDDO; Phase C GFN-xTB; blocked (9).

Workers must: isolation worktree on nobackup; `wf/<slug>`; **rebase resume
tips onto current `dev` first**; open `.claude/specs/<slug>/` before code;
SE-only cmake; offline pins only (`tools/tb_oracle/`); no push.

## Wave 2 land (2026-07-27) — historical

| Slug | Result |
| ---- | ------ |
| `dftb2-dftbplus-energy-pin` | **merged → done** (mio-1-1 H–H absolute E) |
| `dftb2-periodic-gamma` | **merged → done** (Γ lattice-sum energy) |

Roadmap ○→✓: Periodic DFTB (Γ-point foothold). k-mesh / PBC forces / stress /
DFTB0/3/D3 remain ○ or incomplete.

## Wave 1 land (2026-07-27) — historical

| Slug | Result |
| ---- | ------ |
| `dftb2-3ob-organic` | **merged → done** |
| `dftb2-analytic-forces` | **merged → done** |
| `spin-dftb2-o2` | **merged → done** |

## Oracle / build (binding)

| Knob | Value |
| ---- | ----- |
| Family oracles | NDDO → **OpenMOPAC**; DFTB* → **DFTB+**; GFN-xTB → **xtb** |
| Forbidden | OpenMOPAC as DFTB pin; DFTB+ as PM6/PM7 pin; GFN as DFTB3 successor |
| Offline only | `tools/tb_oracle/`; never live oracle in ctest |
| SE-only cmake | `SEMIEMPIRICAL`+`EVALUATE`+TESTS+REGRESSIONS **ON**; HF/PW/GPW/MD/MC/MINIMIZER/VIB/PATH/KSPACE/POTENTIALS/NEIGHBOR/CUDA/TORCH **OFF** |
| Domain suite (merge gate) | `ctest -R '^test_cpu_(am1\|rm1\|dftb2\|xtb_gfn0\|nddo_)'` |
| Worker build | scoped `--target atomiverse <test_*> reg_<slug>` only |

## Next

1. **Worktrees** for the 8 staff (parallelism ≤ staff_count).
2. Resumes: **rebase** `wf/*` onto `dev`, re-open/repair spec if missing, dual gate.
3. New: spec → offline pin → impl → dual gate → commit on `wf/<slug>`.
4. Merge greens → Clearup (same wave) → Harness → GapCensus wave 4.

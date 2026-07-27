# TB / SE method loop — handoff

Last harness pass: **2026-07-27 wave 1** on `dev` (`chore(harness)`).
Last GapCensus: **2026-07-27 wave 2** on `dev` @ `77d80486`.
Merged greens (wave 1): `dftb2-3ob-organic`, `dftb2-analytic-forces`, `spin-dftb2-o2`.
**incomplete_remaining=5** (unmerged TB `wf/*` tips kept `[~]`, behind≈12).
**pending=20** (wave-2 refilled +5 deepen). **blocked=9**.
GapCensus still owns discovery — **queue not exhausted**.

**Workflow order (binding):** Recover incomplete first → Census →
Assign (`incomplete ≫ pending`) → Worktrees (**spec-first**) → Merge → Clearup →
Harness → Tutorial. See `docs/developer/tb-method-loop.md` and
`docs/developer/oracle-capability-loop.md`.

## Wave 1 land (2026-07-27)

| Slug | Component / change | Regression | oracle | Result |
| ---- | ------------------ | ---------- | ------ | ------ |
| `dftb2-3ob-organic` | `DftbSkSet::ThreeOb` + `basis_data/dftb_3ob_tables.h` (3ob-3-1 CHNO) | `regressions/semiempirical/dftb2-3ob-organic` | dftb+ 24.1 | **merged → done** |
| `dftb2-analytic-forces` | Closed-form restricted nuclear F (HF+Pulay+γ+rep) on MioSketch | `regressions/semiempirical/dftb2-analytic-forces` | dftb+ (formula) | **merged → done** |
| `spin-dftb2-o2` | Collinear unrestricted spin channel; O₂ triplet m | `regressions/semiempirical/spin-dftb2-o2` | literature + family dftb+ | **merged → done** |

**Specs closed (deleted):** `dftb2-3ob-organic`, `dftb2-analytic-forces`,
`spin-dftb2-o2` under `.claude/specs/` (history = git log).

**Roadmap ○→✓:** DFTB analytic gradient (restricted closed-form); spin-polarized
DFTB (collinear foothold). 3ob organic SK is prose/product depth on DFTB2 ✓
(not a separate roadmap row; DFTB3 remains ○). **No claim** that Phase A is
closed — DFTB0/3/D3/PBC/stress/TD and mio offline pins remain.

**Perf:** all three `perf_exempt` (no fair DFTB+ wall on pin host for ratio).

## Incomplete retained `[~]` (do not demote)

| Slug | Branch tip | Notes |
| ---- | ---------- | ----- |
| `dftb0-h2` | `wf/dftb0-h2` | Non-SCC foothold commits ahead of `dev` |
| `dftb2-dftbplus-energy-pin` | `wf/dftb2-dftbplus-energy-pin` | Offline pin / SK hygiene |
| `dftb3-water-mio` | `wf/dftb3-water-mio` | DFTB3 energy foothold tip |
| `dftb2-dispersion-d3` | `wf/dftb2-dispersion-d3` | D3(BJ) water-dimer tip |
| `dftb2-periodic-gamma` | `wf/dftb2-periodic-gamma` | Γ-point PBC tip |

Clearup already wiped isolation worktrees; branch tips remain until merge or
explicit drop. **Next wave Recover** must prefer these over new pending staff
(`incomplete ≫ pending`).

## Prior census (wave-1 GapCensus, pre-land)

See `docs/developer/tb-method-queue.md` § GapCensus — 2026-07-27 wave-1.
Refilled Phase A–C pending from roadmap ○ + oracle surfaces; foothold ≠ closed.

## Oracle / build (binding)

| Knob | Value |
| ---- | ----- |
| Family oracles | NDDO → **OpenMOPAC**; DFTB* → **DFTB+**; GFN-xTB → **xtb** |
| Forbidden | OpenMOPAC as DFTB pin; DFTB+ as PM6/PM7 pin; GFN as DFTB3 successor |
| Offline only | `tools/tb_oracle/`; never live oracle in ctest |
| SE-only cmake | `SEMIEMPIRICAL`+`EVALUATE`+TESTS+REGRESSIONS **ON**; HF/PW/GPW/MD/MC/MINIMIZER/VIB/PATH/KSPACE/POTENTIALS/NEIGHBOR/CUDA/TORCH **OFF** |
| Domain suite (merge gate) | `ctest -R '^test_cpu_(am1\|rm1\|dftb2\|xtb_gfn0\|nddo_)'` |
| Worker build | scoped `--target atomiverse <test_*> reg_<slug>` only |

## Public surface landed (wave 1)

```cpp
cpu::semiempirical::DFTB2 dftb;
dftb.sk_set = cpu::semiempirical::DftbSkSet::ThreeOb;  // or MioSketch
dftb.multiplicity = 3;   // collinear open shell (MioSketch)
dftb.with_forces = true; // analytic F when restricted MioSketch
// Frame I/O via Evaluate + Semiempirical; no *Config; no make_*/build_*
```

## GapCensus wave 2 (2026-07-27)

Re-diffed family oracles + roadmap ○ + tree on `dev`. Wave-1 footholds do **not**
close ThreeOb/open-shell forces, k-mesh, DFTB3 forces, or hard mopac RM1 pins —
pending refilled (**+5** slugs). `tools/tb_oracle/README.md` added (was missing
but referenced). No roadmap ○→✓ this wave (nothing new merged).

## Next

1. **Recover** the five `[~]` tips (**rebase onto `dev` first** — behind≈12).
2. Prefer incomplete ≫ new Phase A pending on Assign.
3. Ladder remaining: DFTB0/3 → D3/PBC/k/stress/TD → force deepen → NDDO PM → GFN late.

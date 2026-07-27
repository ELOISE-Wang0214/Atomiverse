# TB / SE method loop — handoff

Last harness pass: **2026-07-27 wave 2** on `dev` (`chore(harness)`).
Last GapCensus: **2026-07-27 wave 2** on `dev` @ `77d80486` (pre-land snapshot;
landed tips updated below).
Merged greens (wave 2): `dftb2-dftbplus-energy-pin`, `dftb2-periodic-gamma`.
**incomplete_remaining=3** (unmerged TB `wf/*` tips kept `[~]`, behind≈39).
**pending=20** (wave-2 refill unchanged this harness). **blocked=9**.
GapCensus still owns discovery — **queue not exhausted**.

## Workflow order (binding)

Recover incomplete first → Census → Assign (`incomplete ≫ pending`) →
Worktrees (**spec-first**) → Merge → Clearup → Harness → Tutorial. See
`docs/developer/tb-method-loop.md` and `docs/developer/oracle-capability-loop.md`.

## Wave 2 land (2026-07-27)

| Slug | Component / change | Regression | oracle | Result |
| ---- | ------------------ | ---------- | ------ | ------ |
| `dftb2-dftbplus-energy-pin` | mio-1-1 H–H static SK (`dftb_mio_hh_sk.h`) + H₂ total-E pin | `regressions/semiempirical/dftb2-dftbplus-energy-pin` | dftb+ 24.1 | **merged → done** |
| `dftb2-periodic-gamma` | Γ-point lattice-sum H⁰/S, screened γ, repulsive (MioSketch) | `regressions/semiempirical/dftb2-periodic-gamma` | literature (+ family dftb+) | **merged → done** |

**Specs closed (deleted):** `dftb2-dftbplus-energy-pin`, `dftb2-periodic-gamma`
under `.claude/specs/` (history = git log).

**Roadmap ○→✓:** Periodic DFTB (Γ-point foothold). H–H absolute energy pin is
hygiene prose (not a separate roadmap row). k-mesh / stress / DFTB0/3/D3 remain
○ or incomplete tips. **No claim** that Phase A is closed.

**Perf:** both wave-2 regs `perf_exempt` / literature twin where host lacked
fair DFTB+ wall (documented in REFERENCE).

## Wave 1 land (2026-07-27) — historical

| Slug | Result |
| ---- | ------ |
| `dftb2-3ob-organic` | **merged → done** |
| `dftb2-analytic-forces` | **merged → done** |
| `spin-dftb2-o2` | **merged → done** |

## Incomplete retained `[~]` (do not demote)

| Slug | Branch tip | Notes |
| ---- | ---------- | ----- |
| `dftb0-h2` | `wf/dftb0-h2` @ `bcff8c4a` | Non-SCC foothold; ahead=1 behind≈39 |
| `dftb3-water-mio` | `wf/dftb3-water-mio` @ `f26ba40b` | DFTB3 energy foothold; ahead=1 behind≈39 |
| `dftb2-dispersion-d3` | `wf/dftb2-dispersion-d3` @ `edabfa64` | D3(BJ) water-dimer tip; ahead=2 behind≈39 |

Clearup already wiped isolation worktrees; branch tips remain until merge or
explicit drop. **Next wave Recover** must prefer these over new pending staff
(`incomplete ≫ pending`). Rebase onto current `dev` first.

## Oracle / build (binding)

| Knob | Value |
| ---- | ----- |
| Family oracles | NDDO → **OpenMOPAC**; DFTB* → **DFTB+**; GFN-xTB → **xtb** |
| Forbidden | OpenMOPAC as DFTB pin; DFTB+ as PM6/PM7 pin; GFN as DFTB3 successor |
| Offline only | `tools/tb_oracle/`; never live oracle in ctest |
| SE-only cmake | `SEMIEMPIRICAL`+`EVALUATE`+TESTS+REGRESSIONS **ON**; HF/PW/GPW/MD/MC/MINIMIZER/VIB/PATH/KSPACE/POTENTIALS/NEIGHBOR/CUDA/TORCH **OFF** |
| Domain suite (merge gate) | `ctest -R '^test_cpu_(am1\|rm1\|dftb2\|xtb_gfn0\|nddo_)'` |
| Worker build | scoped `--target atomiverse <test_*> reg_<slug>` only |

## Public surface landed (through wave 2)

```cpp
cpu::semiempirical::DFTB2 dftb;
dftb.sk_set = cpu::semiempirical::DftbSkSet::ThreeOb;  // or MioSketch
dftb.multiplicity = 3;   // collinear open shell (MioSketch)
dftb.with_forces = true; // analytic F when restricted MioSketch (molecular)
// Γ-point PBC: lattice sums active under finite box (MioSketch)
// Frame I/O via Evaluate + Semiempirical; no *Config; no make_*/build_*
// H2 mio-1-1 absolute E: reg_dftb2-dftbplus-energy-pin
```

## Next

1. **Recover** the three `[~]` tips (**rebase onto `dev` first** — behind≈39).
2. Prefer incomplete ≫ new Phase A pending on Assign.
3. Ladder remaining: DFTB0/3 → D3/k/stress/TD force deepen → NDDO PM → GFN late.
4. GapCensus owns further discovery — **never claim queue exhausted**.

# TB / SE method loop — handoff

Last assign: **2026-07-27 wave 3** on parent `dev` @ `223dcd4c`.
Last GapCensus: **2026-07-27 wave 3**.
Last harness: **wave 3** (2026-07-29) — all **8** assign-w3 staff → **done**
on `dev` @ `30672ccf` (final batch: `dftb3-water-mio`, `dftb2-periodic-forces`,
`dftb2-kpoints`; earlier wave-3 landings: DFTB0, D3, ThreeOb F, spin F, stress).
**incomplete_remaining=0**. **pending=16**. **blocked=9**.
GapCensus still owns discovery — **queue not exhausted**.

## Workflow order (binding)

Recover incomplete first → Census → Assign (`incomplete ≫ pending`) →
Worktrees (**spec-first**) → Merge → Clearup → Harness → Tutorial. See
`docs/developer/tb-method-loop.md` and `docs/developer/oracle-capability-loop.md`.

## Harness wave 3 land (2026-07-29)

| Slug | Result |
| ---- | ------ |
| `dftb0-h2` | **merged → done** (non-SCC DFTB0; DFTB+ 24.1 H–H pin) |
| `dftb3-water-mio` | **merged → done** (third-order SCC; dual-gate re-pin) |
| `dftb2-dispersion-d3` | **merged → done** (molecular D3(BJ) + analytic F) |
| `dftb2-analytic-forces-threeob` | **merged → done** (closed-form ThreeOb F) |
| `dftb2-spin-analytic-forces` | **merged → done** (collinear open-shell F) |
| `dftb2-analytic-stress` | **merged → done** (restricted MioSketch virial) |
| `dftb2-periodic-forces` | **merged → done** (lattice-image Γ-PBC F) |
| `dftb2-kpoints` | **merged → done** (Monkhorst–Pack multi-k energy) |

Specs closed (deleted): `dftb3-water-mio`, `dftb2-periodic-forces`,
`dftb2-kpoints`. TB `wf/*` tips and isolation worktrees: **none**.

Roadmap prose: PBC lattice-image forces + multi-k energy + ThreeOb analytic F
documented as shipped; table rows already ✓ for DFTB0/3/D3/stress/gradient/
periodic (no new table ○→✓ this harness beyond prose deepen). **Still open:**
multi-k F/stress, periodic D3 lattice sum, `dftb3-analytic-forces`, TD-DFTB,
Phase B OpenMOPAC NDDO, Phase C GFN product.

## Assign wave 3 staff (2026-07-27) — historical

| Slug | is_resume | Note | oracle |
| ---- | --------- | ---- | ------ |
| `dftb0-h2` | true | landed | dftb+ |
| `dftb3-water-mio` | true | landed | dftb+ |
| `dftb2-dispersion-d3` | true | landed | dftb+\|literature |
| `dftb2-analytic-forces-threeob` | true | landed | dftb+\|literature |
| `dftb2-spin-analytic-forces` | true | landed | dftb+\|literature |
| `dftb2-analytic-stress` | true | landed | dftb+\|literature |
| `dftb2-periodic-forces` | false | landed | dftb+\|literature |
| `dftb2-kpoints` | false | landed | dftb+ |

## Wave 2 land (2026-07-27) — historical

| Slug | Result |
| ---- | ------ |
| `dftb2-dftbplus-energy-pin` | **merged → done** (mio-1-1 H–H absolute E) |
| `dftb2-periodic-gamma` | **merged → done** (Γ lattice-sum energy) |

## Wave 1 land (2026-07-27) — historical

| slug | Result |
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

1. **GapCensus wave 4** on `dev` (re-diff OpenMOPAC|DFTB+|xtb + roadmap ○ +
   `cpu/semiempirical`); foothold ≠ closed — refill pending if deepen needed.
2. Assign next staff from pending (prefer Phase A `dftb3-analytic-forces` /
   `td-dftb-h2` only after census; then Phase B mopac ladder).
3. No incomplete `[~]` tips to recover.
4. **Queue not exhausted** — GapCensus owns discovery.

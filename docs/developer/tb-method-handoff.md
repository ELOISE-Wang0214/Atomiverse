# TB / SE method loop — handoff

Last assign: **2026-07-27 wave 3** on parent `dev` @ `223dcd4c` (historical).
Last GapCensus: **2026-07-29 wave 4** on parent `dev` @ `2cf4b935`.
Last harness: **wave 3** (2026-07-29) — all **8** assign-w3 staff → **done**
on `dev` @ `30672ccf`…`2cf4b935` (final batch: `dftb3-water-mio`,
`dftb2-periodic-forces`, `dftb2-kpoints`; earlier: DFTB0, D3, ThreeOb F,
spin F, stress). **incomplete_remaining=0**. **pending=20**. **blocked=9**.
`new_slugs_added=4` (wave-4 deepen). GapCensus still owns discovery —
**queue not exhausted**.

## Workflow order (binding)

Recover incomplete first → Census → Assign (`incomplete ≫ pending`) →
Worktrees (**spec-first**) → Merge → Clearup → Harness → Tutorial. See
`docs/developer/tb-method-loop.md` and `docs/developer/oracle-capability-loop.md`.

## GapCensus wave 4 (2026-07-29)

| Item | Result |
| ---- | ------ |
| Parent | `dev` @ `2cf4b935` |
| Incomplete TB `wf/*` | **none** |
| Isolation worktrees (TB) | **none** (only unrelated metalwalls `repin-*` trees) |
| Live SE/TB specs | **none** |
| Pending rewrite | Phase A: keep `dftb3-analytic-forces`, `td-dftb-h2`; **add** `dftb2-multi-k-forces`, `dftb2-dispersion-d3-periodic`, `dftb2-analytic-stress-threeob`, `dftb2-spin-analytic-stress`; Phase B/C unchanged |
| Blocked | **9** (unchanged; scope-out never pending) |
| Roadmap | prose deepen (multi-k F, periodic D3, stress deepens); table ○/✓ **unchanged** (no new product merge this census) |

### Foothold ≠ closed (staff these next)

1. **DFTB3 analytic F** — energy green; forces = re-SCC FD
2. **Multi-k F/stress** — MP energy green; F/stress throw
3. **Periodic D3** — molecular D3 green; lattice sum absent on HEAD
4. **ThreeOb / spin analytic stress** — F green; stress still strain-FD
5. Then TD-DFTB → Phase B OpenMOPAC → Phase C xtb (ladder)

## Harness wave 3 land (2026-07-29)

| slug | Result |
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

1. **Assign wave 4** — staff Phase A deepens (DFTB3 F, multi-k F, periodic D3,
   stress deepens); incomplete empty so no Recover first.
2. Spec-first workers on `wf/<slug>`; SE-only cmake; offline DFTB+ pins only.
3. Do **not** staff Phase C GFN while Phase A TD/DFTB3 F open unless
   `theme=gfn|xtb`.
4. mopac/xtb `tools/tb_oracle/inputs/` still empty — seed with first hard pin.
5. **No push.** Wipe builds after dual-gate; Clearup worktrees same wave.

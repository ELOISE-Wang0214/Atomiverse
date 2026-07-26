# CP2K method loop — worktree → merge → dev

Abstract pattern: [`oracle-capability-loop.md`](oracle-capability-loop.md).  
Scaffold skill: `.agents/skills/create-oracle-loop/`.

## Continuous discovery (binding)

**Every wave starts with GapCensus:** re-diff CP2K GPW/PW/periodic surface +
`docs/roadmap/dft.md` vs ATV on `target_branch`, refill `## pending`, update
roadmap. Foothold ≠ capability closed. Stop `queue_empty` only when
`gaps_remaining=0` — **not** after one seed batch.

## Pipeline

```text
parent (dev)
  recover incomplete
       │
  for each wave:
    GapCensus → Assign → Worktrees → Merge? → Clearup → Harness? → Tutorial?
```

Instance = template + CONFIG/`extra_law` only
(`.grok/workflows/cp2k-method-loop.rhai`).

| Role | `isolation_worktree` | Writes |
| ---- | -------------------- | ------ |
| Recover / GapCensus / assign / merge / harness | **false** (parent) | `dev` only |
| Feature worker | **true** | private worktree only until merge |

Workers must report `worktree_path`, `branch` (`wf/<slug>`), `head_sha`.

**Recover incomplete first** (`## in_progress` `[~]`, open GPW/PW specs,
unmerged `wf/<slug>`). Spec is source of truth. Offline CP2K pins only
(`tools/cp2k_oracle/`); never live CP2K in ctest.

**Scope:** GPW / PW / periodic DFT. **Not** molecular post-HF (PySCF loop) or
classical MD electrochemistry (MetalWalls loop).

## Why not edit dev in parallel?

Parallel agents on one tree thrash shared GPW/PW headers and leave a half-broken
`dev`. Isolation = each feature lands as a branch; merge is ordered and
abortable per slug.

## Args

| Arg | Default | Meaning |
| --- | ------- | ------- |
| `staff_count` / `parallelism` | 10 | Concurrent worktrees |
| `max_waves` | 20 | Assignment rounds this shift |
| `target_branch` | `dev` | Merge destination |
| `do_perf` | true | Fair wall vs CP2K when available |
| `mode` | `exhaust` | Until gaps gone / budget |

## Launch

```text
agent_budget: 1024
args: { "mode":"exhaust", "staff_count":10, "max_waves":20, "target_branch":"dev", "do_perf":true }
```

## Queue

`docs/developer/cp2k-method-queue.md`

## Run log

(autonomous harness appends lines here)

- **2026-07-26 harness wave 7 (parent `dev` @ `8cb47664`):** 9 merged →
  done: `periodic-urpa-si`, `periodic-ugw-kpoints`, `gpw-uks-smearing`,
  `gpw-hirshfeld-uks`, `cpmd-uks`, `second-gen-cpmd-uks`, `ot-kpoints`,
  `ot-dftu`, `dfpt-spin-si`. Incomplete remaining **1** (`[~] gpw-vv10-ice`
  on unmerged `wf/gpw-vv10-ice`); pending residual **0** (census list
  emptied — **not** queue exhausted); blocked **5**. Specs closed
  (deleted) for all 9. Roadmap: **no new ○→✓ table flips**; residual
  prose closed URPA@UKS, multi-k UGW, GPW UKS smear, UKS Hirshfeld,
  CPMD/SGCP UKS, OT+k, OT+DFT+U, collinear-spin DFPT; stock VV10 stays
  open. Clearup already ran (no isolation worktrees; no builds).
  **Queue not exhausted** — next wave Recover incomplete then GapCensus.
  No push.
- **2026-07-26 GapCensus wave 7 (parent `dev` @ `9bdeddf1`):** recover
  incomplete=**1** (`[~] gpw-vv10-ice` on unmerged `wf/gpw-vv10-ice` —
  not demoted); CP2K specs on `dev`=0; isolation worktrees=0;
  builds=none. Prior pending **1** (`periodic-urpa-si`) re-verified;
  **+8** residuals → pending **9** (`periodic-urpa-si`,
  `periodic-ugw-kpoints`, `gpw-uks-smearing`, `gpw-hirshfeld-uks`,
  `cpmd-uks`, `second-gen-cpmd-uks`, `ot-kpoints`, `ot-dftu`,
  `dfpt-spin-si`); in_progress **1**; blocked **5**. No ○/✓ table flips;
  §8/§10 residual prose + index GPW residual line + MD AIMD UKS residual
  notes refreshed. Next: Assign prefers incomplete first (VV10
  docs/block resolution) then surface-disjoint pending (URPA ∥ multi-k
  UGW careful share with GW/RPA kernels; UKS smear ∥ UKS Hirshfeld free
  of each other; CPMD UKS ∥ SGCP UKS careful AIMD surface; OT+k ∥ OT+DFT+U
  careful OT tree share; DFPT-spin free of GPW).
- **2026-07-26 harness wave 6 (parent `dev` @ `0869833a`):** 9 merged →
  done: `gpw-metallic-smearing`, `gpw-xdm-ice`, `gpw-uks-rpbe`,
  `gpw-dftu-kpoint-stress`, `periodic-ump2-kpoints`, `periodic-ugw-si`,
  `gpw-hirshfeld-kpoints`, `dfpt-gga-si`, `xl-bomd-uks`. Incomplete
  remaining **1** (`[~] gpw-vv10-ice` on unmerged `wf/gpw-vv10-ice`);
  pending residual **1** (`periodic-urpa-si` — **not** queue exhausted);
  blocked **5**. Specs closed (deleted) for all 9. Roadmap: **XDM ✓**;
  residual prose closed GPW metallic smear, Γ-UKS RPBE, multi-k DFT+U
  stress, multi-k UMP2, UGW, multi-k Hirshfeld, DFPT-GGA, XL-BOMD UKS;
  stock VV10 / URPA stay open. Clearup already ran (no isolation
  worktrees; no builds). **Queue not exhausted** — next wave Recover
  incomplete then GapCensus. No push.
- **2026-07-26 GapCensus wave 6 (parent `dev` @ `5847a35a`):** recover
  incomplete=**2** (`[~] gpw-metallic-smearing`, `[~] gpw-vv10-ice` on
  unmerged `wf/*` — not demoted); CP2K specs on `dev`=0; isolation
  worktrees=0; builds=none. Prior pending **1** (`gpw-xdm-ice`) re-
  verified; **+8** residuals → pending **9** (`gpw-xdm-ice`,
  `gpw-uks-rpbe`, `gpw-dftu-kpoint-stress`, `periodic-ump2-kpoints`,
  `periodic-ugw-si`, `periodic-urpa-si`, `gpw-hirshfeld-kpoints`,
  `dfpt-gga-si`, `xl-bomd-uks`); in_progress **2**; blocked **5**. No
  ○/✓ table flips; §3/§8/§10–12 residual prose + index GPW residual line
  refreshed. Next: Assign prefers incomplete first (metallic merge gate;
  VV10 docs/block resolution) then surface-disjoint pending (one nonlocal
  XDM; UKS RPBE free; multi-k DFT+U stress ∥ multi-k UMP2 careful share
  with RI kernels; UGW/URPA after Γ UMP2 ✓ free of closed-shell GW/RPA
  surface; Hirshfeld×k ∥ DFPT-GGA ∥ XL-BOMD UKS free).
- **2026-07-26 harness wave 5 (parent `dev` @ `2d43d1c2`):** 6 merged →
  done: `gpw-dftu-kpoints`, `periodic-gw-finite-q`, `gpw-uks-revpbe`,
  `cpmd-rattle-water`, `periodic-ump2-si`, `ot-uks-o2`. Incomplete
  remaining **2** (`[~] gpw-metallic-smearing`, `[~] gpw-vv10-ice` on
  unmerged `wf/*`); pending residual **1** (`gpw-xdm-ice` — **not**
  queue exhausted); blocked **5**. Specs closed (deleted) for all 6.
  Roadmap: no new ○→✓ table flips; residual prose closed multi-k DFT+U,
  finite-q G0W0 \(W\), Γ-UKS revPBE, CPMD electronic RATTLE, UMP2,
  OT-UKS; XDM / stock VV10 / GPW metallic smear stay open. Clearup
  already ran (no isolation worktrees; no builds). **Queue not
  exhausted** — next wave Recover incomplete then GapCensus. No push.
- **2026-07-26 GapCensus wave 5 (parent `dev` @ `4437d80f`):** recover
  clean (incomplete=0, CP2K specs=0, `wf/*`=0, builds=none). Prior pending
  residual **0** (harness emptied list — not capability-closed). **+9**
  residuals → pending **9** (`gpw-dftu-kpoints`, `periodic-gw-finite-q`,
  `gpw-uks-revpbe`, `gpw-xdm-ice`, `gpw-vv10-ice`, `cpmd-rattle-water`,
  `gpw-metallic-smearing`, `periodic-ump2-si`, `ot-uks-o2`); in_progress
  **0**; blocked **5**. No ○/✓ table flips; §10 residual prose + index
  Cell multi-k UKS line refreshed (stress productized; multi-k DFT+U open).
  Next: Assign (prefer surface-disjoint: multi-k DFT+U ∥ finite-q GW
  careful share with RPA kernels ∥ one nonlocal (XDM or VV10)/wave ∥
  UKS revPBE ∥ CPMD RATTLE free of OT-UKS / UMP2).
- **2026-07-26 harness wave 4 (parent `dev` @ `85eaa2d9`):** 9 merged →
  done: `gpw-kpoint-uks-stress`, `periodic-gw-kpoints`, `gpw-vdw-df2-ice`,
  `gpw-rvv10-ice`, `gpw-mbd-periodic`, `pw-spin-gga`, `gpw-uks-blyp`,
  `periodic-rpa-finite-q`, `gpw-dftu-stress`. Incomplete remaining **0**;
  pending residual **0** (census list emptied — **not** queue exhausted);
  blocked **5**. Specs closed (deleted) for all 9. Roadmap: vdW-DF2 ✓
  (LMKLL), rVV10 ✓; multi-k G0W0 residual closed; periodic Ewald MBD,
  PW true spin-GGA, Γ-UKS BLYP, multi-k UKS stress, finite-q inter-k
  dRPA, DFT+U analytic stress productized. Clearup already ran (no
  isolation worktrees; no builds). **Queue not exhausted** — next wave
  GapCensus. No push.
- **2026-07-26 GapCensus wave 4 (parent `dev` @ `decf917f`):** recover
  clean (incomplete=0, CP2K specs=0, `wf/*`=0, builds=none). Prior pending
  residual **0** (harness emptied list — not capability-closed). **+9**
  residuals → pending **9** (`gpw-kpoint-uks-stress`, `periodic-gw-kpoints`,
  `gpw-vdw-df2-ice`, `gpw-rvv10-ice`, `gpw-mbd-periodic`, `pw-spin-gga`,
  `gpw-uks-blyp`, `periodic-rpa-finite-q`, `gpw-dftu-stress`); in_progress
  **0**; blocked **5**. No ○/✓ table flips; §8/§10 residual prose + index
  GPW residual line refreshed. Next: Assign (prefer surface-disjoint:
  stress/DFT+U ∥ nonlocal (one VdW slug/wave) ∥ multi-k GW/RPA careful
  surface share; PW spin-GGA ∥ GPW UKS BLYP free).
- **2026-07-25 harness wave 3 (parent `dev` @ `0a642122`):** 8 merged →
  done: `gpw-rpbe-water`, `gpw-kpoint-uks-forces`, `gpw-uks-stress`,
  `gpw-vdw-df-ice`, `gpw-mbd-hirshfeld`, `periodic-rpa-kpoints`,
  `periodic-gw-si`, `dfpt-incommensurate-si`. Incomplete remaining **0**;
  pending residual **0** (census list emptied — **not** queue exhausted);
  blocked **5**. Specs closed (deleted) for all 8. Roadmap: RPBE ✓,
  vdW-DF ✓ (DRSLL), Periodic GW ✓ (Γ spectral G0W0; multi-k residual);
  Hirshfeld-scaled MBD/TS; multi-k dRPA residual closed; Γ-UKS stress +
  multi-k UKS forces + incommensurate DFPT productized (multi-k UKS stress
  residual). Clearup already ran (no isolation worktrees; no builds).
  **Queue not exhausted** — next wave GapCensus. No push.
- **2026-07-25 GapCensus wave 3 (parent `dev` @ `def040bd`):** recover
  clean (incomplete=0, CP2K specs=0, `wf/*`=0, builds=none). Prior pending
  **3** re-verified true gaps; **+5** residuals (`gpw-kpoint-uks-forces`,
  `gpw-uks-stress`, `periodic-rpa-kpoints`, `gpw-mbd-hirshfeld`,
  `dfpt-incommensurate-si`) → pending **8**; in_progress **0**; blocked
  **5**. No ○/✓ table flips; open-follow-on / forces / stress / RPA / DFPT
  prose + index Cell “strain-FD” stale line refreshed. Next: Assign
  (prefer surface-disjoint staff: RPBE alone; k-UKS forces ∥ UKS stress
  careful surface share; multi-k RPA after Γ RPA ✓ free).
- **2026-07-25 harness wave 2 (parent `dev` @ `f5739506`):** 10 merged →
  done: `gpw-kpoint-forces`, `gpw-revpbe-si`, `gpw-kpoints-uks`,
  `gpw-uks-forces`, `gpw-kpoint-stress`, `gpw-mbd-ice`,
  `periodic-mp2-kpoints`, `periodic-rpa-si`, `dfpt-finite-q-si`,
  `gpw-variable-cell-molopt`. Incomplete remaining **0**; pending residual
  **3** (`gpw-rpbe-water`, `gpw-vdw-df-ice`, `periodic-gw-si`); blocked **5**.
  Specs closed (deleted) for all 10. Roadmap: revPBE ✓, free-atom MBD ✓,
  Periodic RPA ✓ (Γ dRPA@KS), finite-q DFPT productized; multi-k MP2 residual
  closed; GPW k-UKS / k-forces / k-stress / Γ-UKS forces / MOLOPT CELL_OPT
  prose. Clearup already ran (no isolation worktrees; no builds). **Queue not
  exhausted** — next wave GapCensus. No push.
- **2026-07-25 GapCensus wave 2 (parent `dev` @ `507bb5f8`):** recover found
  unmerged `wf/gpw-kpoint-forces` (ahead 1) → `## in_progress` `[~]`. Prior
  pending (−k-forces) re-verified true gaps; **+4** follow-ons
  (`dfpt-finite-q-si`, `periodic-mp2-kpoints`, `gpw-kpoint-stress`,
  `gpw-uks-forces`) → pending **12**; in_progress **1**; blocked **5**. No
  ○/✓ table flips; stale DFPT/k-point intro prose refreshed. Next: Assign
  (prefer resume k-forces + surface-disjoint pending).
- **2026-07-25 harness wave 1 (post GapCensus staff; parent `dev` @
  `83a100df`):** 9 merged → done: `gpw-blyp-water`, `gpw-uks-pbe`,
  `pw-uks-fe-smearing`, `pw-variable-cell-si`, `gpw-ts-ice`,
  `phonon-dispersion-si`, `periodic-wannier-si`, `second-gen-cpmd-water`,
  `periodic-ri-mp2-si`. Incomplete remaining **0**; pending residual **9**;
  blocked **5**. Specs closed (deleted) for all 9. Roadmap: TS ✓, Periodic
  MP2 ✓ (Γ RI product; multi-k residual), MLWF ✓, 2nd-gen CP ✓; phonon
  IFC-Fourier productized; GPW BLYP / UKS-PBE / Fe smear / PW CELL_OPT prose.
  Clearup already ran (no isolation worktrees; no builds). **Queue not
  exhausted** — next wave GapCensus. No push.
- **2026-07-25 harness wave 1 (claude, earlier staff batch):** 5 merged into
  `dev` (head `bd39507a`) → done: `gpw-kpoints-si` (MP k-point GPW RKS),
  `gpw-variable-cell-ice` (vs CP2K CELL_OPT), `dfpt-phonon-si` (Γ DFPT, PW),
  `berry-polarization` (KSV, PW), `periodic-mp2-si` (Γ GPW MP2 **foothold** —
  later closed by `periodic-ri-mp2-si`). Incomplete then: `xl-bomd-water` WIP
  (later merged). Merge gate: cp2k domain suite green. Clearup: 71 merged
  `wf/*` pruned. Specs closed for all 5; roadmap Berry ✓, k-points prose;
  MP2 ○ at that time.
- **2026-07-25 bootstrap (parent `dev` @ `5bdae27d`):** GapCensus recover
  clean (incomplete=0, CP2K specs=0, wf/*=0, builds=none). Pending **14**
  re-verified true gaps; done **30**; blocked **5**. No roadmap ○/✓ flips.
  Ready for Assign. Suggested staff: `gpw-blyp-water`, `pw-uks-fe-smearing`,
  `pw-variable-cell-si`, `gpw-ts-ice`, `phonon-dispersion-si`,
  `periodic-ri-mp2-si`, `second-gen-cpmd-water`, `gpw-variable-cell-molopt`
  (hold revPBE/RPBE / MBD/vdW-DF / RPA/GW for surface order).
- **2026-07-25 GapCensus wave 1 (parent `dev` @ `f6838676`):** recover clean.
  Prior pending **14** re-verified; **+4** follow-on slugs (`gpw-uks-pbe`,
  `gpw-kpoints-uks`, `gpw-kpoint-forces`, `periodic-wannier-si`) → pending
  **18**; done **30**; blocked **5**. No table ○/✓ flips; fixed stale
  `docs/roadmap/index.md` GPW blurb. Ready for Assign.

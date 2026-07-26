# CP2K capability loop — handoff

## Harness wave 9 2026-07-26 (parent `dev` @ `dab2b190`)

| Field | Value |
| ----- | ----- |
| Branch | `dev` (tracks `origin/dev`, ahead of origin by merges) |
| Incomplete / `in_progress` | **2** (`[~] gpw-vv10-ice`, `[~] dfpt-spin-smearing`) |
| Open CP2K specs on `dev` | **0** (PySCF deferred only: `periodic-rhf-gamma`, `rsh-eri-cam-b3lyp`) |
| Specs closed this harness | **8** (deleted on close) |
| Unmerged staff `wf/*` | **2** incomplete kept: `wf/gpw-vv10-ice`, `wf/dfpt-spin-smearing` |
| Isolation worktrees | **0** (clearup already ran; main tree only) |
| Leftover builds | **none** |
| Pending gaps | **0** residual from prior census list (not exhausted — GapCensus owns refill) |
| Blocked | **5** (hybrid / GAPW / PAW / CUDA-GPW / meta-GGA) |
| Done this wave | **8** landings (see table) |
| Oracle | `tools/cp2k_oracle/` offline pins only |
| Loop docs | queue / loop / handoff / roadmap residual prose **synced** |

### Wave staff → done

| Slug | Domain | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) |
| ---- | ------ | ---------- | --------------------------------------------- |
| `ot-uks-smearing` | gpw | `reg_ot-uks-smearing` | **~0.12** |
| `ot-uks-dftu` | gpw | `reg_ot-uks-dftu` | **perf_exempt** |
| `ot-dftu-smearing` | gpw | `reg_ot-dftu-smearing` | **perf_exempt** |
| `gpw-d4-forces` | gpw | `reg_gpw-d4-forces` | **perf_exempt** |
| `gpw-mbd-forces` | gpw | `reg_gpw-mbd-forces` | **perf_exempt** |
| `gpw-uks-multigrid` | gpw | `reg_gpw-uks-multigrid` | **~1.54** |
| `dfpt-finite-q-smearing` | response | `reg_dfpt-finite-q-smearing` | **perf_exempt** |
| `second-gen-cpmd-ot` | dynamics | `reg_second-gen-cpmd-ot` | **perf_exempt** (NVE drift) |

### Incomplete kept `[~]` (do not demote)

| Slug | Branch | Tip | Note |
| ---- | ------ | --- | ---- |
| `gpw-vv10-ice` | `wf/gpw-vv10-ice` | `2f1d90e2` | CP2K 2026.1 has no stock `TYPE VV10` (≠ rVV10); incomplete |
| `dfpt-spin-smearing` | `wf/dfpt-spin-smearing` | `d7576dea` | feat+reg on branch; not merged to `dev` |

### Roadmap ○→✓ this harness

| Row | Doc | Note |
| --- | --- | ---- |
| *(none)* | — | D4/MBD/OT/multigrid table rows already ✓; stock VV10 stays ○ |
| Residual prose | `dft.md` §8/§10–12 + `index.md` + MD AIMD | closed UKS OT+smear; multi-k UKS OT+U; OT+U+smear; D4/MBD forces; UKS multigrid; finite-q DFPT smear; SGCP OT |
| Still open | DFPT spin+smear; stock VV10 | both `[~]` incomplete staff |

### Pending residual

**Pending = 0** residual from GapCensus wave 9 landed staff (emptied — **not**
queue exhausted). **Incomplete = 2** `[~]` (`gpw-vv10-ice`,
`dfpt-spin-smearing`). Next wave starts with **Recover incomplete** then
**GapCensus** (tree/roadmap ○ vs `dev`). Blocked = 5 unchanged.

Do **not** staff PySCF deferred specs. No push from parent.
**Queue not exhausted** — GapCensus owns next discovery.

---

## Harness wave 8 2026-07-26 (parent `dev` @ `a6752916`) — historical

| Field | Value |
| ----- | ----- |
| Branch | `dev` (tracks `origin/dev`, ahead of origin by merges) |
| Incomplete / `in_progress` | **1** (`[~] gpw-vv10-ice`) |
| Open CP2K specs on `dev` | **0** (PySCF deferred only: `periodic-rhf-gamma`, `rsh-eri-cam-b3lyp`) |
| Specs closed this harness | **9** (deleted on close) |
| Unmerged staff `wf/*` | **1** incomplete kept: `wf/gpw-vv10-ice` |
| Isolation worktrees | **0** (clearup already ran; main tree only) |
| Leftover builds | **none** |
| Pending gaps | **0** residual from prior census list (not exhausted — GapCensus owns refill) |
| Blocked | **5** (hybrid / GAPW / PAW / CUDA-GPW / meta-GGA) |
| Done this wave | **9** landings (see table) |
| Oracle | `tools/cp2k_oracle/` offline pins only |
| Loop docs | queue / loop / handoff / roadmap residual prose **synced** |

### Wave staff → done (historical)

| Slug | Domain | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) |
| ---- | ------ | ---------- | --------------------------------------------- |
| `periodic-urpa-kpoints` | gpw | `reg_periodic-urpa-kpoints` | **perf_exempt** |
| `ot-smearing` | gpw | `reg_ot-smearing` | **~0.20** |
| `ot-uks-kpoints` | gpw | `reg_ot-uks-kpoints` | **perf_exempt** |
| `ot-dftu-kpoints` | gpw | `reg_ot-dftu-kpoints` | **perf_exempt** |
| `gpw-d4-stress` | gpw | `reg_gpw-d4-stress` | **perf_exempt** |
| `gpw-vdw-scf-potential` | gpw | `reg_gpw-vdw-scf-potential` | **perf_exempt** |
| `dfpt-metallic-si` | response | `reg_dfpt-metallic-si` | **perf_exempt** |
| `dfpt-spin-incommensurate` | response | `reg_dfpt-spin-incommensurate` | **perf_exempt** |
| `gpw-optb88-vdw` | gpw | `reg_gpw-optb88-vdw` | **perf_exempt** (E_disp / XC combo) |

### Incomplete kept `[~]` (do not demote)

| Slug | Branch | Tip | Note |
| ---- | ------ | --- | ---- |
| `gpw-vv10-ice` | `wf/gpw-vv10-ice` | `d15cd804` | CP2K 2026.1 has no stock `TYPE VV10` (≠ rVV10); incomplete |

### Roadmap ○→✓ this harness (historical)

| Row | Doc | Note |
| --- | --- | ---- |
| optB88-vdW | `dft.md` §8 | exchange + DRSLL nonlocal (`reg_gpw-optb88-vdw`) |
| Residual prose | `dft.md` §8/§10–12 + `index.md` | multi-k URPA; OT+smear; multi-k UKS OT; multi-k OT+DFT+U; D4 stress; vdW SCF pot/forces; metallic DFPT; spin-incommensurate DFPT |
| Still open | stock VV10 | `[~] gpw-vv10-ice` (incomplete staff) |

### Pending residual

**Pending = 0** residual from GapCensus wave 8 staff list (emptied — **not**
queue exhausted). Then GapCensus wave 9 refilled 9 pending. **Incomplete = 1**
`[~]` (`gpw-vv10-ice`). Blocked = 5 unchanged.

Do **not** staff PySCF deferred specs. No push from parent.
**Queue not exhausted** — GapCensus owns next discovery.

---

## Harness wave 7 2026-07-26 (parent `dev` @ `8cb47664`) — historical

| Field | Value |
| ----- | ----- |
| Branch | `dev` (tracks `origin/dev`, ahead of origin by merges) |
| Incomplete / `in_progress` | **1** (`[~] gpw-vv10-ice`) |
| Open CP2K specs on `dev` | **0** (PySCF deferred only) |
| Specs closed this harness | **9** (deleted on close) |
| Pending residual | **0** then refilled by GapCensus wave 8 |
| Done this wave | **9** landings |

### Wave staff → done (historical)

| Slug | Domain | ratio |
| ---- | ------ | ----- |
| `periodic-urpa-si` | gpw | **perf_exempt** |
| `periodic-ugw-kpoints` | gpw | **perf_exempt** |
| `gpw-uks-smearing` | gpw | **~0.017** |
| `gpw-hirshfeld-uks` | gpw | **perf_exempt** (E_disp pin) |
| `cpmd-uks` | dynamics | **perf_exempt** (NVE drift) |
| `second-gen-cpmd-uks` | dynamics | **perf_exempt** (NVE drift) |
| `ot-kpoints` | gpw | **perf_exempt** |
| `ot-dftu` | gpw | **~0.26** |
| `dfpt-spin-si` | response | **perf_exempt** |

---

## Harness wave 6 2026-07-26 (parent `dev` @ `0869833a`) — historical

| Field | Value |
| ----- | ----- |
| Branch | `dev` (tracks `origin/dev`, ahead of origin by merges) |
| Incomplete / `in_progress` | **1** (`[~] gpw-vv10-ice`) |
| Open CP2K specs on `dev` | **0** (PySCF deferred only) |
| Specs closed this harness | **9** (deleted on close) |
| Pending residual | **1** (`periodic-urpa-si`) then refilled by GapCensus wave 7 |
| Done this wave | **9** landings |

### Wave staff → done (historical)

| Slug | Domain | ratio |
| ---- | ------ | ----- |
| `gpw-metallic-smearing` | gpw | **~0.011** |
| `gpw-xdm-ice` | gpw | **perf_exempt** |
| `gpw-uks-rpbe` | gpw | **~2.70** |
| `gpw-dftu-kpoint-stress` | gpw | **perf_exempt** |
| `periodic-ump2-kpoints` | gpw | **perf_exempt** |
| `periodic-ugw-si` | gpw | **perf_exempt** |
| `gpw-hirshfeld-kpoints` | gpw | **perf_exempt** |
| `dfpt-gga-si` | response | **perf_exempt** |
| `xl-bomd-uks` | dynamics | **perf_exempt** |

---

## Harness wave 5 2026-07-26 (parent `dev` @ `2d43d1c2`) — historical

| Field | Value |
| ----- | ----- |
| Branch | `dev` (tracks `origin/dev`, ahead of origin by merges) |
| Incomplete / `in_progress` | **2** (`[~] gpw-metallic-smearing`, `[~] gpw-vv10-ice`) |
| Open CP2K specs on `dev` | **0** (PySCF deferred only) |
| Specs closed this harness | **6** (deleted on close) |
| Pending residual | **1** (`gpw-xdm-ice`) then refilled by GapCensus wave 6 |
| Done this wave | **6** landings |

### Wave staff → done (historical)

| Slug | Domain | ratio |
| ---- | ------ | ----- |
| `gpw-dftu-kpoints` | gpw | **perf_exempt** |
| `periodic-gw-finite-q` | gpw | **perf_exempt** |
| `gpw-uks-revpbe` | gpw | **~3.22** |
| `cpmd-rattle-water` | dynamics | **perf_exempt** |
| `periodic-ump2-si` | gpw | **perf_exempt** |
| `ot-uks-o2` | gpw | **~8.50** |

---

## Harness wave 4 2026-07-26 (parent `dev` @ `85eaa2d9`) — historical

| Field | Value |
| ----- | ----- |
| Branch | `dev` (tracks `origin/dev`, ahead of origin by merges) |
| Incomplete / `in_progress` | **0** (no residual `[~]`) |
| Open CP2K specs | **0** (PySCF deferred only: `periodic-rhf-gamma`, `rsh-eri-cam-b3lyp`) |
| Specs closed this harness | **9** (deleted on close) |
| Unmerged staff `wf/*` | **0** for the 9 wave slugs (clearup already ran) |
| Isolation worktrees | **0** (main tree only) |
| Leftover builds | **none** |
| Pending gaps | **0** residual from prior census list (not exhausted — GapCensus owns refill) |
| Blocked | **5** (hybrid / GAPW / PAW / CUDA-GPW / meta-GGA) |
| Done this wave | **9** landings (see table) |
| Oracle | `tools/cp2k_oracle/` offline pins only |
| Loop docs | queue / loop / handoff / roadmap **synced** |

### Wave staff → done

| Slug | Domain | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) |
| ---- | ------ | ---------- | --------------------------------------------- |
| `gpw-kpoint-uks-stress` | gpw | `reg_gpw-kpoint-uks-stress` | **perf_exempt** |
| `periodic-gw-kpoints` | gpw | `reg_periodic-gw-kpoints` | **perf_exempt** |
| `gpw-vdw-df2-ice` | gpw | `reg_gpw-vdw-df2-ice` | **perf_exempt** (E_disp pin) |
| `gpw-rvv10-ice` | gpw | `reg_gpw-rvv10-ice` | **perf_exempt** (E_disp pin) |
| `gpw-mbd-periodic` | gpw | `reg_gpw-mbd-periodic` | **perf_exempt** (E_disp pin) |
| `pw-spin-gga` | pw | `reg_pw-spin-gga` | **perf_exempt** |
| `gpw-uks-blyp` | gpw | `reg_gpw-uks-blyp` | **~2.48** |
| `periodic-rpa-finite-q` | gpw | `reg_periodic-rpa-finite-q` | **perf_exempt** |
| `gpw-dftu-stress` | gpw | `reg_gpw-dftu-stress` | **perf_exempt** |

### Roadmap ○→✓ this harness

| Row | Doc | Note |
| --- | --- | ---- |
| vdW-DF2 | `dft.md` §8 | LMKLL nonlocal on GPW (`reg_gpw-vdw-df2-ice`) |
| rVV10 | `dft.md` §8 | nonlocal rVV10 on GPW (`reg_gpw-rvv10-ice`); VV10 stock still ○ |
| Multi-k G0W0 residual | `dft.md` §10 prose | closed (`reg_periodic-gw-kpoints`) |
| Periodic Ewald MBD | `dft.md` §8/§10 prose | productized (`reg_gpw-mbd-periodic`) |
| PW true collinear spin-GGA | `dft.md` §10 prose | productized (`reg_pw-spin-gga`) |
| Γ-UKS BLYP / multi-k UKS stress / finite-q dRPA / DFT+U stress | prose + index | productized |

### Pending residual

**Pending = 0** residual from GapCensus wave 4 staff list (emptied — **not**
queue exhausted). Next wave starts with **GapCensus** (tree/roadmap ○ vs
`dev`). Blocked = 5 unchanged.

Do **not** staff PySCF deferred specs. No push from parent.
**Queue not exhausted** — GapCensus owns next discovery.

---

## Harness wave 3 2026-07-25 (parent `dev` @ `0a642122`) — historical

| Field | Value |
| ----- | ----- |
| Branch | `dev` (tracks `origin/dev`, ahead of origin by merges) |
| Incomplete / `in_progress` | **0** (no residual `[~]`) |
| Open CP2K specs | **0** (PySCF deferred only: `periodic-rhf-gamma`, `rsh-eri-cam-b3lyp`) |
| Specs closed this harness | **8** (deleted on close) |
| Unmerged staff `wf/*` | **0** for the 8 wave slugs (clearup already ran) |
| Isolation worktrees | **0** (main tree only) |

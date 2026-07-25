# CP2K capability loop — handoff

## Harness wave 4 2026-07-26 (parent `dev` @ `85eaa2d9`)

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
| Leftover builds | **none** |
| Pending gaps | **0** residual from prior census list (not exhausted — GapCensus owns refill) |
| Blocked | **5** (hybrid / GAPW / PAW / CUDA-GPW / meta-GGA) |
| Done this wave | **8** landings (see table) |
| Oracle | `tools/cp2k_oracle/` offline pins only |
| Loop docs | queue / loop / handoff / roadmap **synced** |

### Wave staff → done

| Slug | Domain | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) |
| ---- | ------ | ---------- | --------------------------------------------- |
| `gpw-rpbe-water` | gpw | `reg_gpw-rpbe-water` | **perf_exempt** (dual-gate vs PBE) |
| `gpw-kpoint-uks-forces` | gpw | `reg_gpw-kpoint-uks-forces` | **perf_exempt** |
| `gpw-uks-stress` | gpw | `reg_gpw-uks-stress` | **perf_exempt** (Γ-UKS; multi-k UKS stress residual) |
| `gpw-vdw-df-ice` | gpw | `reg_gpw-vdw-df-ice` | **perf_exempt** (E_disp pin) |
| `gpw-mbd-hirshfeld` | gpw | `reg_gpw-mbd-hirshfeld` | **perf_exempt** (E_disp + V/V_free) |
| `periodic-rpa-kpoints` | gpw | `reg_periodic-rpa-kpoints` | **perf_exempt** |
| `periodic-gw-si` | gpw | `reg_periodic-gw-si` | **perf_exempt** |
| `dfpt-incommensurate-si` | response | `reg_dfpt-incommensurate-si` | **perf_exempt** |

### Roadmap ○→✓ this harness

| Row | Doc | Note |
| --- | --- | ---- |
| RPBE | `dft.md` §3 | `XCType`/`XcKind::RPBE`; GPW dual-gate |
| vdW-DF | `dft.md` §8 | DRSLL nonlocal on GPW (`reg_gpw-vdw-df-ice`) |
| Periodic GW | `dft.md` §10 | Γ spectral G0W0@PBE product; multi-k residual |
| Hirshfeld volumes | `dft.md` §8 prose | density-scaled TS/MBD (`reg_gpw-mbd-hirshfeld`) |
| Multi-k dRPA residual | `dft.md` §10 prose | closed (`reg_periodic-rpa-kpoints`) |
| Multi-k UKS forces / Γ-UKS stress / incommensurate DFPT | prose + index | productized; multi-k UKS stress residual |

### Pending residual (GapCensus wave 4 — 2026-07-26 @ `decf917f`)

**Pending = 9** were refilled then staffed and closed in harness wave 4 above.

Do **not** staff PySCF deferred specs. No push from parent.
**Queue not exhausted** — next is Assign from pending.

---

## Harness wave 2 2026-07-25 (parent `dev` @ `f5739506`) — historical

| Field | Value |
| ----- | ----- |
| Branch | `dev` (tracks `origin/dev`, ahead of origin by merges) |
| Incomplete / `in_progress` | **0** (no residual `[~]`) |
| Open CP2K specs | **0** (PySCF deferred only: `periodic-rhf-gamma`, `rsh-eri-cam-b3lyp`) |
| Specs closed this harness | **10** (deleted on close) |
| Unmerged staff `wf/*` | **0** for the 10 wave slugs (clearup already ran) |
| Isolation worktrees | **0** (main tree only) |
| Leftover builds | **none** |
| Pending gaps | **3** residual (not exhausted — GapCensus owns refill) |
| Blocked | **5** (hybrid / GAPW / PAW / CUDA-GPW / meta-GGA) |
| Done this wave | **10** landings (see table) |
| Oracle | `tools/cp2k_oracle/` offline pins only |
| Loop docs | queue / loop / handoff / roadmap **synced** |

### Wave staff → done

| Slug | Domain | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) |
| ---- | ------ | ---------- | --------------------------------------------- |
| `gpw-kpoint-forces` | gpw | `reg_gpw-kpoint-forces` | **perf_exempt** |
| `gpw-revpbe-si` | gpw | `reg_gpw-revpbe-si` | **perf_exempt** |
| `gpw-kpoints-uks` | gpw | `reg_gpw-kpoints-uks` | **perf_exempt** |
| `gpw-uks-forces` | gpw | `reg_gpw-uks-forces` | **~2.52** |
| `gpw-kpoint-stress` | gpw | `reg_gpw-kpoint-stress` | **perf_exempt** |
| `gpw-mbd-ice` | gpw | `reg_gpw-mbd-ice` | **perf_exempt** (libMBD E_disp) |
| `periodic-mp2-kpoints` | gpw | `reg_periodic-mp2-kpoints` | **perf_exempt** |
| `periodic-rpa-si` | gpw | `reg_periodic-rpa-si` | **perf_exempt** |
| `dfpt-finite-q-si` | response | `reg_dfpt-finite-q-si` | **perf_exempt** |
| `gpw-variable-cell-molopt` | gpw | `reg_gpw-variable-cell-molopt` | **~7.5** |

### Roadmap ○→✓ this harness

| Row | Doc | Note |
| --- | --- | ---- |
| revPBE | `dft.md` §3 | `XCType`/`XcKind::RevPBE`; GPW dual-gate |
| Many-body dispersion, MBD | `dft.md` §8 | free-atom MBD@rsSCS on GPW |
| Periodic RPA | `dft.md` §10 | Γ dRPA@KS product; multi-k / RI-RPA residual |
| Finite-q DFPT Sternheimer | `dft.md` §11 prose | PW X-point product path |
| Multi-k GPW MP2 residual | `dft.md` §10 prose | closed (`reg_periodic-mp2-kpoints`) |
| GPW k-UKS / k-forces / k-stress / Γ-UKS forces / MOLOPT CELL_OPT | prose + index | productized |

### Pending residual (next Assign / GapCensus)

1. `gpw-rpbe-water` (XcKind GGA expand — after revPBE ✓)
2. `gpw-vdw-df-ice` (after TS ✓ / MBD ✓ order)
3. `periodic-gw-si` (after RPA ✓)

Do **not** staff PySCF deferred specs. No push from parent.
**Queue not exhausted** — next wave starts with GapCensus.

---

## Harness wave 1 2026-07-25 (parent `dev` @ `83a100df`) — historical

| Field | Value |
| ----- | ----- |
| Branch | `dev` (tracks `origin/dev`, ahead of origin by merges) |
| Incomplete / `in_progress` | **0** (no residual `[~]`) |
| Open CP2K specs | **0** (PySCF deferred only: `periodic-rhf-gamma`, `rsh-eri-cam-b3lyp`) |
| Specs closed this harness | **9** (deleted on close) |
| Unmerged staff `wf/*` | **0** for the 9 wave slugs (clearup already ran) |
| Isolation worktrees | **0** (main tree only) |
| Leftover builds | **none** |
| Pending gaps | **9** residual (not exhausted — GapCensus owns refill) |
| Blocked | **5** (hybrid / GAPW / PAW / CUDA-GPW / meta-GGA) |
| Done this wave | **9** landings (see table) |
| Oracle | `tools/cp2k_oracle/` offline pins only |
| Loop docs | queue / loop / handoff / roadmap **synced** |

### Wave staff → done

| Slug | Domain | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) |
| ---- | ------ | ---------- | --------------------------------------------- |
| `gpw-blyp-water` | gpw | `reg_gpw-blyp-water` | **perf_exempt** (dual-gate vs PBE) |
| `gpw-uks-pbe` | gpw | `reg_gpw-uks-pbe` | **0.85** |
| `pw-uks-fe-smearing` | pw | `reg_pw-uks-fe-smearing` | **perf_exempt** |
| `pw-variable-cell-si` | pw | `reg_pw-variable-cell-si` | **perf_exempt** |
| `gpw-ts-ice` | gpw | `reg_gpw-ts-ice` | **perf_exempt** (libMBD E_disp) |
| `phonon-dispersion-si` | response | `reg_phonon-dispersion-si` | **perf_exempt** |
| `periodic-wannier-si` | response | `reg_periodic-wannier-si` | **perf_exempt** |
| `second-gen-cpmd-water` | dynamics | `reg_second-gen-cpmd-water` | **perf_exempt** |
| `periodic-ri-mp2-si` | gpw | `reg_periodic-ri-mp2-si` | **~0.09** |

### Roadmap ○→✓ this harness

| Row | Doc | Note |
| --- | --- | ---- |
| Tkatchenko–Scheffler | `dft.md` §8 | free-atom TS on GPW |
| Periodic MP2 | `dft.md` §10 | Γ RI-GPW product; multi-k residual |
| Maximally localized Wannier | `dft.md` §11 | molecular Boys + periodic PW MV |
| Second-generation Car–Parrinello MD | `molecular-dynamics.md` §AIMD | ASPC + partial SCF |
| Phonon ω(q) IFC Fourier | `dft.md` §11 prose | productized; finite-q DFPT still open |
| GPW BLYP / UKS-PBE / Fe smear / PW CELL_OPT | prose + index | already-✓ rows annotated |

### Pending residual (next Assign / GapCensus)

1. `gpw-variable-cell-molopt` (low priority re-twin)
2. `gpw-revpbe-si` / `gpw-rpbe-water` (XcKind GGA expand — one per wave)
3. `gpw-kpoints-uks` / `gpw-kpoint-forces`
4. `gpw-mbd-ice` / `gpw-vdw-df-ice` (after TS ✓)
5. `periodic-rpa-si` / `periodic-gw-si` (after MP2 ✓)

Do **not** staff PySCF deferred specs. No push from parent.
**Queue not exhausted** — next wave starts with GapCensus.

---

## Bootstrap 2026-07-25 (parent `dev` @ `5bdae27d`) — historical

| Field | Value |
| ----- | ----- |
| Incomplete / `in_progress` | **0** |
| Pending gaps | **14** (later GapCensus → 18) |
| Done | **30** landings prior to this harness wave |

## Parallel wave summary (2026-07-25) — historical

| Status | Slugs |
| ------ | ----- |
| merged → done | `mpi-pw-scaling-4rank`, `gpw-analytic-stress`, `wannier-centers`, `mpi-gpw-scaling-4rank`, `cpmd-water` (+ `xl-bomd-water`) |

## Claude wave 1 summary (harness 2026-07-25) — historical

| Status | Slugs |
| ------ | ----- |
| merged → done | `gpw-kpoints-si`, `gpw-variable-cell-ice`, `dfpt-phonon-si`, `berry-polarization`, `periodic-mp2-si` (Γ foothold; capability later closed by `periodic-ri-mp2-si`) |

## Oracle

Offline CP2K pins live in each regression `REFERENCE.md`.
Pin path: `tools/cp2k_oracle/` (see README). **Never** call live CP2K from ctest.

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

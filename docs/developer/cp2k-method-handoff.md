# CP2K capability loop — handoff

## Harness wave 1 2026-07-25 (parent `dev` @ `83a100df`)

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

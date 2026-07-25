# CP2K capability loop — handoff

## Bootstrap 2026-07-25 (parent `dev` @ `5bdae27d`)

| Field | Value |
| ----- | ----- |
| Branch | `dev` (tracks `origin/dev`) |
| Incomplete / `in_progress` | **0** |
| Open CP2K specs | **0** (PySCF deferred only: `periodic-rhf-gamma`, `rsh-eri-cam-b3lyp`) |
| Unmerged `wf/*` | **0** |
| Isolation worktrees | **0** (only main tree) |
| Leftover builds | **none** |
| Pending gaps | **14** (all re-verified true) |
| Blocked | **5** (hybrid / GAPW / PAW / CUDA-GPW / meta-GGA) |
| Done | **30** landings (waves 1–2 + claude wave 1 + parallel wave + xl-bomd) |
| Oracle | `tools/cp2k_oracle/README.md` + 23 inputs/notes; offline pins only |
| Loop docs | queue / loop / handoff / `oracle-capability-loop.md` **confirmed** |

**Next:** Assign wave from pending (incomplete ≫ pending is vacuous). Suggested
surface-disjoint staff (≤10):

| Priority | Slug | Why now | Surface caution |
| -------- | ---- | ------- | --------------- |
| 1 | `gpw-blyp-water` | §3 BLYP molecular ✓ / GPW `XcKind` still 3-valued | only **one** XcKind GGA expand this wave |
| 2 | `pw-uks-fe-smearing` | follow-on after O₂ k-mesh UKS ✓ | `cpu/planewave` UKS+smear |
| 3 | `pw-variable-cell-si` | dep Nielsen–Martin stress ✓ | minimize + PW stress |
| 4 | `gpw-ts-ice` | §8 TS ○; D3/D4 pin patterns ✓ | new dispersion path |
| 5 | `phonon-dispersion-si` | dep FD + Γ DFPT ✓ | DFPT q≠0 |
| 6 | `periodic-ri-mp2-si` | foothold ≠ closed (Periodic MP2 ○) | GPW MP2 productization |
| 7 | `second-gen-cpmd-water` | dep `cpmd-water` ✓ | dynamics / Cpmd |
| 8 | `gpw-variable-cell-molopt` | declared re-twin; low risk | mostly pin/reg |
| hold | `gpw-revpbe-si` / `gpw-rpbe-water` | after BLYP XcKind lands | same `XcKind` files |
| hold | `gpw-mbd-ice` / `gpw-vdw-df-ice` | after TS foothold | dispersion stack |
| hold | `periodic-rpa-si` / `periodic-gw-si` | after RI-MP2 or staff-heavy | post-HF stack |

Do **not** staff PySCF deferred specs. No push from parent.

---

Last land: **2026-07-25 parallel wave** — 5 staffed `[~]` merged to `dev`:
`mpi-pw-scaling-4rank`, `gpw-analytic-stress`, `wannier-centers`,
`mpi-gpw-scaling-4rank`, `cpmd-water` (+ earlier `xl-bomd-water`).
**in_progress empty.** Pending residual for next census/assign.

## Parallel wave summary (2026-07-25)

| Slug | Domain | Status | Regression | Note |
| ---- | ------ | ------ | ---------- | ---- |
| `mpi-pw-scaling-4rank` | pw | **merged → done** | `reg_mpi-pw-scaling-4rank` | k-RR 4 ranks; speedup 2.54× |
| `gpw-analytic-stress` | gpw | **merged → done** | `reg_gpw-analytic-stress` | fixed-P+Pulay vs FD ~3e-10 |
| `wannier-centers` | scf | **merged → done** | `reg_wannier-centers` | Foster–Boys molecular foothold |
| `mpi-gpw-scaling-4rank` | gpw | **merged → done** | `reg_mpi-gpw-scaling-4rank` | k-RR; energy rank-identical |
| `cpmd-water` | dynamics | **merged → done** | `reg_cpmd-water` | CP E_cons drift 2e-3 |

Specs deleted on close. Incomplete remaining = **0**.

## xl-bomd-water land (2026-07-25)

| Field | Value |
| ----- | ----- |
| Status | **merged → done** |
| API | `cpu::gpw::XlBomd` (Niklasson JCP 130, 214109 Table I); `GpwDft::run(..., guess_P)`; `GPWEvaluator::{xl_bomd, xl_bomd_k}` |
| Unit | `test_cpu_gpw_xl_bomd` **16/16 PASS** |
| Regression | `reg_xl-bomd-water` **PASS** — worst \|ΔE\|/scale = **1.2e-5** (tol 0.05), wall **33 s**, SCF iters 6 |
| A9 guards | `test_cpu_gpw_bomd` / `solver` / `forces` **PASS** |
| perf | **perf_exempt** (no CP2K binary on host; twin input under `tools/cp2k_oracle/inputs/xl_bomd_water_nve.inp`) |
| Roadmap | MD §AIMD **XL-BOMD → ✓** (UKS/OT/k-point XL still ○) |
| Spec | closed / deleted |

**Still `[~]` (staff wave 2):** `wannier-centers`, `cpmd-water`,
`mpi-gpw-scaling-4rank`, `mpi-pw-scaling-4rank`, `gpw-analytic-stress`.

## Claude wave 1 summary (harness 2026-07-25)

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) |
| ---- | ------ | ------ | ---------- | --------------------------------------------- |
| `gpw-kpoints-si` | gpw | **merged → done** | `regressions/gpw/gpw-kpoints-si` | **perf_exempt** |
| `gpw-variable-cell-ice` | gpw | **merged → done** | `regressions/gpw/gpw-variable-cell-ice` | walls ATV 1.9 s / CP2K 3.9 s (GTH chain) |
| `dfpt-phonon-si` | response | **merged → done** | `regressions/response/dfpt-phonon-si` | **perf_exempt** |
| `berry-polarization` | response | **merged → done** | `regressions/pw/berry-polarization` | **perf_exempt** |
| `periodic-mp2-si` | gpw | **merged → done** (Γ **foothold**; capability row stays ○) | `regressions/gpw/periodic-mp2-si` | **perf_exempt** |
| `xl-bomd-water` | dynamics | **merged → done** (2026-07-25 resume) | `regressions/dynamics/xl-bomd-water` | **perf_exempt** |

- **Not staffed wave 1, still `[~]`:** `wannier-centers`, `cpmd-water`,
  `mpi-gpw-scaling-4rank`, `mpi-pw-scaling-4rank` (+ wave-2 promote
  `gpw-analytic-stress`).
- **Merge gate:** cp2k domain suite green on `dev`
  (`ctest -R '^test_cpu_(gpw|pw|scf|xc)_|^test_xc_'`). Dev head `bd39507a`.
- **Specs closed:** all 5 wave-1 slugs + `xl-bomd-water` deleted on close.
- **Clearup:** 71 fully-merged `wf/*` tips pruned; isolation worktrees removed.
- **Roadmap sync:** dft.md — GPW k-point RKS prose ✓, Berry-phase
  polarization ✓, GPW variable-cell productization noted; **Periodic MP2 row
  kept ○** (foothold annotated); DFPT row was already Γ-scoped ✓ by the
  worker merge. optimization.md — GPW CELL_OPT productization note updated.
  molecular-dynamics.md — **XL-BOMD ✓**.

## Assign wave 3 (2026-07-24) — staff 10

| Priority | Slug | Chapter | is_resume | Notes |
| -------- | ---- | ------- | --------- | ----- |
| 1 | `gpw-kpoints-si` | gpw | **true** | only incomplete; no open spec — write-first |
| 2 | `gpw-variable-cell-ice` | gpw | false | dep stress ✓ |
| 3 | `dfpt-phonon-si` | response | false | dep FD phonons ✓ |
| 4 | `xl-bomd-water` | dynamics | false | dep GPW BOMD ✓ |
| 5 | `berry-polarization` | response | false | dep bands ✓ |
| 6 | `wannier-centers` | response | false | dep bands ✓ |
| 7 | `cpmd-water` | dynamics | false | after solid BOMD ✓ |
| 8 | `periodic-mp2-si` | gpw | false | solid Γ GPW energy ✓ |
| 9 | `mpi-gpw-scaling-4rank` | infra | false | industrial; `perf_exempt` OK |
| 10 | `mpi-pw-scaling-4rank` | infra | false | k-point RR; industrial |

- **Held pending (2):** `periodic-rpa-si` (after MP2), `periodic-gw-si` (staff full)
- **Queue after assign:** pending = 2; in_progress = 10; done = 19; blocked = 5
- **Workers:** each NEW slug writes `.claude/specs/<slug>/` first; resume
  `gpw-kpoints-si` also needs a live spec (INDEX: no open CP2K spec dir)
- **Recover 2026-07-24:** all 10 still `[~]` with no branch/spec/regression
  progress; `resume_from=write-spec-first` for the full staff set. Do not
  start brand-new pending (`periodic-rpa-si` / `periodic-gw-si`) until these
  finish or staff slots free.

## Wave 2 summary (harness 2026-07-24)

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) |
| ---- | ------ | ------ | ---------- | --------------------------------------------- |
| `gpw-stress-water` | gpw | **merged → done** | `regressions/gpw/gpw-stress-water` | **4.75** |
| `gpw-bomd-water` | dynamics | **merged → done** | `regressions/gpw/gpw-bomd-water` | **perf_exempt** |
| `pw-analytic-stress-si` | pw | **merged → done** | `regressions/pw/pw-analytic-stress-si` | **perf_exempt** |
| `fd-phonon-supercell` | response | **merged → done** | `regressions/pw/fd-phonon-supercell` | **perf_exempt** |
| `gpw-dftu-nio` | gpw | **merged → done** | `regressions/gpw/gpw-dftu-nio` | **0.21** |
| `bands-dos-si` | pw | **merged → done** | `regressions/pw/bands-dos-si` | **perf_exempt** |
| `pw-uks-kmesh-fe` | pw | **merged → done** | `regressions/pw/pw-uks-kmesh-fe` | **perf_exempt** |
| `gpw-pbesol-si` | gpw | **merged → done** | `regressions/gpw/gpw-pbesol-si` | **perf_exempt** |
| `gpw-d4-water` | gpw | **merged → done** | `regressions/gpw/gpw-d4-water` | **perf_exempt** |
| `gpw-kpoints-si` | gpw | **→ wave 3 resume** | — | incomplete leftover from assign wave 2 |

- **Blocked this wave:** none new (hybrid/GAPW/PAW/CUDA-GPW/meta-GGA stand)
- **Merge-failed this wave:** none
- **Specs closed:** all 9 wave-2 merged slugs (acceptance met on `dev`);
  do not re-open without superseding intent.

## Wave 1 summary (historical)

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) |
| ---- | ------ | ------ | ---------- | --------------------------------------------- |
| `gpw-water-pbe-mt` | gpw | **merged → done** | `regressions/gpw/gpw-water-pbe-mt` | **perf_exempt** |
| `gpw-forces-water` | gpw | **merged → done** | `regressions/gpw/gpw-forces-water` | **3.56** |
| `gpw-uks-o2-gamma` | gpw | **merged → done** | `regressions/gpw/gpw-uks-o2-gamma` | **0.98** |
| `gpw-pbe-d3bj-ice` | gpw | **merged → done** | `regressions/gpw/gpw-pbe-d3bj-ice` | **perf_exempt** |
| `pw-lda-si` | pw | **merged → done** | `regressions/pw/pw-lda-si` | **perf_exempt** |
| `pw-pbe-gga-si` | pw | **merged → done** | `regressions/pw/pw-pbe-gga-si` | **perf_exempt** |
| `pw-smearing-al` | pw | **merged → done** | `regressions/pw/pw-smearing-al` | **0.008** |
| `pw-spin-polarized-o2` | pw | **merged → done** | `regressions/pw/pw-spin-polarized-o2` | **perf_exempt** |
| `pw-bomd-h` | dynamics | **merged → done** | `regressions/pw/pw-bomd-h` | **perf_exempt** |
| `ot-minimizer-water` | gpw | **merged → done** | `regressions/gpw/ot-minimizer-water` | **~9.0** |

## Next after wave 3 (held pending)

1. `periodic-rpa-si` (**after** `periodic-mp2-si` foothold)
2. `periodic-gw-si` (deps clear; staff-held)

## Oracle

Offline CP2K pins live in each regression `REFERENCE.md`.
Pin path: `tools/cp2k_oracle/` (see README). **Never** call live CP2K from ctest.

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

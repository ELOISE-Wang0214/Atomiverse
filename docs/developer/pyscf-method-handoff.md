# PySCF method loop — handoff

Last harness pass: **2026-07-23 wave 6** on `dev` (harness commit series `wf wave 3`).

## Wave 6 summary

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{pyscf}}\) |
| ---- | ------ | ------ | ---------- | ----------------------------------------------- |
| `ccsd-t-h2o-sto3g` | post-hf | **merged → done** | `regressions/scf/ccsd-t-h2o-sto3g` | **0.333** |
| `casci-h2-sto3g` | post-hf | **merged → done** | `regressions/scf/casci-h2-sto3g` | **0.076** |
| `delta-scf-n2` | excited-states / hf | **merged → done** | `regressions/scf/delta-scf-n2` | **7.83** |
| `eom-ee-ccsd` | excited-states | **merged → done** | `regressions/scf/eom-ee-ccsd` | **0.193** |

- **Blocked this wave:** none
- **Merge-failed this wave:** none
- **Queue:** **pending empty** (PySCF-oracle census queue exhausted for this
  shift); in_progress empty; historical blocked list (cam-b3lyp,
  periodic-rhf-gamma, td-dftb, gfn2-xtb, pm6-pm7) unchanged — not re-opened.
- **Roadmap:** CPU ✓ for **CCSD(T)** (`RHF::ccsd_t`), **CASCI**
  (`RHF::casci`), **EOM-EE-CCSD** (`RHF::eom_ee_ccsd`). **ΔSCF / MOM / IMOM**
  table rows were already ✓ (`UHF::delta_scf`). **No ⚡** — ratios noted only;
  no standing perf suite in ctest.
- **Speed notes:** all four fair-benched and within 10×. casci / eom-ee /
  ccsd-t faster than PySCF on these footholds (0.08–0.33×); delta-scf 7.83×
  (IMOM re-SCF heavier but under the 10× gate — **not** blocked-for-speed).
  No blocked-for-speed slugs.

## Wave 5 (prior)

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{pyscf}}\) |
| ---- | ------ | ------ | ---------- | ----------------------------------------------- |
| `adc2-h2co-sto3g` | excited-states | **merged → done** | `regressions/scf/adc2-h2co-sto3g` | **0.267** |
| `mp3-h2-sto3g` | post-hf | **merged → done** | `regressions/scf/mp3-h2-sto3g` | **0.774** |
| `ccsd-h2-sto3g` | post-hf | **merged → done** | `regressions/scf/ccsd-h2-sto3g` | **0.317** |
| `g0w0-h2-sto3g` | excited-states | **merged → done** | `regressions/scf/g0w0-h2-sto3g` | **0.109** |
| `b2plyp-h2o` | dft | **merged → done** | `regressions/scf/b2plyp-h2o` | **1.75** |

## Wave 4 (prior)

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{pyscf}}\) |
| ---- | ------ | ------ | ---------- | ----------------------------------------------- |
| `romp2-o2-sto3g` | post-hf | **merged → done** | `regressions/scf/romp2-o2-sto3g` | **1.37** |
| `roks-grad-o2-sto3g` | dft | **merged → done** | `regressions/scf/roks-grad-o2-sto3g` | **1.117** |
| `h2-b3lyp-rpa-tddft` | excited-states | **merged → done** | `regressions/scf/h2-b3lyp-rpa-tddft` | **0.430** |
| `spin-flip-tda-o2` | excited-states | **merged → done** | `regressions/scf/spin-flip-tda-o2` | **9.56** |
| `rohf-hessian-o2` | hartree-fock | **merged → done** | `regressions/scf/rohf-hessian-o2` | **perf_exempt** (no PySCF `hessian.ROHF`) |

## Wave 3 (prior)

| Slug | Domain | Status | Regression |
| ---- | ------ | ------ | ---------- |
| `sos-mp2-h2o-sto3g` | post-hf | **merged → done** | `regressions/scf/sos-mp2-h2o-sto3g` |
| `mp2-grad-h2o-sto3g` | post-hf | **merged → done** | `regressions/scf/mp2-grad-h2o-sto3g` |
| `uks-tda-o2-lsda` | excited-states | **merged → done** | `regressions/scf/uks-tda-o2-lsda` |
| `d3bj-molecular-pbe` | dft | **merged → done** | `regressions/scf/d3bj-molecular-pbe` |
| `uhf-hessian-o2` | hartree-fock | **merged → done** | `regressions/scf/uhf-hessian-o2` |

## Wave 2 (prior)

| Slug | Domain | Status | Regression |
| ---- | ------ | ------ | ---------- |
| `scs-mp2-h2o-sto3g` | post-hf | **merged → done** | `regressions/scf/scs-mp2-h2o-sto3g` |
| `full-lr-tddft-h2-pbe` | excited-states | **merged → done** | `regressions/scf/full-lr-tddft-h2-pbe` |
| `roks-o2-sto3g` | dft | **merged → done** | `regressions/scf/roks-o2-sto3g` |
| `static-polariz-h2o` | hartree-fock | **merged → done** | `regressions/scf/static-polariz-h2o` |
| `ief-pcm-h2o` | multiscale | **merged → done** | `regressions/scf/ief-pcm-h2o` |

## Wave 1 (prior)

| Slug | Domain | Status | Regression |
| ---- | ------ | ------ | ---------- |
| `frozen-core-mp2-h2o` | post-hf | **merged → done** | `regressions/scf/frozen-core-mp2-h2o` |
| `cisd-h2co-sto3g` | excited-states | **merged → done** | `regressions/scf/cisd-h2co-sto3g` |
| `blyp-water-sto3g` | dft | **merged → done** | `regressions/scf/blyp-water-sto3g` |
| `ghf-grad-h3-sto3g` | hartree-fock | **merged → done** | `regressions/scf/ghf-grad-h3-sto3g` |
| `pcm-grad-h2o` | multiscale | **merged → done** | `regressions/scf/pcm-grad-h2o` |

## Oracle

Pins live in each regression `REFERENCE.md`. Offline pin path:
`tools/setup_pyscf_env.sh` + `tools/pyscf_oracle/pin_scf.py` (union restored
after serial merges in prior `fix(oracle)` / `fix(wf): re-union…` commits;
wave-6 workers each carried their pin path).

## Next assign candidates (pending head)

**None** — `docs/developer/pyscf-method-queue.md` **pending** is empty.
Historical **blocked** rows stay blocked (erf ERI, periodic HF, TD-DFTB, GFN2,
PM6/7). Further census would need new roadmap ○ rows with a PySCF offline
oracle twin (e.g. CASSCF, EOM-IP/EA, CC gradients) before a new assign wave.

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

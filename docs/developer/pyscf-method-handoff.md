# PySCF method loop — handoff

Last harness pass: **2026-07-23 wave 5** on `dev` (harness commit series `wf wave 2`).

## Wave 5 summary

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{pyscf}}\) |
| ---- | ------ | ------ | ---------- | ----------------------------------------------- |
| `adc2-h2co-sto3g` | excited-states | **merged → done** | `regressions/scf/adc2-h2co-sto3g` | **0.267** |
| `mp3-h2-sto3g` | post-hf | **merged → done** | `regressions/scf/mp3-h2-sto3g` | **0.774** |
| `ccsd-h2-sto3g` | post-hf | **merged → done** | `regressions/scf/ccsd-h2-sto3g` | **0.317** |
| `g0w0-h2-sto3g` | excited-states | **merged → done** | `regressions/scf/g0w0-h2-sto3g` | **0.109** |
| `b2plyp-h2o` | dft | **merged → done** | `regressions/scf/b2plyp-h2o` | **1.75** |

- **Blocked this wave:** none
- **Merge-failed this wave:** none
- **Queue:** 4 pending (`ccsd-t-h2o-sto3g`, `casci-h2-sto3g`, `delta-scf-n2`,
  `eom-ee-ccsd` — unblocked now that CCSD energy ships); in_progress empty;
  historical blocked list (cam-b3lyp, periodic-rhf-gamma, td-dftb, gfn2-xtb,
  pm6-pm7) unchanged
- **Roadmap:** CPU ✓ for **MP3** (`RHF::mp3`), **CCSD** (`RHF::ccsd`),
  **G0W0** (`RHF::g0w0`). **ADC(2)** and **B2PLYP** table rows were already ✓;
  B2PLYP prose updated (composed hybrid KS + scaled MP2). EOM-CC blocked note
  refreshed (CCSD energy yes; EOM layer still ○). **No ⚡** — ratios noted only;
  no standing perf suite in ctest.
- **Speed notes:** all five fair-benched and within 10×; g0w0 / adc2 / ccsd
  faster than PySCF on these footholds (0.11–0.32×); mp3 0.77×; b2plyp 1.75×.
  No blocked-for-speed slugs.

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
after serial merges in `fix(oracle): re-union pin_scf/fair_bench…` and
`fix(wf): re-union ccsd/mp3/b2plyp…` / `re-union adc2 after g0w0…`).

## Next assign candidates (pending head)

1. `ccsd-t-h2o-sto3g` — (T) triples on closed-shell CCSD
2. `casci-h2-sto3g` — active-space on FCI kernel
3. `delta-scf-n2` — MOM/IMOM foothold
4. `eom-ee-ccsd` — EOM-EE on the new CCSD energy path

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

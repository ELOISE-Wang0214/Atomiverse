# PySCF method loop — handoff

Last harness pass: **2026-07-24 wave 2 (post re-census)** on `dev`
(`chore(harness): wf wave 2`).

## Wave 2 (post re-census) summary

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{pyscf}}\) | Tutorial |
| ---- | ------ | ------ | ---------- | ----------------------------------------------- | -------- |
| `eom-cc2-h2` | excited-states | **merged → done** | `regressions/scf/eom-cc2-h2` | **0.300** | [methods/eom-cc2-h2](../tutorial/methods/eom-cc2-h2.md) |
| `eom-ip-ccsd` | excited-states | **merged → done** | `regressions/scf/eom-ip-ccsd` | **0.754** | [methods/eom-ip-ccsd](../tutorial/methods/eom-ip-ccsd.md) |
| `uks-full-lr-o2-lsda` | excited-states | **merged → done** | `regressions/scf/uks-full-lr-o2-lsda` | **2.937** | [methods/uks-full-lr-o2-lsda](../tutorial/methods/uks-full-lr-o2-lsda.md) |
| `casscf-h2-sto3g` | post-hf | **merged → done** | `regressions/scf/casscf-h2-sto3g` | **0.0435** | [methods/casscf-h2-sto3g](../tutorial/methods/casscf-h2-sto3g.md) |
| `ccsd-grad-h2o-sto3g` | post-hf | **merged → done** | `regressions/scf/ccsd-grad-h2o-sto3g` | **1.665** | [methods/ccsd-grad-h2o-sto3g](../tutorial/methods/ccsd-grad-h2o-sto3g.md) |
| `rks-hessian-h2o-sto3g` | dft | **merged → done** | `regressions/scf/rks-hessian-h2o-sto3g` | **3.08** (vs FD twin) | [methods/rks-hessian-h2o-sto3g](../tutorial/methods/rks-hessian-h2o-sto3g.md) |
| `ump2-grad-o2-sto3g` | post-hf | **merged → done** | `regressions/scf/ump2-grad-o2-sto3g` | **1.37** | [methods/ump2-grad-o2-sto3g](../tutorial/methods/ump2-grad-o2-sto3g.md) |

- **Blocked this wave:** none
- **Merge-failed this wave:** none
- **Not in this merge set (still in_progress):** `uks-hessian-o2-sto3g`
- **Queue:** pending = 1 (`nevpt2-h2-sto3g` unblocked by CASSCF ✓); in_progress =
  1 retained (`uks-hessian-o2-sto3g`); historical blocked list unchanged
  (cam-b3lyp, periodic-rhf-gamma, td-dftb, gfn2-xtb, pm6-pm7, caspt2).
- **Roadmap:** CPU ✓ for **EOM-IP-CCSD**, **UKS full-LR TDDFT**, **CASSCF**,
  **CCSD nuclear gradient** (`RHF::ccsd_forces`; FD of CCSD energy foothold),
  **RKS nuclear Hessian** (`RKS::hessian`; FD of forces), **UMP2 analytic
  gradient**. **EOM-CC2** was already ✓ from the merge commit prose before
  harness. **No ⚡** — ratios noted only; no standing perf suite in ctest.
- **Speed notes:** all seven fair-benched and within 10×. Fastest foothold
  casscf 0.0435×; eom-cc2 0.300×; eom-ip 0.754×; ump2-grad 1.37×; ccsd-grad
  1.665×; uks-full-lr 2.937×; rks-hessian 3.08× vs FD twin (unfair vs PySCF
  analytic CPKS ~33× — not the product gate). No blocked-for-speed slugs.

## Wave 1 (post re-census)

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{pyscf}}\) |
| ---- | ------ | ------ | ---------- | ----------------------------------------------- |
| `eom-ea-ccsd` | excited-states | **merged → done** | `regressions/scf/eom-ea-ccsd` | **0.068** |
| `adc2x-h2co-sto3g` | excited-states | **merged → done** | `regressions/scf/adc2x-h2co-sto3g` | **0.240** |
| `ip-adc2-h2o-sto3g` | excited-states | **merged → done** | `regressions/scf/ip-adc2-h2o-sto3g` | **0.191** |
| `ea-adc2-h2-sto3g` | excited-states | **merged → done** | `regressions/scf/ea-adc2-h2-sto3g` | **0.213** |
| `sf-tddft-o2` | excited-states | **merged → done** | `regressions/scf/sf-tddft-o2` | **9.51** |
| `mp4-h2-sto3g` | post-hf | **merged → done** | `regressions/scf/mp4-h2-sto3g` | **0.031** |
| `cc2-h2-sto3g` | post-hf | **merged → done** | `regressions/scf/cc2-h2-sto3g` | **0.392** |
| `direct-rpa-h2-sto3g` | excited-states | **merged → done** | `regressions/scf/direct-rpa-h2-sto3g` | **0.760** |
| `bse-h2-sto3g` | excited-states | **merged → done** | `regressions/scf/bse-h2-sto3g` | **0.060** |
| `dft-d4-h2o-pbe` | dft | **merged → done** | `regressions/scf/dft-d4-h2o-pbe` | **2.17** |
| `eom-sf-ccsd-o2` | excited-states | **merged → done** | `regressions/scf/eom-sf-ccsd-o2` | **3.12** |
| `gf2-h2-sto3g` | excited-states | **merged → done** | `regressions/scf/gf2-h2-sto3g` | **0.071** |

## Wave 6 (prior, pre re-census)

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{pyscf}}\) |
| ---- | ------ | ------ | ---------- | ----------------------------------------------- |
| `ccsd-t-h2o-sto3g` | post-hf | **merged → done** | `regressions/scf/ccsd-t-h2o-sto3g` | **0.333** |
| `casci-h2-sto3g` | post-hf | **merged → done** | `regressions/scf/casci-h2-sto3g` | **0.076** |
| `delta-scf-n2` | excited-states / hf | **merged → done** | `regressions/scf/delta-scf-n2` | **7.83** |
| `eom-ee-ccsd` | excited-states | **merged → done** | `regressions/scf/eom-ee-ccsd` | **0.193** |

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

## Wave 2 (prior, first exhaust)

| Slug | Domain | Status | Regression |
| ---- | ------ | ------ | ---------- |
| `scs-mp2-h2o-sto3g` | post-hf | **merged → done** | `regressions/scf/scs-mp2-h2o-sto3g` |
| `full-lr-tddft-h2-pbe` | excited-states | **merged → done** | `regressions/scf/full-lr-tddft-h2-pbe` |
| `roks-o2-sto3g` | dft | **merged → done** | `regressions/scf/roks-o2-sto3g` |
| `static-polariz-h2o` | hartree-fock | **merged → done** | `regressions/scf/static-polariz-h2o` |
| `ief-pcm-h2o` | multiscale | **merged → done** | `regressions/scf/ief-pcm-h2o` |

## Wave 1 (prior, first exhaust)

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
after serial merges in prior `fix(oracle)` / `fix(wf): re-union…` /
`fix(merge): complete CC2+UCCSD…` commits). Wave-2 workers each carried their
pin path; merge commits landed on `dev` before this harness pass.

## Next assign candidates (pending head)

1. `nevpt2-h2-sto3g` — unblocked (CASSCF energy ships)
2. Re-dispatch `uks-hessian-o2-sto3g` if still mid-flight / not merged

Historical **blocked** rows stay blocked (erf ERI, periodic HF, TD-DFTB, GFN2,
PM6/7, CASPT2).

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

# PySCF method loop — handoff

Last harness pass: **2026-07-23 wave 1 (post re-census)** on `dev`
(`chore(harness): wf wave 1`).

## Wave 1 (post re-census) summary

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{pyscf}}\) | Tutorial |
| ---- | ------ | ------ | ---------- | ----------------------------------------------- | -------- |
| `eom-ea-ccsd` | excited-states | **merged → done** | `regressions/scf/eom-ea-ccsd` | **0.068** | `methods/eom-ea-ccsd.md` |
| `adc2x-h2co-sto3g` | excited-states | **merged → done** | `regressions/scf/adc2x-h2co-sto3g` | **0.240** | `methods/adc2x-h2co-sto3g.md` |
| `ip-adc2-h2o-sto3g` | excited-states | **merged → done** | `regressions/scf/ip-adc2-h2o-sto3g` | **0.191** | `methods/ip-adc2-h2o-sto3g.md` |
| `ea-adc2-h2-sto3g` | excited-states | **merged → done** | `regressions/scf/ea-adc2-h2-sto3g` | **0.213** | `methods/ea-adc2-h2-sto3g.md` |
| `sf-tddft-o2` | excited-states | **merged → done** | `regressions/scf/sf-tddft-o2` | **9.51** | `methods/sf-tddft-o2.md` |
| `mp4-h2-sto3g` | post-hf | **merged → done** | `regressions/scf/mp4-h2-sto3g` | **0.031** | `methods/mp4-h2-sto3g.md` |
| `cc2-h2-sto3g` | post-hf | **merged → done** | `regressions/scf/cc2-h2-sto3g` | **0.392** | `methods/cc2-h2-sto3g.md` |
| `direct-rpa-h2-sto3g` | excited-states | **merged → done** | `regressions/scf/direct-rpa-h2-sto3g` | **0.760** | `methods/direct-rpa-h2-sto3g.md` |
| `bse-h2-sto3g` | excited-states | **merged → done** | `regressions/scf/bse-h2-sto3g` | **0.060** | `methods/bse-h2-sto3g.md` |
| `dft-d4-h2o-pbe` | dft | **merged → done** | `regressions/scf/dft-d4-h2o-pbe` | **2.17** | `methods/dft-d4-h2o-pbe.md` |
| `eom-sf-ccsd-o2` | excited-states | **merged → done** | `regressions/scf/eom-sf-ccsd-o2` | **3.12** | `methods/eom-sf-ccsd-o2.md` |
| `gf2-h2-sto3g` | excited-states | **merged → done** | `regressions/scf/gf2-h2-sto3g` | **0.071** | `methods/gf2-h2-sto3g.md` |

- **Blocked this wave:** none
- **Merge-failed this wave:** none
- **Not in this merge set (still in_progress):** eom-ip-ccsd, uks-full-lr-o2-lsda,
  casscf-h2-sto3g, ccsd-grad-h2o-sto3g, rks-hessian-h2o-sto3g,
  ump2-grad-o2-sto3g, uks-hessian-o2-sto3g
- **Queue:** pending = 2 (`nevpt2-h2-sto3g` needs CASSCF; `eom-cc2-h2` unblocked
  by CC2 energy ✓); in_progress = 7 retained; historical blocked list unchanged
  (cam-b3lyp, periodic-rhf-gamma, td-dftb, gfn2-xtb, pm6-pm7, caspt2).
- **Roadmap:** CPU ✓ for **EOM-EA-CCSD**, **EOM-SF-CCSD**, **ADC(2)-x**,
  **IP-ADC(2)**, **EA-ADC(2)**, **Spin-flip TDDFT**, **MP4(SDTQ)**, **CC2**,
  **Direct RPA**, **BSE**, **GF2**, **DFT-D4**. **No ⚡** — ratios noted only;
  no standing perf suite in ctest.
- **Speed notes:** all twelve fair-benched and within 10×. Fastest footholds
  mp4/bse/eom-ea/gf2 (0.03–0.07×); sf-tddft 9.51× (dense SF-RPA vs tiny twin —
  under the 10× gate, **not** blocked-for-speed); eom-sf 3.12×; dft-d4 2.17×
  (2000-eval loop). No blocked-for-speed slugs.

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

## Wave 2 (prior)

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
`fix(merge): complete CC2+UCCSD…` commits). Wave-1 workers each carried their
pin path; merge-fix commits re-unioned multi-channel ADC / CC2+UCCSD.

## Next assign candidates (pending head)

1. `nevpt2-h2-sto3g` — still gated on CASSCF (`casscf-h2-sto3g` in_progress)

**Wave 2 assigned (in_progress, 8):** eom-cc2-h2 (new from pending; CC2 ✓),
plus re-dispatch of wave-1 retain: eom-ip-ccsd, uks-full-lr-o2-lsda,
casscf-h2-sto3g, ccsd-grad-h2o-sto3g, rks-hessian-h2o-sto3g, ump2-grad-o2-sto3g,
uks-hessian-o2-sto3g.

Historical **blocked** rows stay blocked (erf ERI, periodic HF, TD-DFTB, GFN2,
PM6/7, CASPT2).

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

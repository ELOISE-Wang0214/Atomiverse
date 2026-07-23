# PySCF method loop — handoff

Last harness pass: **2026-07-23 wave 4** on `dev`.

## Wave 4 summary

| Slug | Domain | Status | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{pyscf}}\) |
| ---- | ------ | ------ | ---------- | ----------------------------------------------- |
| `romp2-o2-sto3g` | post-hf | **merged → done** | `regressions/scf/romp2-o2-sto3g` | **1.37** |
| `roks-grad-o2-sto3g` | dft | **merged → done** | `regressions/scf/roks-grad-o2-sto3g` | **1.117** |
| `h2-b3lyp-rpa-tddft` | excited-states | **merged → done** | `regressions/scf/h2-b3lyp-rpa-tddft` | **0.430** |
| `spin-flip-tda-o2` | excited-states | **merged → done** | `regressions/scf/spin-flip-tda-o2` | **9.56** |
| `rohf-hessian-o2` | hartree-fock | **merged → done** | `regressions/scf/rohf-hessian-o2` | **perf_exempt** (no PySCF `hessian.ROHF`) |

- **Blocked this wave:** none
- **Merge-failed this wave:** none
- **Queue:** 8 pending; in_progress empty; historical blocked list unchanged
- **Roadmap:** CPU ✓ for **ROMP2** (`ROHF::mp2`), **ROKS analytic gradient**,
  **hybrid full-LR TDDFT** (B3LYP via `RKS::tddft`; pure-KS row already ✓),
  **ROHF nuclear Hessian** (FD of analytic forces — not CPHF). Spin-flip TDA
  table was already ✓; prose/index aligned. **No ⚡** — ratios noted only;
  no standing perf suite in ctest.
- **Speed notes:** all four fair-benched slugs within 10× (sf-tda closest at
  9.56); hybrid RPA faster than PySCF on H₂/STO-3G (0.43×); rohf-hessian
  exempt by design.

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
after serial merges in `fix(oracle): re-union pin_scf/fair_bench…`).

## Next assign candidates (pending head)

1. `adc2-h2co-sto3g` — algebraic-diagrammatic foothold
2. `mp3-h2-sto3g` — third-order PT
3. `ccsd-h2-sto3g` — coupled-cluster doubles foothold
4. `ccsd-t-h2o-sto3g` — needs CCSD
5. `casci-h2-sto3g` / `g0w0-h2-sto3g` — active-space / Green's-function footholds

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

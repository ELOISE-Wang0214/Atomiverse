# PySCF method loop — handoff

Last harness pass: **2026-07-23 wave 3** on `dev`.

## Wave 3 summary

| Slug | Domain | Status | Regression |
| ---- | ------ | ------ | ---------- |
| `sos-mp2-h2o-sto3g` | post-hf | **merged → done** | `regressions/scf/sos-mp2-h2o-sto3g` |
| `mp2-grad-h2o-sto3g` | post-hf | **merged → done** | `regressions/scf/mp2-grad-h2o-sto3g` |
| `uks-tda-o2-lsda` | excited-states | **merged → done** | `regressions/scf/uks-tda-o2-lsda` |
| `d3bj-molecular-pbe` | dft | **merged → done** | `regressions/scf/d3bj-molecular-pbe` |
| `uhf-hessian-o2` | hartree-fock | **merged → done** | `regressions/scf/uhf-hessian-o2` |

- **Blocked this wave:** none
- **Merge-failed this wave:** none
- **Queue:** 11 pending; in_progress empty; historical blocked list unchanged
- **Roadmap:** CPU ✓ for **MP2 analytic gradient** (`RHF::mp2_forces`); SOS-MP2 /
  D3BJ / UHF Hessian tables already ✓ — prose + index caveats aligned for
  UKS-TDA (`UKS::tda`), RMP2 grad, UHF CPHF Hessian; **RHF analytic Hessian**
  inventory fix (was ○ despite `reg_h2o_rhf_freq` on `dev`)

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
after serial merges in `fix(merge): restore multi-feature tree…`).

## Next assign candidates (pending head)

1. `romp2-o2-sto3g` — MP2 on ROHF reference
2. `adc2-h2co-sto3g` — algebraic-diagrammatic foothold
3. `spin-flip-tda-o2` — SF-TDA on UHF/ROHF
4. `rohf-hessian-o2` — ROHF analytic Hessian
5. `mp3-h2-sto3g` / `ccsd-h2-sto3g` — higher correlation footholds

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

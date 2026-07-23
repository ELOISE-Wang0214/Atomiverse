# PySCF method loop — handoff

Last harness pass: **2026-07-23 wave 2** on `dev`.

## Wave 2 summary

| Slug | Domain | Status | Regression |
| ---- | ------ | ------ | ---------- |
| `scs-mp2-h2o-sto3g` | post-hf | **merged → done** | `regressions/scf/scs-mp2-h2o-sto3g` |
| `full-lr-tddft-h2-pbe` | excited-states | **merged → done** | `regressions/scf/full-lr-tddft-h2-pbe` |
| `roks-o2-sto3g` | dft | **merged → done** | `regressions/scf/roks-o2-sto3g` |
| `static-polariz-h2o` | hartree-fock | **merged → done** | `regressions/scf/static-polariz-h2o` |
| `ief-pcm-h2o` | multiscale | **merged → done** | `regressions/scf/ief-pcm-h2o` |

- **Blocked this wave:** none
- **Merge-failed this wave:** none
- **Queue:** 16 pending; in_progress empty; historical blocked list unchanged
- **Roadmap:** CPU ✓ for SCS-MP2, full linear-response TDDFT (`RKS::tddft`),
  ROKS energy (gradient still ○), RHF static polarizability (CPHF prose +
  index); IEF-PCM already ✓ on multiscale table — intro/index caveats aligned

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
after serial merges in `fix(merge): restore multi-feature queue…`).

## Next assign candidates (pending head)

1. `sos-mp2-h2o-sto3g` — opposite-spin only scale of RMP2
2. `mp2-grad-h2o-sto3g` — RMP2 analytic gradient
3. `romp2-o2-sto3g` — MP2 on ROHF reference
4. `adc2-h2co-sto3g` — algebraic-diagrammatic foothold
5. `spin-flip-tda-o2` / `uks-tda-o2-lsda` — open-shell response

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

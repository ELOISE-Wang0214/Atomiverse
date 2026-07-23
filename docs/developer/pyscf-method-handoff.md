# PySCF method loop — handoff

Last harness pass: **2026-07-23 wave 1** on `dev`.

## Wave 1 summary

| Slug | Domain | Status | Regression |
| ---- | ------ | ------ | ---------- |
| `frozen-core-mp2-h2o` | post-hf | **merged → done** | `regressions/scf/frozen-core-mp2-h2o` |
| `cisd-h2co-sto3g` | excited-states | **merged → done** | `regressions/scf/cisd-h2co-sto3g` |
| `blyp-water-sto3g` | dft | **merged → done** | `regressions/scf/blyp-water-sto3g` |
| `ghf-grad-h3-sto3g` | hartree-fock | **merged → done** | `regressions/scf/ghf-grad-h3-sto3g` |
| `pcm-grad-h2o` | multiscale | **merged → done** | `regressions/scf/pcm-grad-h2o` |

- **Blocked this wave:** none
- **Merge-failed this wave:** none
- **Queue:** 21 pending; in_progress empty; historical blocked list unchanged
- **Roadmap:** CPU ✓ for Frozen-core MP2, CIS(D), BLYP, GHF analytic gradient,
  Implicit-solvent analytic gradient; GHF prose + index caveats aligned with
  landed grads / C-PCM free-energy forces / CIS(D)

## Oracle

Pins live in each regression `REFERENCE.md`. Offline pin path:
`tools/setup_pyscf_env.sh` + `tools/pyscf_oracle/pin_scf.py` (union restored
after serial merges in `fix(merge): restore multi-feature queue/oracle…`).

## Next assign candidates (pending head)

1. `scs-mp2-h2o-sto3g` / `sos-mp2-h2o-sto3g` — spin-component scale of RMP2
2. `mp2-grad-h2o-sto3g` — RMP2 analytic gradient
3. `full-lr-tddft-h2-pbe` — Casida A/B beyond TDA
4. `romp2-o2-sto3g` / `roks-o2-sto3g` — open-shell RO ladder
5. `adc2-h2co-sto3g` / response footholds

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.

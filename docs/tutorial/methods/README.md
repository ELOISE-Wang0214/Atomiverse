# Method handbook pages

**Style (current generation): technical report.**  
Dense, sectioned, validation-first. **Not** the narrative Theory→Practice chapters
in `docs/tutorial/ch*.md`. A house narrative style will be applied in a later
rewrite pass.

Produced/updated by workflow `pyscf-method-loop` **Tutorial** phase after each
successful merge wave (arg `write_tutorials`, default true).

## Required sections

| Section | Content |
| ------- | ------- |
| Summary | What / when / API / regression path |
| Method | Equations, assumptions, exclusions |
| Implementation | Minimal public C++; how to run `reg_<slug>` |
| Validation (vs PySCF) | Hard-coded pins from `REFERENCE.md` |
| Performance | Fair wall \(T_{\mathrm{ATV}}/T_{\mathrm{PySCF}}\) |
| Limitations | Scope + known debt |
| References | Paths + literature |

Template: [`_TEMPLATE.md`](_TEMPLATE.md).

## Naming

```text
docs/tutorial/methods/<slug>.md
```

`slug` = queue / regression directory name.

## Data sources (do not invent)

- Accuracy: `regressions/scf/<slug>/REFERENCE.md`
- Speed: `docs/developer/pyscf-method-perf.md` or `tools/pyscf_oracle/fair_bench.py`
- Oracle rules: `docs/tutorial/appendix-pyscf-oracle.md` (PySCF never in ctest)

## Index (campaign methods)

| Method | Page | Regression |
| ------ | ---- | ---------- |
| UMP2 | [ump2-o2-sto3g](ump2-o2-sto3g.md) | `regressions/scf/ump2-o2-sto3g` |
| Hybrid TDDFT-TDA (B3LYP) | [h2-b3lyp-hybrid-tda](h2-b3lyp-hybrid-tda.md) | `…/h2-b3lyp-hybrid-tda` |
| Frozen-core RMP2 | [frozen-core-mp2-h2o](frozen-core-mp2-h2o.md) | `…/frozen-core-mp2-h2o` |
| CIS(D) | [cisd-h2co-sto3g](cisd-h2co-sto3g.md) | `…/cisd-h2co-sto3g` |
| BLYP GGA | [blyp-water-sto3g](blyp-water-sto3g.md) | `…/blyp-water-sto3g` |
| GHF analytic gradient | [ghf-grad-h3-sto3g](ghf-grad-h3-sto3g.md) | `…/ghf-grad-h3-sto3g` |
| C-PCM free-energy forces | [pcm-grad-h2o](pcm-grad-h2o.md) | `…/pcm-grad-h2o` |
| SCS-MP2 | [scs-mp2-h2o-sto3g](scs-mp2-h2o-sto3g.md) | `…/scs-mp2-h2o-sto3g` |
| Full LR-TDDFT (RKS) | [full-lr-tddft-h2-pbe](full-lr-tddft-h2-pbe.md) | `…/full-lr-tddft-h2-pbe` |
| ROKS | [roks-o2-sto3g](roks-o2-sto3g.md) | `…/roks-o2-sto3g` |
| RHF static polarizability | [static-polariz-h2o](static-polariz-h2o.md) | `…/static-polariz-h2o` |
| IEF-PCM | [ief-pcm-h2o](ief-pcm-h2o.md) | `…/ief-pcm-h2o` |
| ROMP2 | [romp2-o2-sto3g](romp2-o2-sto3g.md) | `…/romp2-o2-sto3g` |
| ROKS analytic gradient | [roks-grad-o2-sto3g](roks-grad-o2-sto3g.md) | `…/roks-grad-o2-sto3g` |
| Hybrid full-LR TDDFT (B3LYP RPA) | [h2-b3lyp-rpa-tddft](h2-b3lyp-rpa-tddft.md) | `…/h2-b3lyp-rpa-tddft` |
| Spin-flip TDA (SF-CIS) | [spin-flip-tda-o2](spin-flip-tda-o2.md) | `…/spin-flip-tda-o2` |
| ROHF nuclear Hessian | [rohf-hessian-o2](rohf-hessian-o2.md) | `…/rohf-hessian-o2` |
| ADC(2) | [adc2-h2co-sto3g](adc2-h2co-sto3g.md) | `…/adc2-h2co-sto3g` |
| RMP3 | [mp3-h2-sto3g](mp3-h2-sto3g.md) | `…/mp3-h2-sto3g` |
| CCSD | [ccsd-h2-sto3g](ccsd-h2-sto3g.md) | `…/ccsd-h2-sto3g` |
| CCSD(T) | [ccsd-t-h2o-sto3g](ccsd-t-h2o-sto3g.md) | `…/ccsd-t-h2o-sto3g` |
| G0W0 | [g0w0-h2-sto3g](g0w0-h2-sto3g.md) | `…/g0w0-h2-sto3g` |
| B2PLYP double hybrid | [b2plyp-h2o](b2plyp-h2o.md) | `…/b2plyp-h2o` |
| CASCI | [casci-h2-sto3g](casci-h2-sto3g.md) | `…/casci-h2-sto3g` |
| ΔSCF (IMOM) | [delta-scf-n2](delta-scf-n2.md) | `…/delta-scf-n2` |
| EOM-EE-CCSD | [eom-ee-ccsd](eom-ee-ccsd.md) | `…/eom-ee-ccsd` |

Existing pages may still use an older Theory/Practice layout; next workflow run
or a bulk rewrite should move them to the technical-report template above.

# Method handbook pages

**Style (current generation): technical report.**  
Dense, sectioned, validation-first. **Not** the narrative Theory→Practice chapters
in `docs/tutorial/ch*.md`. A house narrative style will be applied in a later
rewrite pass.

Produced/updated by workflow **Tutorial** phases after each successful merge
wave:

- `pyscf-method-loop` — molecular post-HF / response (arg `write_tutorials`)
- `cp2k-method-loop` — GPW / PW / periodic DFT vs offline CP2K pins

## Required sections

| Section | Content |
| ------- | ------- |
| Summary | What / when / API / regression path |
| Method | Equations, assumptions, exclusions |
| Implementation | Minimal public C++; how to run `reg_<slug>` |
| Validation | Hard-coded pins from `REFERENCE.md` (PySCF or CP2K) |
| Performance | Fair wall \(T_{\mathrm{ATV}}/T_{\mathrm{oracle}}\) (or `perf_exempt`) |
| Limitations | Scope + known debt |
| References | Paths + literature |

Template: [`_TEMPLATE.md`](_TEMPLATE.md).

## Naming

```text
docs/tutorial/methods/<slug>.md
```

`slug` = queue / regression directory name.

## Data sources (do not invent)

**PySCF campaign (molecular):**

- Accuracy: `regressions/scf/<slug>/REFERENCE.md`
- Speed: `docs/developer/pyscf-method-perf.md` or `tools/pyscf_oracle/fair_bench.py`
- Oracle rules: `docs/tutorial/appendix-pyscf-oracle.md` (PySCF never in ctest)

**CP2K capability (GPW / PW / periodic):**

- Accuracy: `regressions/{gpw,pw}/<slug>/REFERENCE.md`
- Speed / ratios: `docs/developer/cp2k-method-handoff.md` (fair wall or `perf_exempt`)
- Oracle: `tools/cp2k_oracle/README.md` (offline pins only; CP2K never in ctest)
- Queue / design: `docs/developer/cp2k-method-queue.md`, `cp2k-method-loop.md`

## Index — CP2K capability wave 1 (2026-07-23)

| Method | Page | Regression |
| ------ | ---- | ---------- |
| GPW Γ-RKS H₂O PBE (MT) | [gpw-water-pbe-mt](gpw-water-pbe-mt.md) | `regressions/gpw/gpw-water-pbe-mt` |
| GPW analytic forces H₂O | [gpw-forces-water](gpw-forces-water.md) | `regressions/gpw/gpw-forces-water` |
| Γ-UKS GPW O₂ | [gpw-uks-o2-gamma](gpw-uks-o2-gamma.md) | `regressions/gpw/gpw-uks-o2-gamma` |
| GPW PBE-D3(BJ) ice dimer | [gpw-pbe-d3bj-ice](gpw-pbe-d3bj-ice.md) | `regressions/gpw/gpw-pbe-d3bj-ice` |
| PW LDA Si (k-mesh) | [pw-lda-si](pw-lda-si.md) | `regressions/pw/pw-lda-si` |
| PW true PBE GGA Si | [pw-pbe-gga-si](pw-pbe-gga-si.md) | `regressions/pw/pw-pbe-gga-si` |
| Metallic Fermi–Dirac Al | [pw-smearing-al](pw-smearing-al.md) | `regressions/pw/pw-smearing-al` |
| PW collinear spin O₂ | [pw-spin-polarized-o2](pw-spin-polarized-o2.md) | `regressions/pw/pw-spin-polarized-o2` |
| PW BOMD H (NVE foothold) | [pw-bomd-h](pw-bomd-h.md) | `regressions/pw/pw-bomd-h` |
| OT SCF water (GPW PBE) | [ot-minimizer-water](ot-minimizer-water.md) | `regressions/gpw/ot-minimizer-water` |

## Index — PySCF campaign methods

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
| Spin-flip full TDDFT (SF-RPA) | [sf-tddft-o2](sf-tddft-o2.md) | `…/sf-tddft-o2` |
| ROHF nuclear Hessian | [rohf-hessian-o2](rohf-hessian-o2.md) | `…/rohf-hessian-o2` |
| ADC(2) | [adc2-h2co-sto3g](adc2-h2co-sto3g.md) | `…/adc2-h2co-sto3g` |
| ADC(2)-x | [adc2x-h2co-sto3g](adc2x-h2co-sto3g.md) | `…/adc2x-h2co-sto3g` |
| IP-ADC(2) | [ip-adc2-h2o-sto3g](ip-adc2-h2o-sto3g.md) | `…/ip-adc2-h2o-sto3g` |
| EA-ADC(2) | [ea-adc2-h2-sto3g](ea-adc2-h2-sto3g.md) | `…/ea-adc2-h2-sto3g` |
| RMP3 | [mp3-h2-sto3g](mp3-h2-sto3g.md) | `…/mp3-h2-sto3g` |
| RMP4(SDTQ) | [mp4-h2-sto3g](mp4-h2-sto3g.md) | `…/mp4-h2-sto3g` |
| CCSD | [ccsd-h2-sto3g](ccsd-h2-sto3g.md) | `…/ccsd-h2-sto3g` |
| CC2 | [cc2-h2-sto3g](cc2-h2-sto3g.md) | `…/cc2-h2-sto3g` |
| CCSD(T) | [ccsd-t-h2o-sto3g](ccsd-t-h2o-sto3g.md) | `…/ccsd-t-h2o-sto3g` |
| G0W0 | [g0w0-h2-sto3g](g0w0-h2-sto3g.md) | `…/g0w0-h2-sto3g` |
| BSE@G0W0 | [bse-h2-sto3g](bse-h2-sto3g.md) | `…/bse-h2-sto3g` |
| B2PLYP double hybrid | [b2plyp-h2o](b2plyp-h2o.md) | `…/b2plyp-h2o` |
| CASCI | [casci-h2-sto3g](casci-h2-sto3g.md) | `…/casci-h2-sto3g` |
| ΔSCF (IMOM) | [delta-scf-n2](delta-scf-n2.md) | `…/delta-scf-n2` |
| EOM-EE-CCSD | [eom-ee-ccsd](eom-ee-ccsd.md) | `…/eom-ee-ccsd` |
| EOM-SF-CCSD | [eom-sf-ccsd-o2](eom-sf-ccsd-o2.md) | `…/eom-sf-ccsd-o2` |
| DFT-D4 (PBE params) | [dft-d4-h2o-pbe](dft-d4-h2o-pbe.md) | `…/dft-d4-h2o-pbe` |
| EOM-EA-CCSD | [eom-ea-ccsd](eom-ea-ccsd.md) | `…/eom-ea-ccsd` |
| One-shot GF2 (G0F2) | [gf2-h2-sto3g](gf2-h2-sto3g.md) | `…/gf2-h2-sto3g` |
| Direct RPA | [direct-rpa-h2-sto3g](direct-rpa-h2-sto3g.md) | `…/direct-rpa-h2-sto3g` |

Existing pages may still use an older Theory/Practice layout; next workflow run
or a bulk rewrite should move them to the technical-report template above.

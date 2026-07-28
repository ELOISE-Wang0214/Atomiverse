# PySCF fair-bench (post wave-1/2 stop)

Gate: r = T_ATV / T_PySCF ≤ 10. ATV Release regs; PySCF 2.14.0.

| slug | twin | T_pyscf [s] | T_atv [s] | ratio | ≤10× | note |
| ---- | ---- | ----------: | --------: | ----: | ---- | ---- |
| `ump2-o2-sto3g` | ump2 | 2.507 | 0.5712 | 0.2278 | ✓ | UMP2 twin |
| `h2-b3lyp-hybrid-tda` | rks_tda | 3.677 | 4.656 | 1.266 | ✓ | hybrid TDA |
| `frozen-core-mp2-h2o` | mp2 | 3.099 | 0.226 | 0.07291 | ✓ | MP2 twin |
| `scs-mp2-h2o-sto3g` | mp2 | 2.052 | 0.02558 | 0.01246 | ✓ | SCS vs plain MP2 |
| `blyp-water-sto3g` | rks | 6.567 | 14.87 | 2.264 | ✓ | RKS/BLYP |
| `full-lr-tddft-h2-pbe` | rks_tda | 5.456 | 1.451 | 0.2659 | ✓ | TDA twin for full RPA path |
| `roks-o2-sto3g` | rohf | 0.4506 | 0.5898 | 1.309 | ✓ | ROHF twin for ROKS |
| `static-polariz-h2o` | rhf | 0.6333 | 0.1286 | 0.203 | ✓ | RHF twin; polariz extra ATV |
| `ief-pcm-h2o` | rhf | 0.6617 | 6.006 | 9.077 | ✓ | vacuum RHF twin; IEF extra ATV |
| `cisd-h2co-sto3g` | tda | 0.8924 | 0.2184 | 0.2447 | ✓ | TDA twin; CIS(D) extra ATV |
| `ghf-grad-h3-sto3g` | uhf | 0.3195 | 0.02813 | 0.08803 | ✓ | proxy UHF/H2 vs GHF/H3 |
| `pcm-grad-h2o` | rhf | 0.6745 | 42.44 | 62.92 | ✗ | vacuum RHF twin; PCM forces extra |

**Summary:** ≤10× ✓=11  ✗=1  other=0.

## Caveats

- Reg binaries include full checks; PySCF side is twin kernel only.
- Extra ATV physics (PCM/CIS(D)/polariz/IEF) vs vacuum twin inflates ratio.
- loop-3 cancelled after wave 2; scheduler off; stopped.


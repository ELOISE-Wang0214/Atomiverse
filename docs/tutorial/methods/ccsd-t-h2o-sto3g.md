# CCSD(T) — H₂O / STO-3G

*(Technical-report draft; narrative rewrite later.)*

| Field | Value |
| ----- | ----- |
| Slug | `ccsd-t-h2o-sto3g` |
| Regression | `regressions/scf/ccsd-t-h2o-sto3g/` |
| Public API | `cpu::scf::RHF::ccsd_t(Frame&)` → `CCSDTResult` |
| Roadmap | [Post-HF §3 Coupled cluster](../../roadmap/post-hf.md) |

## Summary

Closed-shell **CCSD(T)** correlation energy on H₂O / STO-3G after a converged
RHF reference. Product surface is `RHF::ccsd_t(frame)`: drives a fresh RHF,
spin-orbital Gauss–Stanton CCSD amplitudes, then the non-iterative (T)
correction, and returns
`CCSDTResult{e_corr, e_tot, e_corr_ccsd, e_t, iterations, converged}` in
Hartree. Offline oracle pin is PySCF 2.14.0 `cc.CCSD` + `ccsd_t()` (never live
in ctest). End-to-end gate: `regressions/scf/ccsd-t-h2o-sto3g` → binary
`reg_ccsd-t-h2o-sto3g`. Fair-bench ratio
\(r = T_{\mathrm{ATV}}/T_{\mathrm{PySCF}} = 0.333\) (within 10×).

## Method

**Ansatz.** CCSD cluster operator \(T = T_1 + T_2\), then a non-iterative
triples correction \(E_{(T)}\) (Raghavachari 1989; Watts–Gauss–Bartlett 1993):

\[
E_{\mathrm{corr}}^{\mathrm{CCSD(T)}} = E_{\mathrm{CCSD}} + E_{(T)},
\qquad
E_{\mathrm{tot}} = E_{\mathrm{RHF}} + E_{\mathrm{corr}}.
\]

**Spin treatment.** Same spin-orbital Gauss–Stanton intermediates as
`RHF::ccsd`; (T) uses the spin-orbital formula of JCP **98**, 8718 (1993)
(PySCF `gccsd_t`, energy-matched to restricted `ccsd_t`).

**Key approximations / scope of this foothold:**

- Single-reference, closed-shell RHF only (no UCCSD(T) / ROCCSD(T)).
- Energy only — no \(\Lambda\), no densities, no analytic gradient.
- Canonical MO integrals; no DF/CD/F12/local CC.
- H₂O/STO-3G: 7 AO; 5 occupied + 2 virtual spatial MO.

## API

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::STO3G;
rhf.tol = 1e-12;
const cpu::scf::CCSDTResult ct = rhf.ccsd_t(frame);
// ct.e_corr = E_CCSD + E_(T); ct.e_t; ct.e_corr_ccsd [Ha]
```

## Source map

| Path | Role |
|---|---|
| `include/atomiverse/cpu/scf/ccsd.h` | Public `CCSDResult` / `CCSDTResult` |
| `include/atomiverse/cpu/scf/mean_field_method.h` | `RHF::ccsd` / `RHF::ccsd_t` |
| `include/atomiverse/cpu/scf/internal/ccsd_pack.h` | Private pack entries |
| `src/cpu/scf/ccsd.cpp` | Gauss–Stanton residual + (T) |
| `regressions/scf/ccsd-t-h2o-sto3g/` | E2E consumer + REFERENCE.md |

## Performance

| Twin | PySCF median (s) | ATV median (s) | ratio | ≤10× |
|---|---:|---:|---:|:---:|
| `ccsd_t` H₂O/STO-3G | 0.0933 | 0.0311 | **0.333** | ✓ |

```bash
export OMP_NUM_THREADS=1
~/work/.x86_64/bin/python tools/pyscf_oracle/fair_bench.py \
    --method ccsd_t --geom water --basis sto-3g --repeat 3 -- \
    ./build/regressions/reg_ccsd-t-h2o-sto3g
```

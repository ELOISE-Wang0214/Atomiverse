# Reference — plane-wave Nielsen–Martin analytic stress (diamond Si)

## Method

Cell stress from **Nielsen–Martin** Hellmann–Feynman cell derivatives after
**one SCF**: density is affinely remapped under symmetric strain (fixed
fractional coordinates), \(V_H\) / \(V_{xc}\) / \(V_\mathrm{loc}\) / Ewald
are rebuilt, one Davidson diagonalization is performed per strained cell,
and \(\sigma = (1/\Omega)\,\partial E_\mathrm{NSCF}/\partial\varepsilon\)
is formed by central differences. This is the numerical HF form of
Nielsen & Martin, *Phys. Rev. B* **32**, 3780 (1985) — **one SCF class**
cost, not 12 full SCFs.

Public path: `cpu::pw::pw_stress` → Voigt [eV/Å³].
Full multi-SCF strain-FD oracle: `cpu::pw::pw_stress_fd` (unit tests /
offline pin only).

## System (same twin as `pw-lda-si`)

| Setting | Value |
| ------- | ----- |
| Structure | diamond Si, 2-atom fcc primitive |
| Lattice \(a\) | 5.43 Å |
| Basis | (0,0,0), \((a/4,a/4,a/4)\) |
| PP | `tests/data/pw/Si_lda_nc.upf` (PseudoDojo ONCVPSP-LDA, \(Z_\mathrm{val}=4\)) |
| XC | LSDA (Slater + PW92) |
| \(E_\mathrm{cut}\) | 15 Ry (tutorial; incomplete vs production ~60 Ry) |
| k-mesh | \(2\times2\times2\) |
| SCF | `tol=10^{-6}`, `max_iter=50` |

## Offline dual-gate pin (multi-SCF strain-FD)

Measured **2026-07-24** on Arrhenius epyc (Release, `OMP_NUM_THREADS=16`)
via `pw_stress_fd` (\(h=2\times 10^{-3}\)) on this exact twin:

| Component | Value [eV/Å³] |
| --------- | ------------- |
| \(\sigma_{xx},\sigma_{yy},\sigma_{zz}\) | ≈ **−0.05755** (mean; components equal within \(3\times 10^{-5}\)) |
| \(\sigma_{xy},\sigma_{xz},\sigma_{yz}\) | ≈ **0** (crystal symmetry; measured \(\lesssim 0.02\) residual at incomplete ecut/k) |
| Live dual-gate \(\max_v\lvert\sigma^\mathrm{an}-\sigma^\mathrm{FD}\rvert\) | **\(7.3\times 10^{-5}\)** eV/Å³ |

The regression hard-codes the diagonal pin −0.05755 and off-diagonal 0 with
tol **0.05 eV/Å³** (and a soft symmetry tol **0.025 eV/Å³**). The pin is
**not** a freeze of the analytic output: it is the multi-SCF strain-FD
oracle measured offline; analytic agreed to \(7\times 10^{-5}\) when both
were run live.

## CP2K note

CP2K Quickstep bulk-Si decks use **GTH-PADE + DZVP-GTH**, not ONCV
plane-wave. Absolute `STRESS_TENSOR` values are **not** comparable to this
NC-UPF twin. No CP2K binary was available on the pin host this shift for a
fair PW wall twin. Optional skeleton: `tools/cp2k_oracle/inputs/si_lda_pw_note.inp`.

## Literature / PseudoDojo context

| Quantity | Value | Source |
| -------- | ----- | ------ |
| \(a_0\) (EOS) | ≈ 5.468 Å | PseudoDojo `Si.djrepo` deltafactor \(v_0\) |
| \(B_0\) | 89.96 GPa | same |
| Twin \(a\) | 5.43 Å (compressed vs \(a_0\)) | regression |

Compressed cell → compressive hydrostatic \(\sigma_{xx}<0\) (our
\(\sigma=(1/\Omega)\partial E/\partial\varepsilon\) convention), consistent
with the FD pin sign.

## Performance

| Item | Value |
| ---- | ----- |
| `perf_exempt` | **true** |
| Reason | No CP2K / QE / ABINIT binary on the pin host for a fair NC-UPF PW stress wall twin |
| ATV analytic wall (informative) | ≈ **30 s** (1 SCF + 12 NSCF diagonalizations; Release, 16 threads) |
| Live dual-gate (analytic + full FD) | ≈ 116 s offline only — not in ctest |
| `ratio_atv_over_cp2k` | n/a (exempt) |
| `within_10x` | n/a (exempt; dual-gate correctness only) |

When a plane-wave oracle binary is available: time the same geometry, append
\(r = T_\mathrm{atv}/T_\mathrm{oracle}\) here, require \(r \le 10\).

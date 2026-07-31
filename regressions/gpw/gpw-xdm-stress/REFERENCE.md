# Reference — Γ GPW analytic cell stress with free-atom XDM

## Method

Γ-point **GPW** (PBE) on a periodic water dimer with **free-atom XDM**
(exchange-hole dipole moment) pairwise dispersion — Becke–Johnson combination
rules, QE NC-PBE BJ damping \(a_1=0.4283\), \(a_2=2.4690\) Å, pure real-space
lattice sum with rcut **50 Å** (same kernel as `gpw-xdm-ice` /
`gpw-xdm-forces`).

**Product cell stress** is the Hellmann–Feynman fixed-density-matrix path:

1. One SCF at the unstrained cell → \(P\), energy-weighted \(W\).
2. Central FD of the non-self-consistent Lagrangian \(E_L[P;W;\varepsilon]\)
   under orthorhombic diagonal strain at fixed fractional coordinates, with
   the geometry-only \(E_{\mathrm{disp}}^{\mathrm{XDM}}(\{\mathbf R\},\mathbf h)\)
   **rebuilt at every strained cell** (same `xdm_pbe` call + rcut as the
   energy path).
3. \(\sigma_{aa} = -\Omega_0^{-1}\,\partial E_L/\partial\varepsilon_{aa}\)
   (CP2K `STRESS_TENSOR` sign; \(E\) in Hartree, \(\Omega_0\) in Bohr³).

Rebuilding at each strain is the whole point: holding \(E_{\mathrm{disp}}\)
fixed would differentiate a different functional than the SCF converged (the
0.53 Ha/Bohr ZnO failure mode). The fold itself shipped with the
dispersion-stress plumbing; this case adds the **dedicated pins** (queue row:
energy folds, dedicated pin missing).

Public path: `cpu::gpw::GPWEvaluator` with `xdm = true`,
`compute_stress = true`, `stress_use_fd = false` → Frame `kStress` [eV/Å³];
compared here via `last_stress_ha_bohr3()`.

## Oracles

| Gate | Oracle |
| ---- | ------ |
| A. \(E_{\mathrm{disp}}\) unstrained | offline **dual pure-Python** free-atom XDM lattice sum — bit-identical provenance with the shipped `gpw-xdm-ice` pin |
| B. \(E_{\mathrm{disp}}\) **under strain** (±h per axis, fixed fractional coords, rcut 50 Å) | fresh dual pure-Python pins at each of the 6 strained cells (below) — the strain derivative being folded is externally validated point by point |
| C. \(\sigma\) | **multi-SCF strain-FD** (`gpw_stress_fd_impl`): full re-SCF at each strained cell, XDM folded through the energy path |
| D. XDM actually folded, with the pinned magnitude | \(\sigma(\mathrm{xdm}) - \sigma(\mathrm{plain})\) equals the dual-python strain-FD term \(-\Delta E_{\mathrm{disp}}/(2h\Omega_0)\) per axis |

### Why not CP2K QS?

**CP2K 2026.1** stock `PAIR_POTENTIAL TYPE` is only
`DFTD2|DFTD3|DFTD3(BJ)|DFTD4` — there is **no** stock `TYPE XDM`, so no CP2K
stress twin exists. Hard pins are the dual pure-Python energy sums (same
policy and oracle family as `gpw-xdm-ice` / `gpw-xdm-forces`), plus the
multi-SCF strain-FD dual gate for the tensor — same structure as
`gpw-mbd-stress` / `gpw-mbd-kpoint-stress`.

## System

Same geometry as `gpw-xdm-ice` / `gpw-mbd-stress`:

| Item | Value |
|------|-------|
| System | water dimer, **6 atoms**, order **O H H O H H** |
| Cell | 12 × 12 × 12 Å, `PERIODIC XYZ` |
| Coordinates (Å) | O (0,0,0); H (0.758602, 0, ±0.504284); O (3,0,0); H (3.758602, 0, ±0.504284) |
| Basis | `DZVP-MOLOPT-SR-GTH` |
| PP | O `GTH-PBE-q6`, H `GTH-PBE-q1` → **16** valence e⁻ (`n_occ = 8`) |
| XC | PBE |
| Dispersion | free-atom XDM, QE NC-PBE \(a_1=0.4283\), \(a_2=2.4690\) Å, rcut 50 Å |
| Cutoff | 200 Ry (E_disp is geometry-only — cutoff does not enter the pins) |
| Strain step | \(h = 2\times10^{-3}\) |

## Pins

Source: offline **dual pure-Python** free-atom XDM lattice sum
(`tools/cp2k_oracle/inputs/gpw-xdm-stress/pin_xdm_strain.py`, python3 stdlib
only; independent implementation sharing only the cited free-atom data:
CRC 88th Ed. α, homoatomic \(C_6/C_8/C_{10}\) from Yan et al. PRA **54**,
2824 and Chu & Dalgarno JCP **121**, 4083). Never launched at ctest time.
Output lines (2026-07-31):

```
E_disp unstrained            = -2.1148157612879029e-03 Ha
E_disp x +h              = -2.1110849131940963e-03 Ha
E_disp x -h              = -2.1185431863446984e-03 Ha
sigma_XDM_xx (FD, dual-py)  = -1.5989612853e-07 Ha/Bohr3
E_disp y +h              = -2.1147774136190569e-03 Ha
E_disp y -h              = -2.1148546246483736e-03 Ha
sigma_XDM_yy (FD, dual-py)  = -1.6553087314e-09 Ha/Bohr3
E_disp z +h              = -2.1146664855098696e-03 Ha
E_disp z -h              = -2.1149651549478578e-03 Ha
sigma_XDM_zz (FD, dual-py)  = -6.4031024180e-09 Ha/Bohr3
```

Provenance cross-check at pin time: the unstrained value **bit-reproduces**
the shipped `gpw-xdm-ice` periodic pin (−2.1148157612879029×10⁻³ Ha).

## Tolerances

| Check | Tol |
| ----- | --- |
| A/B \(\lvert E_{\mathrm{disp}} - \mathrm{pin}\rvert\) (unstrained + 6 strained) | 1×10⁻¹² Ha (deterministic geometry-only sum) |
| C \(\sigma_{\mathrm{diag}}\) analytic vs multi-SCF FD | 2×10⁻⁵ Ha/Bohr³ (Γ dispersion-stress reg tier) |
| C \(\lvert\sigma\rvert\) nonzero | > 1×10⁻¹⁰ |
| D \(\lvert(\sigma_{\mathrm{xdm}}-\sigma_{\mathrm{plain}}) - \sigma^{\mathrm{XDM}}_{\mathrm{dual\text{-}py}}\rvert\) | 1×10⁻⁸ Ha/Bohr³ |
| E off-diagonals (orthorhombic) | exactly 0 |
| E `kStress` unit conversion | 1×10⁻⁹ |
| Product wall | ≤ 60 s |

## Tolerance justification

- **A/B 1e−12 Ha:** the product kernel and the dual pure-Python sum are two
  implementations of the same deterministic geometry-only formula on shared
  cited data; the shipped `gpw-xdm-ice` dual matched to ~1e−16.
- **C 2e−5 Ha/Bohr³:** analytic fixed-\(P\) vs multi-SCF FD differ by the
  frozen density-matrix approximation plus SCF residual at `tol = 1e-6`, not
  by the dispersion treatment — same tier as `gpw-d4-stress` /
  `gpw-mbd-stress`.
- **D 1e−8 Ha/Bohr³:** both sides are central FDs at the same \(h\) of the
  same dispersion energy (XDM does not enter the Fock operator, so the
  electronic parts cancel); the gate only has to exclude "silently not
  folded" (σ_XDM_xx ≈ 1.6e−7) and a wrong strained-cell rebuild.

## Dual gate

- **Correctness:** the gates above, plus unit gates in
  `tests/cpu/gpw/test_cpu_gpw_stress.cpp` — `AnalyticMatchesFdWithXdm`
  (cheap MT H₂, sibling parity with D4/D3BJ/MBD) and
  `XdmChangesGammaDiagonal`.
- **Speed:** no XDM-capable oracle binary exists to time against (CP2K stock
  has no XDM), so this case is **perf_exempt** — same policy as
  `gpw-xdm-ice` / `gpw-xdm-forces` / `gpw-mbd-stress`.

## Measured results (2026-07-31, first green)

| Gate | Tol | Measured |
| ---- | --- | -------- |
| B \(E_{\mathrm{disp}}\) under strain vs dual-py (6 cells) | 1e−12 Ha | TBD |
| A \(E_{\mathrm{disp}}\) product vs pin | 1e−12 Ha | TBD |
| C \(\sigma_{\mathrm{diag}}\) analytic vs multi-SCF FD | 2e−5 Ha/Bohr³ | TBD |
| D shift vs dual-py FD term | 1e−8 Ha/Bohr³ | TBD |
| E off-diagonals | exactly 0 | TBD |
| E `kStress` conversion | 1e−9 | TBD |
| Product wall | ≤ 60 s | TBD |

## Notes

- Sibling XDM products: `gpw-xdm-ice` (energy), `gpw-xdm-forces` (ionic
  forces). Sibling dispersion stress: `gpw-d3bj-stress`, `gpw-d4-stress`,
  `gpw-mbd-stress`, `gpw-mbd-kpoint-stress`.
- Multi-k XDM stress folds through the same k-summed Lagrangian
  (`energy_fixed_P_k` dispersion branch) but has no dedicated pin case yet.

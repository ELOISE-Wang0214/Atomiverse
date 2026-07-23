# Reference — PW band structure + DOS (bulk silicon)

## Method

Spin-unpolarized plane-wave Kohn–Sham DFT with LDA (Slater + PW92) on the PW
real-space grid. Norm-conserving ONCVPSP LDA Si UPF. Self-consistent density
from a uniform Monkhorst–Pack mesh, then **non-self-consistent**
diagonalization on a high-symmetry **k-path** (bands) and a denser MP mesh
(DOS, Gaussian broadening).

Public path:

```cpp
cpu::pw::PW method;          // members: upf, ecut, k_mesh, xc=LSDA, …
cpu::pw::Bands bands;        // members: n_bands, kfrac / set_diamond_path
auto br = bands.compute(method, frame);
cpu::pw::Dos dos;            // members: k_mesh, n_bands, sigma_ry, …
auto dr = dos.compute(method, frame);
```

No public `*Config` bag; no factories.

## Pseudopotential

| Item | Value |
| ---- | ----- |
| File | `tests/data/pw/Si_lda_nc.upf` (in-tree) |
| Family | PseudoDojo ONCVPSP-LDA-SR-noCC v0.4 |
| Cite | Hamann, *Phys. Rev. B* **88**, 085117 (2013); van Setten et al., *Comput. Phys. Commun.* **226**, 39 (2018) |

Same PP family as `regressions/pw/pw-lda-si`.

## Literature pin (indirect Γ→X gap)

Bulk Si is an **indirect** gap semiconductor: valence-band maximum at
\(\Gamma\), conduction-band minimum near \(X\).

| Quantity | Value | Source |
| -------- | ----- | ------ |
| Experimental \(E_g\) | ~1.12 eV | textbook |
| Converged LDA \(E_g\) (Γ\(_v\)→ near-\(X_c\)) | **~0.45 eV** | Yin & Cohen, *Phys. Rev. B* **26**, 5668 (1982) |

The regression pins the **textbook Si transition**
\[
E_{\Gamma\to X}
  = \varepsilon_{c}(X) - \varepsilon_{v}(\Gamma)
\]
(not the global \(\max\varepsilon_v - \min\varepsilon_c\) over a sparse path,
which at incomplete \(E_\mathrm{cut}\) can spuriously go slightly negative
from path-point scatter while every **local** same-\(k\) gap stays open).

## Twin geometry and settings

| Setting | Value |
| ------- | ----- |
| Structure | diamond Si, 2-atom fcc primitive cell |
| Lattice constant \(a\) | 5.43 Å |
| Basis atoms | (0,0,0), \((a/4,a/4,a/4)\) |
| SCF \(E_\mathrm{cut}\) | 15 Ry |
| SCF k-mesh | \(2\times2\times2\) (`PW::k_mesh = 2`) |
| XC | `XCType::LSDA` |
| SCF tol | \(10^{-6}\) |
| Bands | diamond path Γ–X–W–L–Γ–K, `n_seg=4`, `n_bands=8` |
| DOS | MP `k_mesh=4`, `n_bands=8`, \(\sigma=0.05\) Ry, 200 bins |

## Acceptance (incomplete basis)

1. SCF for both `Bands` and `Dos` converges.
2. Eigenvalues per k sorted ascending; 8 bands.
3. Every **local** (same-\(k\)) gap \(\varepsilon_{n_\mathrm{occ}}(k)-\varepsilon_{n_\mathrm{occ}-1}(k) > 0\).
4. Path includes \(\Gamma\) and \(X\); **indirect character**
   \(\varepsilon_c(X) < \varepsilon_c(\Gamma)\).
5. **Γ→X gap window** (eV):
   \[
   0.05 \le E_{\Gamma\to X}^\mathrm{ATV} \le 2.50
   \]
   with \(E_{\Gamma\to X}=(\varepsilon_c(X)-\varepsilon_v(\Gamma))\times 13.605693\).
   Measured on this host (2026-07-24, Release): \(\approx 0.48\) eV —
   near Yin–Cohen ~0.45 eV despite tutorial ecut/k.
6. Direct \(\Gamma\) gap \(\ge\) Γ→X − 0.05 eV.
7. DOS occupied integral \(n_\mathrm{elec\_to\_vbm}\) within **1.0** of
   \(n_e = 8\).

## CP2K note

CP2K Quickstep `BAND_STRUCTURE` / DOS use **GTH + Gaussian** (DZVP-GTH), not
ONCV plane waves. Absolute \(\varepsilon_n\) are **not** comparable.
Character (indirect Si, \(\varepsilon_c(X)<\varepsilon_c(\Gamma)\)) and the
Yin–Cohen LDA Γ→X window are the scientific pins. Offline twin input:
`tools/cp2k_oracle/inputs/si_bands_dos_note.inp`.

## Performance

| Item | Value |
| ---- | ----- |
| `perf_exempt` | **true** |
| Reason | No CP2K binary / module on the build host for a fair wall twin this shift |
| ATV twin wall (informative) | ~22 s (Release, 8 OMP threads, ecut 15 / k=2 + path + DOS k=4) |
| `ratio_atv_over_cp2k` | n/a (exempt) |
| `within_10x` | n/a (exempt; correctness gate only) |

When a CP2K binary is available: time the same geometry with the oracle input,
append \(r = T_\mathrm{atv}/T_\mathrm{cp2k}\) here, require \(r \le 10\).

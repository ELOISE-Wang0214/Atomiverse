---
slug: dos-spin-metallic-si
status: approved
created: 2026-07-26
domain: pw / bands_dos
queue: cp2k-method-loop
wave: 12
depends_on: [dos-spin-si, dos-metallic-si]
unblocks: [dos-spin-metallic-si]
sibling: bands-spin-metallic-si
---

# Spec: `dos-spin-metallic-si` — collinear spin + metallic smear DOS (PW)

**Slug:** `dos-spin-metallic-si`  
**Domain:** Plane-wave DOS post-processing (`cpu::pw::Dos`)  
**Queue:** CP2K method loop (wave 12)  
**Status:** approved — implement this shift  
**Sibling:** `bands-spin-metallic-si` (shares `bands_dos`; this slug enables
**Dos only** — leave `Bands::compute` spin+smear reject for the bands worker)

## Summary

Close the remaining `Dos::compute` foothold that throws when
`method.spin_polarized && method.smearing != none`. Product path:

1. Collinear UKS SCF with Fermi–Dirac (or MP/Gaussian) metallic occupations
   — shared \(\mu\) over \(\alpha+\beta\) spin-orbitals, \(f\in[0,1]\).
2. Freeze \(V_{\mathrm{eff},\sigma}\); NSCF on a denser MP mesh per spin.
3. Gaussian DOS with spin-orbital weight 1: \(D=D_\alpha+D_\beta\).
4. Recompute \(\mu\) on the DOS mesh (`spin_degen=1`); occupied integral
   \(\int^{E_F} D(E)\,dE \approx n_e\).

Twin: diamond Si (gapped), \(M=0\) seed, FD smear — same geometry class as
`dos-spin-si` / `dos-metallic-si`. Absolute \(D(E)\) not comparable to CP2K
GTH-GPW; character pins + path enablement.

## Domain basis

Collinear metallic occupations (QE / CP2K SMEAR + LSD):

\[
\sum_{\sigma=\uparrow,\downarrow}\sum_{n\mathbf{k}}
  w_{\mathbf{k}}\, f\!\left(\frac{\varepsilon_{n\mathbf{k}\sigma}-\mu}{\sigma}\right)
  = N_e,
\quad f\in[0,1].
\]

DOS post-processing (Gaussian broadening \(\sigma_{\mathrm{DOS}}\)):

\[
D_\sigma(E)=\sum_{n\mathbf{k}} w_{\mathbf{k}}\,
  g_\sigma\,
  \frac{1}{\sqrt{2\pi}\,\sigma_{\mathrm{DOS}}}
  \exp\!\Big(-\frac{(E-\varepsilon_{n\mathbf{k}\sigma})^2}{2\sigma_{\mathrm{DOS}}^2}\Big),
\quad g_\sigma=1.
\]

Occupied integral uses \(E_F=\mu\) from the DOS mesh (not SCF \(\mu\)):

\[
n_{\mathrm{elec}}(E_F)=\int_{-\infty}^{E_F} D(E)\,dE,\quad D=D_\uparrow+D_\downarrow.
\]

## Public API (no Config / no factories)

```cpp
cpu::pw::PW method;
method.spin_polarized = true;
method.multiplicity = 1;           // seeds ρ↑/ρ↓; M free under smear
method.smearing = Smearing::fermi_dirac;
method.degauss = 0.02;             // Ry
method.extra_bands = 4;
// ecut, k_mesh, upf, xc=LSDA …

cpu::pw::Dos dos;
dos.k_mesh = 4;
dos.n_bands = 8;
dos.sigma_ry = 0.05;
const DosResult dr = dos.compute(method, frame);
// dr.nspin==2; e_fermi_ry; n_elec_to_vbm ≈ 8 (integral to E_F)
// dos / dos_alpha / dos_beta
```

## Files to touch

| Path | Change |
| ---- | ------ |
| `include/.../bands_dos.h` | Document spin+smear DOS product path |
| `src/cpu/planewave/bands_dos.cpp` | Remove Dos throw; collinear smear SCF + μ on DOS mesh |
| `tests/cpu/planewave/test_pw_bands_dos.cpp` | Unit: pure-math / result defaults only (no e2e SCF) |
| `regressions/pw/dos-spin-metallic-si/` | E2E regression + REFERENCE.md |
| `regressions/CMakeLists.txt` | Register `reg_dos-spin-metallic-si` |
| `tools/cp2k_oracle/inputs/dos-spin-metallic-si.note` | Offline CP2K twin note |

**Out of scope:** `Bands::compute` spin+smear (sibling `bands-spin-metallic-si`);
CUDA; GPW DOS.

## Acceptance

See `acceptance.md`. Green ⇔ dual gate + commit on `wf/dos-spin-metallic-si`.

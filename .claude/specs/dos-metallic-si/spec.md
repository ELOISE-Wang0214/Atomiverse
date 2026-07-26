# dos-metallic-si — metallic occupation DOS (PW)

| Field | Value |
| ----- | ----- |
| Slug | `dos-metallic-si` |
| Method domain | `pw` |
| Branch | `wf/dos-metallic-si` |
| Parent | `dev` |
| Oracle | CP2K offline (`tools/cp2k_oracle/`) |
| Status | in progress |

## Problem

`cpu::pw::Dos::compute` throws when `method.smearing != Smearing::none`:

```text
Dos::compute: metallic smearing not supported for DOS foothold
```

Insulator DOS (`reg_bands-dos-si`) is closed. Metallic occupation DOS is an open
roadmap residual (`docs/roadmap/dft.md` § bands/DOS; queue `dos-metallic-si`).

## Goal

Productize **Gaussian-broadened DOS after a metallic-smearing SCF** on the same
public `cpu::pw::Dos` surface (no new `*Config` bag, no factories):

1. SCF uses restricted Fermi–Dirac (or MP/Gaussian) occupations when
   `method.smearing != none` (`run_scf_smear_restricted`).
2. NSCF on a denser MP mesh freezes \(V_\mathrm{eff}\) and diagonalizes
   `Dos::n_bands` states per k.
3. Chemical potential \(\mu\) is recomputed on the DOS mesh from the NSCF
   spectrum (`find_fermi_level`); report as `DosResult::e_fermi_ry`.
4. Occupied integral \(\int_{-\infty}^{\mu} D(E)\,dE \approx n_e\)
   (field remains `n_elec_to_vbm` for API continuity; for metallic runs it is
   the integral to \(E_F\)).
5. Fixed-occupation path (`smearing == none`) is **unchanged** (still VBM
   integral + gap).

## Non-goals

- Collinear-spin DOS (`bands-spin-si` / spin-resolved DOS — separate slug).
- Absolute eigenvalue match to CP2K GTH-GPW (basis/PP differ).
- Changing `Bands::compute` metallic path (optional later; keep throw or
  leave as-is — this slug is DOS-only).

## Public API

```cpp
cpu::pw::PW method;
method.smearing = Smearing::fermi_dirac;
method.degauss = 0.02;       // Ry
method.extra_bands = 4;
// … ecut, k_mesh, upf, xc …

cpu::pw::Dos dos;
dos.k_mesh = 4;
dos.n_bands = 8;
dos.sigma_ry = 0.05;
const DosResult dr = dos.compute(method, frame);
// dr.e_fermi_ry, dr.n_elec_to_vbm ≈ n_e, dr.dos, dr.scf_converged
```

## Twin (regression)

Diamond Si, 2-atom fcc primitive, \(a=5.43\) Å, ONCV LDA Si UPF, LSDA,
SCF \(E_\mathrm{cut}=15\) Ry, SCF k = \(2\times2\times2\), FD `degauss=0.02` Ry,
`extra_bands=4`, DOS mesh \(4^3\), \(\sigma=0.05\) Ry. Same geometry class as
`bands-dos-si` / `dfpt-metallic-si`.

## Acceptance

See `acceptance.md`.

## Files (expected)

| Path | Action |
| ---- | ------ |
| `include/.../bands_dos.h` | `DosResult::e_fermi_ry`; docs for metallic path |
| `src/.../bands_dos.cpp` | smear SCF + \(\mu\) + integral to \(E_F\) |
| `tests/cpu/planewave/test_pw_bands_dos.cpp` | unit: integral-to-mu helper / metal path math |
| `regressions/pw/dos-metallic-si/` | e2e + `REFERENCE.md` |
| `tools/cp2k_oracle/inputs/dos-metallic-si.note` | offline twin note |
| `regressions/CMakeLists.txt` | register `reg_dos-metallic-si` |

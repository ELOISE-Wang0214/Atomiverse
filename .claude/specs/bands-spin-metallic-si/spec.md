# Spec — bands-spin-metallic-si

## Goal

Enable **collinear spin-polarized** plane-wave band structure with **metallic
occupation smearing** on one public path:

```cpp
cpu::pw::PW method;
method.spin_polarized = true;
method.multiplicity = 1;           // seed M; free under smear
method.smearing = Smearing::fermi_dirac;
method.degauss = 0.02;             // Ry
method.extra_bands = 4;
// … ecut, k_mesh, upf, xc …

cpu::pw::Bands bands;
bands.n_bands = 8;
bands.set_diamond_path(4);
const BandsResult br = bands.compute(method, frame);
// br.nspin == 2; br.eigenvalues_ry (α), br.eigenvalues_beta_ry (β)
// br.e_fermi_ry from SCF; br.magnetization free under smear
```

Today `Bands::compute` throws when both `spin_polarized` and
`smearing != none` are set. Restricted metallic bands
(`bands-metallic-si`) and fixed-occupation collinear spin bands
(`bands-spin-si`) already ship.

## Out of scope

- Spin+smear **DOS** (`dos-spin-metallic-si` — sibling slug).
- Absolute ONCV eigenvalue match to CP2K GTH+Gaussian.
- Non-collinear magnetism, SOC, hybrid XC on the path.

## Design

1. **SCF:** collinear UKS + shared chemical potential \(\mu\) over spin
   orbitals (`f ∈ [0,1]` per `(σ,k,b)`), same algorithm as `PW::run`
   spin+smear. Seed \(N_↑,N_↓\) from `multiplicity`; \(M\) free under
   smear. Leave \(\rho_σ\), \(V_{\mathrm{xc},σ}\) on the problem for
   NSCF freeze.
2. **NSCF:** freeze \(V_{\mathrm{eff},σ}\) per spin, rebuild k-tables on
   the high-symmetry path, Davidson diagonalize each channel (same as
   fixed-occ spin bands).
3. **Reporting:** `e_fermi_ry` = SCF \(\mu\); gap from integer-fill
   reference occupations (multiplicity seed) when path has empty
   states — Si twin stays gapped with near-integer FD fill.
4. **API:** members only; no public `*Config`; no factories.

## Files

| Path | Change |
| ---- | ------ |
| `include/.../internal/restricted_scf.h` | declare `run_collinear_smear_scf` |
| `src/.../internal/restricted_scf.cpp` | implement serial spin+smear SCF |
| `include/.../bands_dos.h` | docs: spin+smear product |
| `src/.../bands_dos.cpp` | `Bands::compute` accepts spin+smear |
| `tests/cpu/planewave/test_pw_bands_dos.cpp` | drop reject; gate/acceptance unit |
| `regressions/pw/bands-spin-metallic-si/` | e2e + REFERENCE |
| `tools/cp2k_oracle/inputs/bands-spin-metallic-si.note` | offline twin |
| `regressions/CMakeLists.txt` | register reg |

## Twin

Diamond Si, \(a=5.43\) Å, ONCV LDA Si UPF, ecut 15 Ry, SCF k=2,
FD `degauss=0.02` Ry, `extra_bands=4`, `multiplicity=1`, diamond path
`n_seg=4`, `n_bands=8`. Scientific pins: path runs; open gap; \(E_F\)
finite near gap; \(\alpha\approx\beta\) and \(|M|\) small; Γ→X indirect
character (Yin–Cohen window). CP2K note is character/wall twin only
(GTH≠ONCV).

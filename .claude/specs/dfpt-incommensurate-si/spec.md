# Spec — `dfpt-incommensurate-si`

**Loop:** CP2K method loop  
**Slug:** `dfpt-incommensurate-si`  
**Branch:** `wf/dfpt-incommensurate-si`  
**Domain:** response / plane-wave DFPT  
**Oracle:** offline CP2K DFPT / PHONON twin note (or literature + self-consistent
energy pin when CP2K binary unavailable)  
**Dep:** `dfpt-finite-q-si` ✓ (commensurate supercell path); `dfpt-phonon-si` ✓

## Intent

Ship **incommensurate-q DFPT phonons** for `cpu::pw::DfptPhonon`: when
`q_frac` is **not** commensurate with a small supercell (product ≤ 8 along
axes ≤ 6), run the **direct finite-q Sternheimer** assembly already present as
`internal::dfpt_q_force_constants` (Baroni et al. RMP 2001), mass-weight the
Hermitian \(C(q)\), and return \(\omega(q)\) [cm⁻¹] plus Re/Im force constants
[eV/Å²].

Today the public facade throws
`"direct incommensurate Sternheimer is a follow-on"` for non-rational \(q\).
The internal engine (`dfpt_q_force_constants`) is already implemented for
arbitrary fractional \(q\); this slug **wires the public path** and pins a
minimal diamond-Si twin.

## Non-goals

- No GPW DFPT, no metallic smearing, no UKS, no GGA response kernel.
- No public `*Config` bag; no `make_*` / `build_*` factories.
- No live CP2K in ctest — offline pins only in `REFERENCE.md`.
- No change to the proven **commensurate** supercell+IFC path for rational \(q\).
- No CUDA port.

## Public API (members only)

```cpp
cpu::pw::DfptPhonon dfpt;
dfpt.q_frac = {0.3, 0.0, 0.0};   // incommensurate under nmax=6
auto r = dfpt.compute(method, frame);
// r.frequencies [cm⁻¹], r.force_constants / _imag [eV/Å²]
```

### Path selection in `DfptPhonon::compute`

| `q_frac` | Path |
| -------- | ---- |
| \(\approx\mathbf{0}\) | `dfpt_gamma_force_constants` (unchanged) |
| rational with \(n_a n_b n_c\le 8\), \(n_i\le 6\) | supercell Γ-DFPT + IFC Fourier (unchanged) |
| otherwise (incommensurate or oversize supercell) | **`dfpt_q_force_constants` + Hermitian mass-weighted diagonalization** |

## Acceptance

1. **Unit:** `test_dfpt` — incommensurate `q_frac` no longer fails at the
   commensurate gate; existing guards (empty UPF, non-LSDA, smearing, spin,
   bad mix) still throw **before** SCF for incommensurate \(q\).
2. **Regression:** `regressions/response/dfpt-incommensurate-si/` ≤ 60 s:
   diamond Si, NC-UPF, LSDA, ecut 10 Ry, k=1 or 2, \(q=(0.3,0,0)\).
   - 6 finite frequencies; literature-biased \|ω\| windows (Giannozzi PRB 1991
     Γ–X branch shape, large minimal-settings bias allowed).
   - Hermitian asymmetry of \(C(q)\) ≤ 1e-3 eV/Å².
   - SCF + response flags recorded; soft modes ≤ 2 if ecut/k limited.
3. **REFERENCE.md** cites method, settings, offline pin values, perf_exempt
   or wall ratio.
4. **Architecture:** no public Config; Frame I/O only; energy units eV on Frame
   when written (spectra stay on `DfptPhononResult`).
5. **Perf:** \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\le 10\) or `perf_exempt`
   with reason (no CP2K DFPT wall on host).

## Files (expected)

| Path | Change |
| ---- | ------ |
| `include/.../dfpt_phonon.h` | Docs: incommensurate path shipped |
| `src/cpu/planewave/dfpt_phonon.cpp` | Fall through to `dfpt_q_force_constants` |
| `tests/cpu/planewave/test_dfpt.cpp` | Incommensurate gate unit tests |
| `regressions/response/dfpt-incommensurate-si/` | reg + REFERENCE |
| `regressions/CMakeLists.txt` | `atv_add_regression` + UPF def |
| `tools/cp2k_oracle/inputs/dfpt-incommensurate-si.note` | Offline twin note |

## Science

- Baroni, de Gironcoli, Dal Corso, Giannozzi, *Rev. Mod. Phys.* **73**, 515 (2001).
- Giannozzi et al., *Phys. Rev. B* **43**, 7231 (1991) — Si phonon dispersion.

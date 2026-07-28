# Spec: dftb2-dispersion-d3

## Goal

Add **Grimme DFT-D3** (Becke–Johnson damping) as an opt-in, geometry-only
dispersion correction on `cpu::semiempirical::DFTB2`, for a water-dimer
product gate under SE-only cmake.

## Scope

**IN**

- Public member `DFTB2::d3` (bool, default `false`) — no `*Config` bag.
- When `d3 == true`, evaluate two-body D3(BJ) from static in-tree Grimme
  tables (reuse `gpw::dftd3_pbe_bj` / `D3Params` engine; compile those two
  sources under `ATV_BUILD_SEMIEMPIRICAL` when HF is off).
- Fold \(E_{\mathrm{disp}}\) [eV] into total energy written to
  `frame.meta[kEnergy]` and into FD force path when `with_forces`.
- Accessor `last_e_disp()` [eV] (0 when `d3` off).
- Damping parameters are the published **DFTB3-D3(BJ)** set used by DFTB+
  for tight-binding D3 (s6, s8, a1, a2) — public members with those defaults
  so a later method-specific set can override without a Config bag.
- Unit test(s) under `tests/cpu/semiempirical/`.
- Regression `regressions/semiempirical/dftb2-dispersion-d3/` with
  `REFERENCE.md` (`oracle: literature` via offline **Simple-DFTD3** /
  s-dftd3 pin of pure \(E_{\mathrm{disp}}\); DFTB+ not required for the
  geometry-only correction).

**OUT**

- Three-body ATM / C9.
- Periodic lattice-sum D3 (molecular isolated water dimer only).
- Changing DFTB2 SK / SCC electronics.
- Analytic D3 forces separate from existing re-SCC(+D3) FD of total E.
- HF / GPW / MD / MINIMIZER product flags.

## Public API

```cpp
cpu::semiempirical::DFTB2 dftb;
dftb.d3 = true;           // opt-in D3(BJ)
// optional scale overrides (DFTB3-D3(BJ) defaults):
// dftb.d3_s6, dftb.d3_s8, dftb.d3_a1, dftb.d3_a2
eval.on(evaluate::Phase::evaluate).add(std::move(dftb));
// frame.meta[kEnergy] = E_elec+rep + E_disp  [eV]
// dftb.last_e_disp() [eV]
```

## Science

- Grimme et al., *J. Chem. Phys.* **132**, 154104 (2010) — D3 model.
- BJ damping: Grimme et al., *J. Comput. Chem.* **32**, 1456 (2011).
- DFTB-oriented D3(BJ) parameters (DFTB+ / Grimme DFTB3 list):
  \(s_6=1\), \(s_8=0.5883\), \(a_1=0.5719\), \(a_2=3.6017\) Bohr.
- Molecular Grimme CN (not CP2K tanh); non-periodic.

## Geometry (water dimer)

Isolated (cubic box ≥ 20 Å, min-image unused for D3 non-periodic path):

| Atom | x [Å] | y [Å] | z [Å] |
| ---- | ----- | ----- | ----- |
| O | 0.000 | 0.000 | 0.000 |
| H | 0.757 | 0.000 | 0.586 |
| H | −0.757 | 0.000 | 0.586 |
| O | 0.000 | 0.000 | 2.900 |
| H | 0.757 | 0.000 | 3.486 |
| H | −0.757 | 0.000 | 3.486 |

## Acceptance (see acceptance.md)

1. Spec + REFERENCE + unit + regression green under SE-only cmake.
2. \(E_{\mathrm{disp}}\) matches offline s-dftd3 pin within stated tol.
3. \(E_{\mathrm{tot}}(d3) - E_{\mathrm{tot}}(off) = E_{\mathrm{disp}}\).
4. No public Config / factories; SE-only flags unchanged.

## CMake note

When `ATV_BUILD_SEMIEMPIRICAL` and not `ATV_BUILD_HF`, append
`src/cpu/gpw/internal/dftd3_params.cpp` and
`src/cpu/gpw/internal/dftd3.cpp` so the D3 engine links without turning HF on.
When HF is already on, do not double-list the sources.

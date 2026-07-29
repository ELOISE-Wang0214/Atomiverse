# TB / SE offline oracle pins

Family-locked **offline** references for the Atomiverse semiempirical /
tight-binding loop (`cpu/semiempirical`). **Never** launch OpenMOPAC, DFTB+,
or xtb from ctest — pin numbers and input twins live here; regressions hard-code
values from `REFERENCE.md`.

## Family matrix (binding)

| Family | Primary oracle | Forbidden as energy pin |
| ------ | -------------- | ----------------------- |
| NDDO (MNDO / AM1 / PM3 / RM1 / PM6 / PM7) | **OpenMOPAC** | DFTB+, xtb |
| DFTB0 / DFTB2 / DFTB3 / spin / PBC / TD-DFTB | **DFTB+** | OpenMOPAC, xtb |
| GFN-xTB (GFN0/1/2, IPEA, …) | **xtb** | OpenMOPAC, DFTB+; **not** a DFTB3 successor |
| Unified NDDO+TB API shape (no energy) | SCINE Sparrow | — (architecture only) |
| Evaluate / MD embed patterns | CP2K | — (architecture only) |

`REFERENCE.md` for every SE regression must state:

```text
oracle: mopac|dftb+|xtb|literature
oracle_version: <package version or formula citation>
```

## Layout

```text
tools/tb_oracle/
  README.md                 # this file
  inputs/
    dftbplus/               # DFTB+ HSD / geo.gen twins + notes
      ch4_3ob_dftb2/        # 3ob CH4 pure SCC-DFTB2 pin (wave-1)
      …                     # more twins land with worker tips / merges
    mopac/                  # (grow with NDDO OpenMOPAC pins)
    xtb/                    # (grow with GFN-xTB pins)
```

## Current offline twins (on `dev`)

| Path | Family | Notes |
| ---- | ------ | ----- |
| `inputs/dftbplus/ch4_3ob_dftb2/` | DFTB+ | CH4 3ob-3-1 under pure DFTB2 (`SCC=Yes`, no `ThirdOrder`) |
| `inputs/dftbplus/h2_mio_scc/` | DFTB+ | H₂ mio-1-1 absolute E pin twin (`reg_dftb2-dftbplus-energy-pin`) |
| `inputs/dftbplus/dftb2-periodic-gamma/` | DFTB+ | Γ-PBC C-chain geometry twin (ATV Hamiltonian pin) |
| `inputs/dftbplus/dftb2-periodic-forces/` | DFTB+ | Off-eq Γ-PBC C-chain force twin (analytic vs re-SCC FD) |
| `inputs/dftbplus/dftb2-kpoints/` | DFTB+ | 4×1×1 MP multi-k C-chain twin (`reg_dftb2-kpoints`; energy-only) |
| `inputs/dftbplus/dftb0_h2.note` | DFTB+ | H₂ non-SCC DFTB0 pin twin (`reg_dftb0-h2`; shares the mio-1-1 H₂ transcript) |
| `inputs/dftbplus/dftb3_h2o_mio_note.md` | DFTB+ | H₂O SCC-DFTB3 formula / module boundary (`reg_dftb3-water-mio`) |
| `inputs/dftbplus/dftb2_analytic_forces_note.md` | DFTB+ | Formula-level force module boundary (MioSketch) |
| `inputs/dftbplus/dftb2_spin_analytic_forces_note.md` | DFTB+ | Collinear open-shell force module boundary (O₂) |
| `inputs/dftbplus/dftb2_analytic_forces_threeob_note.md` | DFTB+ | ThreeOb closed-form F vs re-SCC FD (formula) |
| `inputs/dftbplus/dftb2_analytic_stress_note.md` | DFTB+ | Cell-stress module boundary (pair virial; formula) |
| `inputs/dftbplus/dftb2-dispersion-d3/` | s-dftd3 | Water-dimer **molecular** D3(BJ) geometry twin (`reg_dftb2-dispersion-d3`) |
| `inputs/mopac/` | OpenMOPAC | **absent on `dev`** — seed with AM1/RM1/PM* pins (Phase B) |
| `inputs/xtb/` | xtb | **absent on `dev`** — seed with GFN pins (Phase C, late) |

Every twin above is on `dev` as of harness wave 3 (2026-07-29; includes
DFTB3 note + kpoints/periodic-forces twins). GapCensus wave-4: still need
periodic-D3 / multi-k-force twins when those pending slugs staff; mopac/xtb
trees still empty.

## How to pin (operator / worker)

1. Run the family oracle **once** offline on the pin host (or cite a published
   formula when the binary is unavailable).
2. Record exact geometry, method flags, package version, and numbers in
   `regressions/semiempirical/<slug>/REFERENCE.md`.
3. Copy minimal input twins under `tools/tb_oracle/inputs/<family>/…` so the
   pin is reproducible without network or live oracle in CI.
4. Hard-code tolerances in `reg_<slug>.cpp` — ctest must never shell out to
   mopac / dftb+ / xtb.

## Temps

Scratch only: `/tmp/tb_oracle_*` `/tmp/mopac_*` `/tmp/dftbplus_*` `/tmp/xtb_*`
`/tmp/atv-*`. Wipe after pin generation; never commit oracle binaries.

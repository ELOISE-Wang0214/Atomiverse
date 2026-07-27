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
| `inputs/dftbplus/dftb2_analytic_forces_note.md` | DFTB+ | Formula-level force module boundary (no host binary) |

Additional twins may exist only on unmerged `wf/*` tips (e.g. H2 mio SCC,
DFTB0 notes, D3 dimer XYZ, Γ-point PBC). Merge brings them onto `dev`.

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

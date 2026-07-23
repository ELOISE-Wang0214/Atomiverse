# CP2K offline oracle (CP2K capability loop)

Mirrors `tools/pyscf_oracle/`: produce **hard-coded pins** for regressions.
**ctest never launches CP2K.**

## Contract

1. Human or agent runs CP2K **offline** with a fixed input under `inputs/`.
2. Record energy / forces / stress into `REFERENCE.md` under the regression
   (cite CP2K version, functional, cutoff, basis/PP, cell).
3. Regression binary asserts Atomiverse vs those numbers within stated tol.
4. Optional fair wall: same geometry; write \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\)
   into REFERENCE (not ctest).

## Twin used by `gpw-forces-water`

- Input: `inputs/h2o_forces_mt.inp` (isolated H2O, PBE, DZVP-MOLOPT-GTH, MT, 400 Ry)
- Pin destination: `regressions/gpw/gpw-forces-water/REFERENCE.md`
- CP2K data: `$CP2K_DATA_DIR` or EasyBuild share path `…/share/cp2k/data`

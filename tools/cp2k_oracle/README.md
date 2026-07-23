# CP2K offline oracle (CP2K capability loop)

Mirrors `tools/pyscf_oracle/`: produce **hard-coded pins** for regressions.
**ctest never launches CP2K.**

## Contract

1. Human or agent runs CP2K **offline** (login/interactive/batch) with a fixed input.
2. Record energy / forces / stress / frequencies into `REFERENCE.md` under the
   regression (cite CP2K version, functional, cutoff, basis/PP, cell, k-mesh).
3. Regression binary asserts Atomiverse vs those numbers within stated tol.
4. Optional fair wall: same geometry, comparable settings; write
   \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) into REFERENCE (not ctest).

## Suggested layout (grow as needed)

```text
tools/cp2k_oracle/
  README.md          ← this file
  inputs/            ← minimal CP2K .inp twins (optional, not required at ctest)
  pin_notes.md       ← how pins were taken
  fair_bench.py      ← optional wall twin (later)
```

## Version pin

Prefer a single CP2K release string in REFERENCE (e.g. `CP2K 2024.1` or module
version on Arrhenius). Do not mix versions across a regression without noting it.

## Relation to tutorials

- GPW narrative: `docs/tutorial/ch09-gpw.md` (Quickstep parity story).
- Loop product pages: `docs/tutorial/methods/<slug>.md` from `cp2k-method-loop`.

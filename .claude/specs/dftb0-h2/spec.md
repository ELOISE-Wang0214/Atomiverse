# Spec — `dftb0-h2`

## Goal

Ship **non-SCC DFTB / DFTB0** as a concrete `cpu::semiempirical::DFTB0`
`Semiempirical` component that validates the SK → H/S → generalized eigen →
band + repulsive pipeline **without** SCC charge iteration. Foothold system:
**H2**. Family oracle: **DFTB+** (not OpenMOPAC, not xtb).

## Why

Ladder step 1 of the TB method loop. DFTB2 already ships SCC; DFTB0 proves the
non-SCC electronic stack and is the correct baseline before analytic forces /
true mio pin work.

## API (public)

```cpp
cpu::semiempirical::DFTB0 dftb0;  // members only — no *Config, no make_*
dftb0.with_forces = true;        // optional; default true
// Evaluate path
eval.on(evaluate::Phase::evaluate).add(std::move(dftb0));
```

- Parameters: public members (`with_forces`). No public `DFTB0Config`.
- Static SK/repulsive tables: existing `basis_data/dftb_mio_tables.h` (H/C/N/O
  mio-style sketch — same tables as DFTB2; not bit-identical to DFTB.org
  mio-1-1).
- Units: energy eV, positions Å, forces eV/Å on Frame.

## Physics

Non-SCC DFTB (DFTB0 / “DFTB1” electronic):

1. Onsite H0 / unit S from element table
2. Two-centre SK H0 / S from pair table (same radial forms as DFTB2)
3. Generalized eigen of H0 C = S C ε
4. Closed-shell occupation: `n_occ = n_elec/2`
5. `E = Σ_i 2 ε_i + E_rep` with `E_rep = Σ_{A<B} A exp(-B R)`
6. Forces: central FD of total E (skeleton; analytic SK gradient is later)

No Hubbard, no γ, no Δq, no SCC residual.

## Regression

`regressions/semiempirical/dftb0-h2/`

- H2 at R = 0.74 Å (and a second R for non-flatness)
- Hard-coded pins in `REFERENCE.md` with `oracle: dftb+` and `oracle_version`
- ATV asserts vs REFERENCE within stated tol via Evaluate + DFTB0
- SE-only cmake (no MD/MINIMIZER/HF/CUDA)

## Out of scope

- SCC (DFTB2), DFTB3, spin, PBC, analytic gradient product depth
- True DFTB.org mio-1-1 dump / DFTB+ bit-identity energy (→ `dftb2-dftbplus-energy-pin`)
- OpenMOPAC or xtb as reference

## Acceptance

See `acceptance.md`.

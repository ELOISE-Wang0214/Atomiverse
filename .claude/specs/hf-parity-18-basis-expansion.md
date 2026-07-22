---
title: Static-table basis expansion (6-311G*, cc-pVTZ/QZ, def2-TZVPP) + g angular momentum
slug: hf-parity-18-basis-expansion
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P4. Run `/mol:spec
  hf-parity-18-basis-expansion` to promote. NO runtime parsers (permanent
  non-goal) — static C++ tables only, transcribed once.
---

# Basis breadth + g angular momentum

## Summary
Expand the static basis tables: **6-311G\***, **cc-pVTZ/QZ**, **def2-TZVPP**,
heavier def2 rows — still **no runtime parsers**. Raise orbital-basis angular
momentum to **g** (`kMaxL` is f today) if a regression needs it. Parent:
hf-deferred-backlog P4.

## Domain basis
Bases are static C++ data headers transcribed from published tables (the GPW
tables set the precedent). g functions need the l=4 (g) Cartesian/spherical
handling through the integral + derivative + c2s paths (`kMaxL` currently 4 for
DF aux in places, f for orbitals). Reference: EMSL/BSE published exponents;
Weigend/Ahlrichs def2.

## Design (direction)
- Add basis rows to the static tables (constructor-lookup only, throw on
  unsupported). Extend integral/derivative/c2s coverage to g if needed.

## Files (best-effort)
- `include/atomiverse/cpu/scf/basis_data/` tables, `basis_set.h` enum
- `cart.h` (`kMaxL`), integral + derivative + cart2sph paths for g
- `regressions/scf/<basis>_energy/` (new)

## Tasks (provisional)
- [ ] Add the new basis rows (static tables) + enum values.
- [ ] (If needed) g-shell integral/derivative/c2s coverage.
- [ ] Regressions: energies in the new bases vs pinned PySCF.

## Testing
Unit: table lookup + a spherical-overlap pin per new basis (seconds). Regression:
one energy per major new basis ≤ 60 s.

## Out of scope
Runtime basis/ECP file parsers (permanent non-goal). CFMM-K (non-goal).

## Depends on
Independent; g-shell work reuses the shipped spherical c2s path.

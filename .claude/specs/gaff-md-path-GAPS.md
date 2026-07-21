# Capability gaps — GAFF / flexible PEO MD path

Inventory for the `gaff-md-path` chain. Target experiment:
molpack-paper `figures/flexible` (neat PEO-25, GAFF2/AM1-BCC, N=100,
Gudla–Zhang pre-production protocol). Specs close the **engine path**, not
paper-length multi-replica production.

| # | Gap | Severity | Spec |
|---|-----|----------|------|
| 1 | `load_lammps_data` is `atom_style molecular` only — no charges for GAFF `full` data | **Blocker** | 01-io-full |
| 2 | `PairLJCut` / `PairCoulCut` specials only hard-exclude scale=0; no Amber 1–4 ×0.5 | **Blocker** | 02-pair-specials |
| 3 | `CoulLong` has **no** special_bonds — bonded neighbors get full real-space Coulomb | **Blocker** | 03-coul-long-specials |
| 4 | kcal/mol (LAMMPS real) ↔ eV scattered; Fourier dihedral (single-term) not mapped in library helpers | Medium | 04-ff-compose |
| 5 | No GAFF energy/force cross-check vs reference pins | **Scientific gate** | 05-xcheck |
| 6 | No staged min→NVT→NPT smoke on charged molecular topology | Integration | 06-protocol-smoke |
| 7 | Roadmap text understates multi-type LJ / overstates missing pieces | Docs | 07-docs |

## Already available (not gaps)

- MD driver NVE/NVT/NPT + temperature ramp (`NoseHooverNPT`)
- Harmonic bond/angle multi-type vectors
- `DihedralCharmm` (covers PEO single-term `dihedral_style fourier`)
- `PairLJCut` multi-type + Lorentz–Berthelot
- `CoulLong` + `KSpacePME` (ions; NaCl regression)
- Minimizer LBFGS/CG
- Improper harmonic (unused by neat PEO-25)

## Explicitly deferred (after this chain)

- Multi-term `dihedral_style fourier` component
- Full `peo25.ff` text / label-based parser
- CUDA Amber specials parity
- 30 ns × multi-replica paper parity (density / \(R_g\))
- Soft-core pair styles / Packmol integration

## Start here

```text
/mol:impl gaff-md-path-01-io-full
```

Worktree: `.claude/worktrees/gaff-md-path` on branch `feat/gaff-md-path`
(based on clean `dev` HEAD; does not include uncommitted longrange WIP from
the primary checkout).

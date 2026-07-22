---
title: MOM/IMOM/non-Aufbau ΔSCF + fractional occupation / Fermi smearing
slug: hf-parity-13-delta-scf-occupation
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P2. Run `/mol:spec
  hf-parity-13-delta-scf-occupation` to research + promote. science.required.
---

# ΔSCF occupation control (MOM/IMOM) + fractional occupation

## Summary
Add the **MOM / IMOM / non-Aufbau ΔSCF** product path on
`SlaterDeterminant::with_occupation` (core-hole, state following), plus
**fractional occupation / Fermi smearing** for near-degeneracy. Parent:
hf-deferred-backlog P2. `regressions/scf/n2_ionization` already exercises a
ΔSCF-adjacent case.

## Domain basis
MOM: at each SCF step select occupied MOs by maximum overlap with the previous
step (not by energy), enabling excited-determinant / core-hole states. Fermi
smearing: fractional occupations f_i = 1/(1+e^{(ε_i−μ)/kT}) with μ from electron
count. Reference: Gilbert/Besley/Gill MOM; Mermin finite-T HF.

## Design (direction)
- Occupation-selection policy plugged into the SCF density build (overlap-based
  for MOM; smearing function for fractional). Reuse `with_occupation`.

## Files (best-effort)
- `src/cpu/scf/` SCF density/occupation path, `slater_determinant`
- `regressions/scf/{core_hole_mom,fermi_smearing}/` (new)

## Tasks (provisional)
- [ ] RED: MOM follows a chosen non-Aufbau state on a tiny system (unit).
- [ ] MOM/IMOM overlap-based occupation selection.
- [ ] Fractional occupation / Fermi smearing.
- [ ] Regressions: core-hole ΔSCF + a smeared near-degenerate case vs pinned PySCF.

## Testing
Unit: occupation selection (seconds). Regressions vs pinned PySCF ≤ 60 s.

## Out of scope
Ensemble DFT; state-averaged CASSCF.

## Depends on
`SlaterDeterminant::with_occupation` (in tree).

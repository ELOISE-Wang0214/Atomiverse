---
title: Retire transitional debt (SCFSolver / SCFConfig / IntegralPack / AtomList)
slug: hf-parity-21-retire-scf-debt
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P5. Architectural — likely /mol:refactor
  territory, big blast radius. Run `/mol:spec hf-parity-21-retire-scf-debt` to
  scope carefully first. Do LAST (many links use IntegralPack internally today).
---

# Retire transitional SCF debt

## Summary
Retire the transitional public surface — `SCFSolver` / `SCFConfig` / public
`IntegralPack` / `AtomList` / `NDDOSystem` — toward the target: concrete methods
only (`RHF`/`UHF`/`RKS`/`UKS`/AM1/PW/GPW), Frame as sole I/O, integrals as
private workspace. Parent: hf-deferred-backlog P5 + the CLAUDE.md debt note.

## Domain basis
CLAUDE.md target shape: a Driver is a staged loop, a method is concrete,
parameters are members, `IntegralPack`/AO maps are private workspace rebuilt
from Frame. The debt table lists exactly what must go.

## Design (direction — refactor, staged)
- Move `IntegralPack` etc. behind the concrete methods (private); delete the
  public bags and factory-shaped entry points. Grep-gated by
  `no_factory_scrub`. Behavior-preserving under the full test gate; likely
  several `/mol:refactor` passes, not one impl.

## Files (best-effort)
- `include/atomiverse/cpu/scf/` public headers, `methods.h`, callers across the tree
- every regression/test that currently constructs `IntegralPack`/`SCFSolver`

## Tasks (provisional)
- [ ] Inventory every public use of the debt types (grep).
- [ ] Move integrals/AO maps to private method workspace.
- [ ] Delete `SCFSolver`/`SCFConfig`/public `IntegralPack`/`AtomList`/`NDDOSystem`.
- [ ] Green full suite + architecture scrubs after each staged pass.

## Testing
Behavior-preserving: full ctest + `no_factory_scrub` + `test_scf_phase_scrub`
green throughout. No new e2e.

## Out of scope
New features (this is pure debt retirement). New public config bags (banned).

## Depends on
Do LAST — every earlier link that touches SCF internals should land first so this
refactor rewrites a stable surface once.

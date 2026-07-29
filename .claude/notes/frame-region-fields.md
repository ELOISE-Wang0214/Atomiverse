# Frame region / selection fields

**Status:** binding (2026-07-29). Router: `CLAUDE.md` iron law
`frame-region-fields`.

## Rule

Membership and partition (**which** atoms / bonds / region) live on the
**Frame**. Methods **read** at `build`/`launch`; they do not own parallel
region bags or primary `set_links` / `set_mask` / index-list APIs.

Missing or invalid **required** fields throw with a how-to. No silent
whole-system default that hides a forgotten mask.

## What stays on the object

- Scalar / table **parameters** (cutoffs, tolerances, link **scale** if any)
- **Slots** for other objects: `set_evaluator`, `set_engine`, `set_integrator`, …

## Derived vs declared

| Declared on Frame (user) | Derived (engine) |
|--------------------------|------------------|
| QM mask `kAtomQmZ` (`atoms.qm_z`) | QM/MM link cuts = bonds with mixed QM/MM endpoints |
| Future: freeze / electrode / shell columns as introduced | Any pure function of mask + topology + geometry |

**No bonds ⇒ no link cuts.** Valid for LJ / nonbonded-only QM/MM.

## Forbidden user paths

- `links()` return-ref assignment, `set_links`, ctor `vector` of region indices
  as the main partition story
- Parallel `*System` / side tables that duplicate Frame membership

## QM/MM (illustrative)

User: fill `kAtomQmZ` (+ normal topology when covalent).  
Engine: optional auto link H from cross-boundary bonds; three connection
evaluators compose under one `MD` via `set_evaluator`.

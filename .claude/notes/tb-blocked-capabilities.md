# TB / SE blocked capabilities — why they are not pending (2026-07-29)

Durable harness for rows that are **not ordinary pending**. A row leaves a hold
list only when its gate is *actually* removed.

## Permanent rules (never capabilities to implement)

| Id | Rule |
| -- | ---- |
| `mopac-as-dftb-oracle` | Forbidden: pin DFTB to OpenMOPAC |
| `gfn-ff` | Scope-out: classical GFN-FF → MD/FF surface |
| `om3-nddo` | No family-locked pin: OM3 ∉ OpenMOPAC primary surface. Expanding to Thiel OMx is an oracle-set decision, not a missing ATV feature gate. |

## Closed from blocked (2026-07-29)

| Slug | Resolution |
| ---- | ---------- |
| `se-analytic-hessian` | Shipped `reg_se-analytic-hessian` (VIB+DFTB2 H₂) |
| `dftb-plus-u` | Design note + energy-only Dudarev on `DFTB2::u_eff` + `reg_dftb-plus-u` |
| `dftb-gamma-units-fix` | Fixed γ units (was pending head, not blocked) |

## Still blocked (science multi-wave)

| Slug | Gate | Unblock |
| ---- | ---- | ------- |
| `noncollinear-spin-dftb` | Spinor Hamiltonian (2-component density), not a parameter on collinear | Spec + multi-wave after collinear solid |
| `spin-orbit-dftb` | SOC/SK design | After noncollinear |

## Phase C dependents (pending, not blocked)

`stda-xtb` / `gfn-xtb-analytic-hessian` wait on production GFN ground state
(`gfn0-beyond-toy` → `gfn2-…`). Staff under Phase C; do not re-list as blocked.

Related: [[dftb-plus-u]], [[dftb-gamma-units-defect]].

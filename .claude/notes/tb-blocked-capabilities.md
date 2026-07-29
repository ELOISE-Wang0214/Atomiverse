# TB / SE blocked capabilities — why they are not pending (2026-07-29)

Durable harness copy of the `## blocked` set in
`docs/developer/tb-method-queue.md`. The queue is a worklist that census
rewrites every wave; **this file is the reason each row is blocked**, so an
agent does not re-staff it, re-discover the reason, or quietly promote it to
pending.

A row leaves `blocked` only when its stated gate is *actually* removed — not
because a wave was short of work.

## Blocked on a missing predecessor (unblock by landing the predecessor)

| Slug | Gate |
| ---- | ---- |
| `noncollinear-spin-dftb` | Needs collinear spin-DFTB as the base. The **collinear** foothold now ships (dual α/β, Mulliken m, analytic F and σ), but noncollinear needs a spinor Hamiltonian — a different data layout, not a parameter. Still blocked. |
| `dftb-plus-u` | Needs the SCC base (ships) **plus a Hubbard-U design note** that does not exist. Do not improvise the +U double-counting correction. |
| `spin-orbit-dftb` | Needs a specialized SOC / SK design. Not a ladder head; schedule only after noncollinear spin. |
| `stda-xtb` | Needs a **production** GFN ground state. `XtbGfn0` is an H/C/O toy, so sTDA has nothing to sit on. |

## Blocked by the SE-only build (unblock by changing the loop's cmake)

| Slug | Gate |
| ---- | ---- |
| `se-analytic-hessian` | Needs `ATV_BUILD_VIB` / `cpu::vib::Hessian`. The TB loop's binding cmake keeps **VIB OFF**, so the capability cannot even link in this loop. |
| `gfn-xtb-analytic-hessian` | Same VIB gate. xtb's semi-numeric Hessian is product-later regardless. |

These two are blocked by a **loop policy**, not by physics. If the operator
widens the SE-only cmake, both become pending immediately — but widening it
also re-admits the HF/PW/GPW build cost the loop exists to avoid.

## Blocked by the oracle rule (do not "solve" by switching oracles)

| Slug | Gate |
| ---- | ---- |
| `om3-nddo` | OM3 is **not** in OpenMOPAC's primary surface. The NDDO family is oracle-locked to OpenMOPAC, so there is no family-legal pin. Thiel's OMx codes are outside this loop's oracle set. |
| `mopac-as-dftb-oracle` | **Forbidden pattern**, never pending. Pinning a DFTB number to OpenMOPAC crosses families: different Hamiltonian, different parameterization. Listed in `blocked` precisely so it stays visible as a rule rather than being re-invented as a shortcut. |

## Blocked as out of scope (never pending)

| Slug | Gate |
| ---- | ---- |
| `gfn-ff` | Classical force field. The TB loop is the **electronic** SE/TB loop; GFN-FF belongs to the classical MD surface. Scope-out, not "later". |

## The rule these encode

**Family-locked oracles** (also in CLAUDE.md): NDDO → OpenMOPAC; DFTB0/2/3 /
spin / PBC / TD-DFTB → DFTB+; GFN-xTB → xtb. A capability whose family oracle
cannot produce the number is **blocked**, not re-pinned against a different
family and not silently downgraded to `oracle: literature`. Where no oracle
binary and no published value exist, a case may ship **structural / formula
identities** (analytic-vs-FD, extensivity, symmetry limits) — and must say so
in its REFERENCE `Gaps`. See [[oracle-pin-provenance]] for what a pin is
allowed to be.

Related: [[dftb-gamma-units-defect]] — not blocked, but it makes every polar
MioSketch pin provisional until it lands.

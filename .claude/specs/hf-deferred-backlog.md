---
title: "HF deferred backlog — residual Jensen / PySCF gaps"
slug: hf-deferred-backlog
status: approved
created: 2026-07-21
chain: none
note: >
  Parking lot only. Do NOT run /mol:impl on this umbrella.
  When picking work, split one theme into a fresh child spec chain.
supersedes: >
  Closed 2026-07-21 (operator manual attestation after review):
  hf-capability-gaps-01…21, scf-driver-refactor-06-criterion-unify.
  Shipped code remains; unfinished science/docs/hard-gates live below.
---

# HF deferred backlog — residual Jensen / PySCF gaps

> **Decomposed 2026-07-22:** this inventory is now expanded into the 21-link
> **`hf-parity` chain** (`hf-parity-01-ecp-forces` … `-21-retire-scf-debt`, all
> `[draft]` — see INDEX.md). This file stays the priority map; the chain holds
> the per-theme scope/design. When a chain link ships, strike its matching row
> below (as the P0 spherical-gradient bullet already is). Do **not** `/mol:impl`
> this umbrella — work the chain links.

## Summary

After the 2026 HF capability campaign (`hf-capability-gaps-01`–`21`) and the
SCF criterion unify refactor, the **teaching-core HF surface is in tree**:
RHF/UHF/ROHF/GHF (real), basis tables through aug-cc-pVDZ + def2-ECP energy,
guess/DIIS/SOSCF/stability, multipoles/charges/Koopmans, RI-JK, CPHF Hessian
(RHF), TDHF/TDA (RHF), embed/C-PCM/x2c-1e/CP energy paths, FMM kernel + CFMM-J
hook, longrange rename.

This file is the **single residual inventory** vs Jensen (*Introduction to
Computational Chemistry*) and PySCF everyday HF tooling. It is **not** an
implementation plan and **must not** be driven by `/mol:impl` or
`/mol:impl-all`. When free time appears, **carve one theme into a new
`/mol:spec` chain**, implement that chain, and strike the matching rows here
(or supersede this inventory).

## Domain basis

- Jensen textbook arc: HF equations → basis → open shell → properties →
  response / vibrations → environment (solvent, field, relativity) → scaling.
- PySCF user defaults: full analytic forces, DF + screening, PCM family,
  response (polar / TD), broad basis/ECP libraries, post-HF stack, optional
  PBC and GPU paths.
- Engine invariants (unchanged by this backlog): Frame I/O only; no runtime
  basis/ECP file parsers; no factory soup; concrete methods `RHF`/`UHF`/…;
  unit tests never e2e; regressions ≤ 60 s with cited literature / pinned
  offline PySCF values.

## How to use this inventory

1. Pick **one** priority block (P0–P3) or a single bullet.
2. `/mol:spec <theme>` → real design + acceptance; link back here in Summary.
3. `/mol:impl` that child only. On close, edit this file: mark the bullet
   `[x]` or delete it; if the inventory empties, `/mol:close hf-deferred-backlog --manual`.
4. Never treat remaining rows as a single mega-impl.

## Residual checklist

Priority is product / textbook impact, not historical chain order.

### P0 — Forces and geometry (largest user gap)

- [x] Analytic gradients for **spherical** bases (shipped 2026-07-22, child spec
      `scf-spherical-gradients`: RHF def2-family forces via the c2s sandwich on
      the derivative integrals; H2O def2-SVP regression vs pinned PySCF 1e-6
      Ha/Bohr). UHF/UKS spherical forces follow the same path if needed.
- [x] **ECP** nuclear derivatives (shipped 2026-07-22, `hf-parity-01-ecp-forces`:
      **semi-numerical** dU/dR — FD of the cheap ECP 1e block, everything else
      analytic; I₂ def2-SVP+ECP matches PySCF analytic gradient to 1e-9 Ha/Bohr).
      Full analytic dU/dR recurrences remain a future upgrade.
- [ ] Analytic forces for **ROHF** and **GHF** (explicit throw / out of scope
      in the closed campaign).
- [ ] Forces through environment members that currently throw:
      **C-PCM**, **x2c**, **point-charge embed / external field**, **RI-JK**,
      **CFMM-J**.
- [ ] Validate FD / filtered-subspace nuclear gradients where still unpinned.
- [ ] IR intensities (dipole derivatives) on top of CPHF Hessian.

### P1 — Hard gates left open when 01–21 were closed

Code is present; science / docs ledgers were not fully verified. Re-open as
tight child specs only if values still fail or docs are dishonest.

- [ ] **CFMM-J hard gate**: elementwise vs Schwarz direct J; converged RHF
      energy identity (~1e-7 Ha); alkane crossover table; ch19 honest
      “PySCF has no molecular FMM” statement.
- [ ] **x2c-1e**: pVp vs pinned PySCF `int1e_pnucp`; decoupling identities;
      Ar (or similar) energy vs `x2c1e` (~1e-5 Ha); ch20 relativity section
      honesty (no picture-change on properties).
- [ ] **C-PCM**: solvation free-energy vs pinned PySCF (~0.1 kcal/mol class);
      document non-electrostatic absence + cavity discontinuity.
- [ ] **Embed / field**: finite-field dipole ≡ analytic dipole; point-charge
      H₂O vs pinned `qmmm.mm_charge` bookkeeping (MM–MM).
- [ ] Any remaining **docs/Doxygen** stubs from the closed chain that still
      claim Handler language, old `SCFSolver` user path, or missing Results
      tables.

### P2 — Response, occupation, open-shell spectroscopy

- [ ] Static **polarizability** (electric-field CPHF); optional hyperpolarizability.
- [ ] **UHF** CPHF / analytic Hessian; RKS/UKS second derivatives (XC kernel).
- [ ] **UHF-TDHF** / spin-flip TDA; **TDDFT** (RKS/UKS A/B).
- [ ] Excited-state gradients / NACs (only if a concrete consumer needs them).
- [ ] **MOM / IMOM / non-Aufbau** ΔSCF product path on
      `SlaterDeterminant::with_occupation` (core-hole, state following).
- [ ] Fractional occupation / Fermi smearing for near-degeneracy.

### P3 — Solvent / relativity / analysis depth

- [ ] Beyond C-PCM: IEF-PCM / SS(V)PE / SMD-class non-electrostatics; optional
      **ddCOSMO** linear-scaling solver.
- [ ] Surface switching (SWIG / York–Karplus) if geometry+PCM becomes common.
- [ ] Polarizable / periodic MM embedding; Ewald-corrected external charges.
- [ ] x2c **picture-change** property operators; optional DKH; SOC / complex
      GHF only if a real application demands it (not “parity for parity”).
- [ ] SO-ECP; lanthanide 4f-in-core tables if heavy-element work appears.
- [ ] Population depth: NPA/NBO, IAO/IBO, RESP (Mulliken/Löwdin/MK/CHELPG exist).

### P4 — Basis breadth and large-system HF

- [ ] Static-table basis expansion: 6-311G*, cc-pVTZ/QZ, def2-TZVPP, heavier
      def2 rows (still **no** runtime parsers).
- [ ] Angular momentum **g+** if a regression needs it (`kMaxL` still f for
      orbital bases; DF aux already to l=4 in places).
- [ ] **CFMM-K** is intentionally out forever unless a new design appears —
      keep K on direct/RI; do not pretend CFMM solves exchange.
- [ ] Molecular **k-point / PBC HF** — only if periodic mean-field on Gaussians
      becomes a product goal (PW/GPW already cover periodic DFT).

### P5 — Correlation and backends (adjacent to HF)

- [ ] RI-MP2 product path / DF-ao2mo completeness (MP2 teaching path exists).
- [ ] CCSD / CCSD(T) / CASSCF — only under a dedicated science campaign, not
      as HF parity.
- [ ] **CUDA SCF gradients** and any HF capability that must mirror CPU
      (today CUDA SCF is energy/density-oriented).
- [ ] Retire transitional debt: public `SCFSolver` / `SCFConfig` /
      `IntegralPack` / `AtomList` user surface (target: concrete methods only).

### Permanent non-goals (do not re-open as “gaps”)

- Runtime basis or ECP file parsers.
- Turning `SCFSolver` into a method-enum bag.
- E2E full SCF inside `tests/` / ctest.
- Exact one-to-one PySCF feature parity as a success metric.
- x2c + ECP double counting (keep mutual exclusion).

## Design (when a child is carved)

- Child specs own Files / Tasks / Testing; this inventory only prioritizes.
- Prefer extending concrete method members + private workspace; no new
  public Config bags.
- Every delivered force or property path needs one `regressions/scf/...`
  textbook case with literature or offline-pinned PySCF provenance
  (appendix protocol).
- Keep unit tests second-scale; no external engine at build/test time.

## Files

- This inventory: `.claude/specs/hf-deferred-backlog.md` (+ acceptance sibling).
- Child work lands under `include/atomiverse/cpu/scf/`, `src/cpu/scf/`,
  `tests/cpu/scf/`, `regressions/scf/`, `docs/tutorial/ch*.md` as usual.

## Tasks

- [ ] Keep this inventory current when HF work lands elsewhere (edit in place).
- [ ] When starting work: spawn a **child** spec; do not implement under this slug.
- [ ] When the checklist is empty or abandoned: close this slug with `--manual`.

## Testing

No new tests for the inventory itself. Child specs define their own unit +
regression gates.

## Out of scope

- Implementing any checklist row under this slug.
- Re-opening the closed `hf-capability-gaps-*` files (history is the git log).
- Non-HF domains (classical MD, PW/GPW feature work) except where listed as
  explicit HF-adjacent backends.

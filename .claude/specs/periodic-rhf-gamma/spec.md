---
slug: periodic-rhf-gamma
status: approved
created: 2026-07-24
revised: 2026-07-24
domain: scf / hartree-fock
queue: pyscf-engine-gaps
deferred: true
depends_on: []
unblocks: [periodic-rhf-gamma]
grilled: false
sibling: rsh-eri-cam-b3lyp
---

# Spec: `periodic-rhf-gamma` — Γ-point periodic RHF foothold

**Slug:** `periodic-rhf-gamma`  
**Domain:** Gaussian-basis Hartree–Fock under PBC (CPU)  
**Queue:** PySCF engine gaps (deferred — not auto-impl)  
**Status:** approved, **implement later**  
**Unblocks queue row:** `periodic-rhf-gamma`  
**Sibling:** `rsh-eri-cam-b3lyp` (independent; no ordering dep)

## Summary

Deliver **Γ-point periodic restricted Hartree–Fock** energy as a product
foothold: Frame in, energy out, offline pin against PySCF `pbc.scf`. Today
molecular RHF/RKS are solid, and GPW already does **Γ-point lattice sums for
one-electron** operators (S, T, local PP), but **four-index exact exchange
under PBC is missing** — PW/GPW product paths are semi-local DFT only
(`docs/roadmap/hartree-fock.md` §2). This spec adds the minimal periodic HF-X
+ closed-shell SCF so the blocked queue row can graduate; k-point sampling,
UHF, and analytic forces are follow-ons.

## Domain basis

**Γ-point periodic orbitals** (real Bloch phase \(e^{i\mathbf{k}\cdot\mathbf{R}}=1\)):

\[
\chi_\mu^\mathrm{PBC}(\mathbf{r})
  = \sum_{\mathbf{T}}
    \varphi_\mu(\mathbf{r}-\mathbf{R}_\mu-\mathbf{T}),
\]

with lattice vectors \(\mathbf{T}\) of the simulation cell. One-electron
matrix elements become lattice sums already present in GPW:

\[
S_{\mu\nu}
  = \sum_{\mathbf{T}}
    \langle \varphi_\mu | \varphi_\nu(\cdot-\mathbf{T}) \rangle.
\]

**Exact exchange** at Γ (real AOs) needs a lattice-summed four-center
Coulomb operator (naive form; production engines use cutoff / CFMM /
Ewald-split exchange — v1 may use a **truncated lattice sum** with
documented cutoff policy matching the PySCF twin settings):

\[
K_{\mu\nu}
  = \sum_{\lambda\sigma} P_{\lambda\sigma}
    \sum_{\mathbf{T}}
    (\mu\lambda | \nu^\mathbf{T}\sigma^\mathbf{T})
    \quad\text{(schematic; image placement must be consistent)}.
\]

Total RHF energy (closed shell) in Hartree:

\[
E = \mathrm{Tr}(P H_\mathrm{core})
  + \tfrac12 \mathrm{Tr}\big(P(J-K)\big)
  + E_\mathrm{nn}^\mathrm{PBC}.
\]

Nuclear repulsion under PBC must use the **same electrostatic convention**
as the electronic twin (Ewald / minimum-image policy as pinned).

References:

- McMurchie–Davidson / Rys molecular ERI (existing engine) as the
  per-image building block.
- PySCF PBC HF: Sun et al., *WIREs Comput. Mol. Sci.* (2018);
  use PySCF 2.14.x `pbc.gto` + `pbc.scf.RHF` as **offline oracle only**.
- Roadmap: `docs/roadmap/hartree-fock.md` §2; GPW lattice helpers
  `cpu/gpw/internal/periodic.h`, `one_electron.h`.

Units: internal Ha / Bohr; Frame `kEnergy` in eV.

## Design

### Entities

| Piece | Location | Role |
| ----- | -------- | ---- |
| Lattice image enum | `cpu/gpw/internal/periodic.h` (pattern) | reuse / extract shared Γ image ranges |
| Molecular ERI | `cpu/scf/rys_eri.*` | per-image (μν\|λσ) building block |
| Periodic HF exchange | `cpu/scf/` new internal module | lattice-summed K (and J if not already PBC-aware) |
| Periodic one-electron | scf or shared with GPW | S, T, V_nuc / PP under Γ PBC |
| Concrete method | e.g. `cpu::scf::RHF` members **or** thin `PeriodicRHF` | Frame in/out; **no Config bag** |
| SCF loop | existing RHF Driver-shaped solve | reuse DIIS / density path once Fock is PBC-correct |

### Public API (binding shape — v1)

Prefer **members on concrete RHF**, not a second world-model type:

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::STO3G;   // or GTH/MOLOPT only if required by twin
rhf.periodic = true;                    // public member; default false
// Frame.box provides orthorhombic (v1) cell [Å]; positions in cell
rhf.tol = 1e-8;
rhf.max_iter = 100;
rhf.run(frame);                         // frame.meta.get(kEnergy) [eV]
// optional diagnostic: last energy [Ha] via existing last_result style if present
```

Alternative (only if RHF molecular path must stay branch-free): a concrete
`cpu::scf::PeriodicRHF` with the same member style — still **no** `*Config`,
**no** `make_periodic_rhf`.

v1 constraints (binding):

- **Γ-point only**, closed-shell **RHF**, **orthorhombic** cell.
- Truncation / Ewald policy for exchange **documented in REFERENCE.md** and
  mirrored in the PySCF pin input (same cell, basis, cutoff).
- Single-rank CPU; no MPI domain split required for the foothold.

### Reuse decision

| Candidate | Tag | Decision |
| --------- | --- | -------- |
| `cpu/gpw/internal/periodic.h` lattice image counts | reuse / generalize | **reuse** image enumeration pattern; **generalize** into a shared scf/gpw helper only if both call sites need the same header — do not duplicate three copies of image loops |
| GPW `one_electron` Γ lattice S/T | pattern | **pattern** — periodic HF one-electron may live under `cpu/scf/internal/` but follow the same Σ_T contraction style; do not call into GPW DFT solver |
| `rys_eri_block` | reuse | **reuse** as molecular quartet kernel inside the image loop (or pre-shifted shell centers) |
| Molecular `RHF` / `hamiltonian_builder` / DIIS | reuse / generalize | **generalize** Fock assembly for PBC K (and J); keep SCF iteration, DIIS, density build |
| CFMM J (`cfmm_j`) | pattern | optional later acceleration; **v1** may use direct lattice-sum J for tiny cells |
| PySCF `tools/pyscf_oracle` | reuse | **reuse** offline pin workflow |

## Files to create or modify

- `include/atomiverse/cpu/scf/rhf.h` (or `periodic_rhf.h` if split) — public members
- `src/cpu/scf/...` SCF / Fock path for PBC
- `include/atomiverse/cpu/scf/internal/` + `src/cpu/scf/internal/` — periodic exchange (new)
- Possibly thin shared lattice helper if extracted from GPW (document path)
- `tests/cpu/scf/test_periodic_rhf_gamma_energy.cpp` or split unit tests for K assembly
- `tests/cpu/scf/test_periodic_exchange_lattice.cpp` (new) — image sum sanity
- `regressions/scf/periodic-rhf-gamma/` (new) + `REFERENCE.md`
- `tools/pyscf_oracle/` pin notes / script
- `docs/roadmap/hartree-fock.md` — Γ-point periodic RHF ○→✓
- `docs/developer/pyscf-method-queue.md` — blocked → done

## Tasks

- [ ] Write failing unit tests for Γ lattice-summed exchange on a 1–2 electron toy (finite image count, hand or offline pin) (`tests/cpu/scf/test_periodic_exchange_lattice.cpp`)
- [ ] Implement periodic exchange (and Coulomb) assembly reusing `rys_eri_block` + lattice images (`cpu/scf/internal/…`)
- [ ] Write failing unit/integration tests for one-electron PBC S/Hcore consistency (positive definite S; hermiticity)
- [ ] Wire `RHF` (or `PeriodicRHF`) `periodic=true` path: build PBC Fock → existing SCF loop → energy on Frame
- [ ] Offline PySCF `pbc.scf.RHF` pin (tiny cell: e.g. H₂ in large cubic box **or** 2-atom prototype matching REFERENCE); write `REFERENCE.md`
- [ ] Add regression example `regressions/scf/periodic-rhf-gamma/` (public API; hard-coded golden; no live PySCF)
- [ ] Wire regression CMake under `ATV_BUILD_HF`; update queue + roadmap; Doxygen
- [ ] Run full check + test suite

## Testing strategy

**Unit:**

1. Exchange lattice sum: adding a redundant outer image shell beyond cutoff does not change K beyond tol (cutoff stability).
2. Molecular limit: very large cell + `periodic=true` recovers molecular RHF energy within a documented loose tol **or** explicit `periodic=false` bit-identical to current RHF (prefer the latter gate on the non-PBC path).

**Regression:**

- Tiny orthorhombic system, STO-3G (or twin basis), Γ RHF.
- \(|E_\mathrm{atv}-E_\mathrm{pyscf}|\le 1\times 10^{-5}\) Ha (foothold; tighten if twin is molecular-in-box) — record exact tol in REFERENCE.
- Public headers only.

**Perf:** `perf_exempt` acceptable for first foothold (no standing HF-PBC suite).

## Out of scope

- k-point sampling (Monkhorst–Pack), complex Bloch orbitals
- Periodic UHF / GHF / ROKS
- Analytic nuclear gradients / stress for periodic HF
- Hybrid DFT under PBC (needs this HF-X + XC grid PBC — separate)
- GPW/PW exact exchange (different basis)
- CUDA
- Public Config bags / factories
- Live PySCF in ctest

## Notes for implementers

- Correctness risk is **electrostatic convention mismatch** (Ewald vs truncated) between J, K, and \(E_\mathrm{nn}\). Pin **all three** consistently with PySCF settings in REFERENCE.
- Do not block on CAM/RSH; molecular attenuated ERI is a different queue (`rsh-eri-cam-b3lyp`).
- Deferred: do **not** auto-run `/mol:impl` until scheduled.

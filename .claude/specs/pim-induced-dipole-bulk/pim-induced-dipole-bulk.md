# Spec — `pim-induced-dipole-bulk`

**Slug:** `pim-induced-dipole-bulk`  
**Branch:** `wf/pim-induced-dipole-bulk`  
**Map:** MetalWalls capability **E4** (polarization — bulk PIM)  
**Queue:** `docs/developer/metalwalls-method-queue.md` wave-2  
**Created:** 2026-07-24  

## Intent

Ship classical **Polarizable Ion Model (PIM)** self-consistent induced
dipoles for a **bulk** toy molten-salt / ionic-liquid geometry as a
`cpu::Potential` Component on `cpu::md::MD` / `cpu::Evaluate` — not a new
Driver category, not electrode-joint polarization (that is E5 /
`dipole-electrode-joint`, held until this lands).

Physics (Applequist / Thole / MetalWalls bulk PIM form):

\[
\boldsymbol{\mu}_i=\alpha_i\,\mathbf{E}_i^{\mathrm{tot}},\quad
\mathbf{E}_i^{\mathrm{tot}}
  =\mathbf{E}_i^{\mathrm{perm}}
  +\sum_{j\neq i}\mathbf{T}_{ij}\cdot\boldsymbol{\mu}_j
\]

with permanent field from Frame charges (min-image, cutoff), mutual
induction via the dipole tensor \(\mathbf{T}\) with **Thole** exponential
damping (default) and optional **Tang–Toennies (TT)** short-range
damping. Polarization energy:

\[
E_{\mathrm{pol}}
  =-\tfrac12\sum_i\boldsymbol{\mu}_i\cdot\mathbf{E}_i^{\mathrm{perm}}
  \quad\text{[eV via `units::coulomb`]}
\]

Forces use Hellmann–Feynman at SCF (treat converged \(\mu\) as fixed
multipoles): charge–dipole + dipole–dipole contributions only
(**not** permanent charge–charge — that stays `PairCoulCut` / PME).

## Public API (no Config bags / factories)

```cpp
// include/atomiverse/cpu/potentials/pim_induced_dipole.h
cpu::PimInducedDipole pim(
    /*cutoff [Å]=*/12.0,
    /*alpha_by_type [Å³]=*/{0.0, 1.0, 1.5},  // index = atom type
    /*thole_a=*/0.39,                         // 0 → undamped T
    /*pol_tol [e·Å]=*/1e-10,
    /*pol_max_iter=*/100);
// optional public members after construction:
pim.tt_b = 0.0;   // [Å⁻¹]; >0 enables Tang–Toennies f_n(b r)
```

- Inherits `cpu::Potential` (many-body SCC — **not** `Pair`).
- `name()` → `"pim_induced_dipole"`.
- Writes \(E_{\mathrm{pol}}\) into `frame.meta.get(kEnergy)`; forces into
  `kAtom_fx/fy/fz`.
- Accessors: `last_pol_energy()`, `last_mu()`, `last_pol_iters()`,
  `last_pol_converged()`.
- O(N²) all-pairs foothold (NeighborConsumer withdrawn), same shape as
  `PairBuck` / early coul footholds. Production NL deferred.
- Frame I/O only; parameters = ctor args + public members.

## Oracle

| Field | Value |
| ----- | ----- |
| **Primary** | **metalwalls** (MW-only classical PIM; LAMMPS uses Drude / permanent dipole, not this SCF bulk PIM) |
| **Host reality** | No MetalWalls binary on agent PATH → offline pin is the **published multipole + Thole formula** MetalWalls implements (same as `gaussian-electrode-charge` / charging-work literature path) |
| **REFERENCE** | `regressions/md/pim-induced-dipole-bulk/REFERENCE.md` must state `oracle: metalwalls` and `oracle_version: …` |
| **ctest** | No live MW / LAMMPS |

## Geometry / observables (offline pin)

Toy bulk dimer of opposite charges (molten-salt foothold), large cubic box
so min-image ≡ free-space:

| Site | type | \(q\) [e] | \(\alpha\) [Å³] | \(\mathbf{r}\) [Å] |
| ---- | ---- | --------- | --------------- | ------------------ |
| A | 1 | +1 | 1.0 | (0, 0, 0) |
| B | 2 | −1 | 1.5 | (3, 0, 0) |

- Box \(L=50\) Å cubic, PBC on, cutoff 20 Å  
- `thole_a = 0.39`, `tt_b = 0`  
- Observables: \(E_{\mathrm{pol}}\) [eV], per-atom forces [eV/Å]  
- Dual gate: energy **and** forces vs offline pin; forces also checked
  vs central finite difference of \(E_{\mathrm{pol}}\) (internal)
- Tolerances: stated in REFERENCE.md (target \(|\Delta E|<10^{-9}\) eV,
  \(|\Delta F|<10^{-8}\) eV/Å or FD match)

## Tasks

1. Write this spec + `acceptance.md` (done first).
2. Implement `cpu::PimInducedDipole` under `cpu/potentials/` (header+cpp;
   auto-globbed by `src/cpu/potentials/*.cpp`).
3. Unit tests `tests/cpu/potentials/test_pim_induced_dipole.cpp` (SCC
   one-site, two-site energy, Thole on/off, TT path smoke, FD force).
4. Offline pin + `REFERENCE.md` + e2e
   `regressions/md/pim-induced-dipole-bulk/`.
5. Wire `atv_add_regression(md pim-induced-dipole-bulk)` under
   `ATV_BUILD_MD`.
6. MD-only configure/build; run unit + regression until PASS.
7. Perf gate: fair wall vs MW if binary appears; else `perf_exempt` with
   reason (no MW twin on host).
8. Commit on `wf/pim-induced-dipole-bulk` only.

## Out of scope

- Electrode-coupled dipole solve (E5 `dipole-electrode-joint`)
- PME / Ewald reciprocal for multipoles (real-space cutoff only)
- Drude oscillators / extended Lagrangian
- CUDA / Kokkos
- Multi-rank halo for induced multipoles
- Permanent charge–charge energy (use existing coulomb pairs)
- AIM / DAIM (E11)

## Build (binding)

MD-only flags only — see `docs/developer/metalwalls-method-loop.md`.

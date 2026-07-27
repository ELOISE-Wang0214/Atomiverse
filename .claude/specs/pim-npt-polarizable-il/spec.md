# pim-npt-polarizable-il — Bulk polarizable IL under NPT (PIM product twin)

Workflow: metalwalls-method-loop, wave 5. Queue row:
`pim-npt-polarizable-il` | Bulk polarizable IL under NPT (PIM product twin) |
polarizable | small IL box NPT | E, P, short NPT | **oracle=metalwalls**;
map §9 PIM-NPT ~ (bulk SCF PIM foothold only).

## Problem

`cpu::PimInducedDipole` (map E4, wave ≤4 foothold) writes polarization energy
and Hellmann–Feynman forces to Frame but contributes **no configurational
virial** to `frame.meta.get(kStress)`. Every barostat in the tree
(`NoseHooverNPT`, `ParrinelloRahmanNPT`, Berendsen `Barostat`, `NphIntegrator`)
reads `P_int = trace(kStress + kStressKinetic)/(3V)`. A PIM + NPT product run
therefore barostats on a pressure that is blind to the polarization
interaction — the MetalWalls PIM-NPT product twin (map §9 "Bulk polarizable
IL / PIM-NPT" = `~`) is not reachable honestly.

## Deliverable

1. **Component:** `PimInducedDipole::launch` additionally accumulates the
   pairwise Hellmann–Feynman virial
   \(W = \tfrac12\sum_{i\neq j} \mathbf r_{ij}\otimes\mathbf f_{ij}\)
   (min-image \(\mathbf r_{ij}\), the same per-pair force it already sums into
   `kAtom_f*`) into `frame.meta.get(kStress)` — Voigt order xx yy zz xy xz yz,
   units eV, LAMMPS sign (positive trace ⇔ compressive), matching
   `PairCoulCut` / `EAM` convention (unconditional write; cost is negligible
   next to the O(N²) SCC). No new public surface, no Config bag.
2. **Regression (product twin):** `regressions/md/pim-npt-polarizable-il/` —
   toy polarizable molten-salt IL (alternating ±1 charges, polarizable anions,
   `PairLJCut` repulsion + `PairCoulCut` + `PimInducedDipole`; **not**
   `PairBorn` — found gap: `PairBorn`/`PairBuck` write no `kStress`, reported
   upstream for a queue row) under `cpu::md::MD` + `NoseHooverNPT`, short run.
   Out-of-library public-API consumer, ≤ 60 s.
3. **Offline pin:** independent Python implementation of the published
   Thole-PIM formulas (MetalWalls JOSS 5:2373 (2020) bulk multipole form;
   Thole Chem. Phys. 59, 341 (1981)) evaluates the deterministic start
   configuration; `E_pol` [eV] and the PIM virial diagonal / polarization
   pressure are hard-coded pins. Script + notes under
   `tools/metalwalls_oracle/` (`inputs/` + `pin_notes_pim-npt-polarizable-il.md`).
   No oracle binary at ctest/regression time.

## Acceptance criteria

1. [x] `PimInducedDipole` accumulates its configurational virial into
   `frame.meta.get(kStress)` (Voigt, eV) consistent with the forces it writes;
   Newton-pair form `½ Σ_{i≠j} r_ij ⊗ f_ij` over the same damped kernels.
2. [x] Unit (`tests/cpu/potentials/test_pim_induced_dipole.cpp`): polarizable
   ion-pair dimer — `kStress[0]` equals the closed pair form
   `(x_0−x_1)·F_{0,x}` of the forces the component itself reports; all other
   Voigt components ~0; `kStress` untouched when every α = 0.
3. [x] Unit: undamped (`thole_a = 0`, `tt_b = 0`) non-symmetric periodic
   configuration — virial trace matches the central finite difference
   `−dE/ds` under uniform scaling of coordinates **and** box (s = 1 ± h),
   relative tolerance ≤ 1e-5 (forces are exact gradients in the undamped
   case, so this pins the virial to the energy surface).
4. [x] Regression `reg_pim-npt-polarizable-il`:
   a. start-config `E_pol` and PIM virial-trace polarization pressure match
      the hard-coded offline-Python pins (independent implementation of the
      published formulas) to ≤ 1e-8 eV / ≤ 1e-10 eV/Å³;
   b. short NPT run (NoseHooverNPT, P_target = 0): SCC converged every step,
      finite conserved quantity, bounded volume, and the start state relaxes
      in the pressure-relief direction — final design: **stretched** start
      (d = 2.9 Å, P₀ ≈ −0.12 eV/ų tension) contracts onto the stiff r⁻¹²
      wall (V_min < 0.97·V₀; an offline scan showed the compressed pin
      density would expand ballistically through the soft side), with the
      PIM virial live in the barostat pressure (leg-2 additivity);
   c. `REFERENCE.md` states oracle: metalwalls, oracle_version, citations,
      units; runtime ≤ 60 s.
5. [x] Registered: `atv_add_regression(md pim-npt-polarizable-il)` →
   targets `reg_pim-npt-polarizable-il` + `run_reg_pim-npt-polarizable-il`.
6. [x] Dual gate green under the MD-only configure (HF/PW/GPW/CUDA/Torch OFF):
   scoped build, `ctest -R '^(test_cpu_pim_induced_dipole)$'`, and
   `run_reg_pim-npt-polarizable-il`.
7. [x] perf: `perf_exempt=true` — no MetalWalls binary on PATH; LAMMPS has no
   SCF bulk-PIM twin (its polarizable path is Drude/QEq, different
   algorithm), so no fair wall-clock twin exists. Same exemption as
   `pim-induced-dipole-bulk`.

## Out of scope

Extended-Lagrangian dipole dynamics (separate queue row
`extended-lagrangian-dipoles`); neighbor-list PIM production path (deferred
in-header); Ewald/PME dipole long range; per-pair damped-gradient (dλ/dr)
force terms (foothold keeps the documented undamped-gradient approximation —
the virial must be consistent with the *implemented* forces, which is what a
barostat mechanically requires).

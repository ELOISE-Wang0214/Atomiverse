---
slug: charging-work-energy
status: open
created: 2026-07-23
track: metalwalls-method-loop
chapter: E1
grilled: false
oracle: literature
---

# charging-work-energy

## Summary

Book the **charging-work** term in the classical MD energy channel:

\[
E_{\mathrm{charge}} = -\sum_{i\in\mathrm{electrode}} \Psi_i\, Q_i
\quad[\mathrm{eV}]
\]

when electrode sites carry applied potentials \(\Psi_i\) [V] and charges \(Q_i\) [e].
This is the MetalWalls total-energy term that “includes the work necessary to
charge the electrodes” (JOSS eq. 6) and Scalfi’s \(-\,q^{\mathrm{T}}\Psi\)
contribution in the constant-potential ensemble.

**Public surface:** `cpu::electrodes::ChargingWork` — a `cpu::Potential`
Component bound on `cpu::md::MD` / `cpu::Evaluate` via `add_potential`. No
public `*Config` bag, no factory. Parameters = constructor arguments
(electrode atom types + applied potentials). Writes
`frame.meta.get(kEnergy)` and `frame.meta.get(kEnergyChargingWork)`.

**Does not** solve for electrode charges (that is `const-pot-bo-parallel-plate`).
Reads `kAtom_charge` as written by a prior const-pot Component, by the user,
or by a fixed-charge control.

## Domain basis

MetalWalls (Marin-Laflèche et al., *JOSS* **5**, 2373 (2020), DOI
10.21105/joss.02373) defines the total potential energy under constant
potential as

\[
V(Q) = V_c(Q) - \sum_{i=1}^{N} \Psi_i Q_i
\tag{6}
\]

where \(V_c\) is the Coulomb energy of the charge distribution and the second
term is the charging work. Scalfi et al., *PCCP* **22**, 10480 (2020)
(arXiv:1911.09351) identify \(q^{\mathrm{T}}\Psi = Q_L\Delta\psi\) as the work
exchanged with the external circuit when charging under fixed \(\Psi\); the
Boltzmann weight uses \(-\,q^{\mathrm{T}}\Psi\).

**Units (ATV metal / LAMMPS metal):** \(Q\) in elementary charge \(e\),
\(\Psi\) in volts; \(1\,e\cdot 1\,\mathrm{V} = 1\,\mathrm{eV}\).

## Design

| Symbol | Role |
| ------ | ---- |
| `cpu::electrodes::ChargingWork` | Potential; \(E = -\sum \Psi(t_i) Q_i\) over owned atoms |
| ctor `(electrode_types, potential_V)` | type→Ψ map; length match required |
| `kEnergyChargingWork` | serializable meta [eV]; cleared by `ZeroEnergy` |
| `kEnergy` | accumulates charging work with other potentials |
| private type→Ψ lookup | built at `build()` |

- Kind: `Potential` (energy source for MD/Evaluate validation).
- No forces from this term (BO charges / frozen electrodes: nuclear forces
  come from \(V_c\), not from the linear \(-\Psi\cdot Q\) bookkeeping).
- Owned atoms only (`nlocal`); MPI energy reduce is the existing
  `Reduce` path.
- Run **after** charge solve in the Plan when const-pot is present
  (`force` order 0 is fine if const-pot binds earlier / at `pre_force`).

### Architecture constraints (binding)

- Component on `cpu::md::MD` — **not** a new Driver category.
- No charge solve inside `KSpace*`.
- No public `*Config` / `make_*` / `build_*`.
- Map slice: **E1** (charging work row of
  `.claude/notes/metalwalls-capability-map.md` §1).

## Geometry (pin)

Symmetric parallel-plate toy (fixed charges — proxy for post-const-pot \(Q\)):

| Region | type | \(N\) | \(Q_i\) [e] | \(\Psi\) [V] |
| ------ | ---: | ----: | ----------: | -----------: |
| Left wall | 1 | 4 | \(+0.10\) | \(+0.50\) |
| Right wall | 2 | 4 | \(-0.10\) | \(-0.50\) |
| Electrolyte (optional inert) | 0 | 0–2 | 0 | — |

Positions: left wall \(z=0\), right wall \(z=10\) Å, cubic or orthorhombic box
large enough that min-image is irrelevant (bookkeeping is charge-only).

**Analytic pin (MetalWalls eq. 6 second term):**

\[
E_{\mathrm{charge}}
  = -\bigl(4\cdot 0.5\cdot 0.1 + 4\cdot(-0.5)\cdot(-0.1)\bigr)
  = -0.4~\mathrm{eV}.
\]

Equivalent Scalfi form: \(Q_L = 0.4\), \(\Delta\psi = 1.0\,\mathrm{V}\),
\(q^{\mathrm{T}}\Psi = Q_L\Delta\psi = 0.4\), \(E = -q^{\mathrm{T}}\Psi = -0.4\) eV.

## Offline pin / oracle

- **Oracle priority:** formula is MetalWalls-only energy bookkeeping (no
  LAMMPS twin for \(-\sum\Psi_i Q_i\) as a named energy term). Primary
  citation: MetalWalls JOSS eq. (6); cross-check Scalfi PCCP 2020 eq. (12).
- **No live MetalWalls/LAMMPS at ctest.** Pin is the closed-form evaluation
  above, recorded in `regressions/md/charging-work-energy/REFERENCE.md`.
- **oracle:** `literature` (MW JOSS + Scalfi); formula identity with MW.
- **oracle_version:** MetalWalls JOSS 5(53) 2373 (2020) eq.6; Scalfi et al.
  PCCP 22 10480 (2020) / arXiv:1911.09351 eq.12
- **Tolerance:** absolute \(1\times 10^{-12}\) eV (exact arithmetic).

## Tasks

- [x] Spec + acceptance (this file + `.acceptance.md`)
- [x] Offline pin → `regressions/md/charging-work-energy/REFERENCE.md`
- [x] Meta key `kEnergyChargingWork`; Frame insert + `ZeroEnergy` clear
- [x] `cpu::electrodes::ChargingWork` Potential (header + cpp)
- [x] CMake: source under `ATV_BUILD_MD`; unit test; regression under MD
- [x] Unit tests: parallel-plate pin; type mismatch throw; ZeroEnergy clears
- [x] E2e regression public API only vs pin
- [x] MD-only build green; perf note (exempt: pure O(N) bookkeeping, no oracle binary)

## Testing strategy

- **Unit:** fixed \(Q,\Psi\) frame → \(E = -0.4\) eV; empty electrode list → 0;
  ctor length mismatch throws; `kEnergyChargingWork` matches `last_work()`.
- **Regression:** `cpu::Evaluate` + `ChargingWork` on parallel-plate geometry;
  assert `kEnergy` and `kEnergyChargingWork` within \(10^{-12}\) eV of pin.

## Out of scope

- BO / CG electrode charge solve (`const-pot-bo-parallel-plate`)
- Gaussian electrode kernels, Thomas–Fermi, piston
- Capacitance observers, potential-profile analysis
- CUDA path
- Forces from \(\partial Q/\partial R\) (const-pot / \(V_c\) territory)

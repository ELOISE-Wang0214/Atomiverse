# Reference — bulk polarizable IL under NPT (PIM product twin, map §9)

oracle: metalwalls  
oracle_version: MetalWalls JOSS **5**, 2373 (2020) bulk PIM multipole form;
Thole, *Chem. Phys.* **59**, 341 (1981) exponential damping; offline numpy
pin (no MetalWalls binary on agent host) —
`tools/metalwalls_oracle/inputs/pim_npt_polarizable_il_pin.py`, numpy 2.5.1.
NPT integrator: MTK (Tuckerman, *J. Phys. A* **39**, 5629 (2006); Shinoda,
*Phys. Rev. B* **69**, 134103 (2004)), LAMMPS `fix npt iso` convention.  
perf_exempt: true  
perf_exempt_reason: No MetalWalls binary on PATH; LAMMPS has no SCF bulk-PIM
twin (its polarizable paths are Drude/QEq — a different algorithm), so no
fair wall-clock twin exists. Same exemption as `pim-induced-dipole-bulk`.

## What this regression pins

The **PIM-NPT product twin**: the polarizable-ion-model virial must reach
the shared stress meta (`kStress`) every barostat reads, so a bulk
polarizable ionic melt barostats on the **full** pressure including the
polarization contribution
\(W=\tfrac12\sum_{i\neq j}\mathbf r_{ij}\otimes\mathbf f_{ij}\)
(Hellmann–Feynman pair forces of the published Thole-PIM form, Voigt, eV,
LAMMPS sign).

## Geometry — deterministic toy polarizable IL

4×4×4 simple-cubic rock-salt lattice, d = 2.35 Å (compressed start), cubic
box L = 9.4 Å, N = 64: +1 cations (type 0, α = 0) / −1 anions (type 1,
α = 3.0 ų) alternating as (−1)^(ix+iy+iz); deterministic displacement
`disp[a] = 0.15·sin(1.7·(i+1) + 2.3·a)` Å breaks centro-symmetry (a perfect
rock-salt lattice has E_pol ≡ 0). PIM: cutoff 4.0 Å, thole_a = 0.39,
tt_b = 0. k = 14.399645 eV·Å/e².

## Pin values (hard-coded; offline independent numpy implementation)

| Observable | Value | Units | Tol |
| --- | --- | --- | --- |
| E_pol | −8.695749870555906 | eV | 1e-8 |
| W_xx | −10.66184373345355 | eV | 1e-7 |
| W_yy | −12.125914211727187 | eV | 1e-7 |
| W_zz | −9.657852793291474 | eV | 1e-7 |
| W_xy | −0.2325271422922529 | eV | 1e-7 |
| W_xz | −2.2469237837238 | eV | 1e-7 |
| W_yz | +1.375774440282161 | eV | 1e-7 |
| P_pol = tr(W)/(3V) | −0.013021203851134545 | eV/ų | 1e-10 |

These are literature-formula evaluations (published PIM/Thole equations,
independent implementation), not a golden dump of Atomiverse output.

## Legs

1. **PIM pin:** `Evaluate` + `PimInducedDipole` alone on the start config —
   E_pol, all six Voigt virial components, and P_pol vs the table above.
2. **Stress additivity:** trace(kStress) of the LJ+coul/cut+PIM stack minus
   the LJ+coul/cut stack equals the pinned PIM trace — the polarization
   pressure reaches the meta the barostat reads.
3. **Short NPT product run:** same lattice **stretched** to d = 2.9 Å
   (start tension P₀ ≈ −0.12 eV/ų — a stretched ionic melt contracts onto
   the stiff r⁻¹² wall, which is self-limiting; the compressed pin-leg
   density would expand ballistically through the soft side).
   `cpu::md::MD` + `PairLJCut{0.01 eV, 3.0 Å}` + `PairCoulCut` (both
   virial-correct) + `PimInducedDipole` +
   `NoseHooverNPT{T=800 K, τ_T=100 fs, P=0, τ_P=500 fs, dt=1 fs}` +
   `Virial` (kinetic stress) + `InitVelocity` + an out-of-library
   `TrackVolume` Criterion, 250 steps. Invariants: per-step SCC convergence
   (the component throws otherwise), finite MTK conserved quantity,
   contraction happened (V_min < 0.97·V₀), bounded
   (0.5·V₀ < V_min, V_max < 1.2·V₀).

Runtime ≪ 60 s (N = 64, O(N²) SCC per step).

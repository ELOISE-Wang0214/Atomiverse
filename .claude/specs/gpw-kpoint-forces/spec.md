# gpw-kpoint-forces

## Goal

Ship **analytic ionic forces for closed-shell GPW RKS at Monkhorst–Pack
k-points** — the follow-on named in `gpw_solver.h` after Γ analytic forces
(`gpw-forces-water` ✓) and k-point energy (`gpw-kpoints-si` ✓).

Today `gpw_forces_analytic` / `assemble_gpw_forces` are Γ-only (real P/W,
real lattice sums). The k-path SCF stores only `eps_k` / `k_frac` and
discards P_k; forces with `k_mesh ≠ {1,1,1}` either ignore the mesh
(`copy_gpw` in `gpw_forces.cpp` drops k fields) or are wrong.

## Scope

**In**

- Analytic forces F_A = −∂E/∂R_A for periodic closed-shell RKS GPW with
  `k_mesh` / explicit `k_points` (same branch constraints as the energy
  foothold: single grid, `local_pp_on_grid`, no UKS/OT/D3).
- k-summed Pulay (W·dS), core (P·dT, P·dV_nl), grid (V_eff against Bloch
  χ derivatives), ionic Hellmann–Feynman, local-PP HF (grid route).
- Store occupied MO coefficients per k on `GpwResult` so forces rebuild
  P_k / W_k after SCF without a second solve.
- Unit: analytic vs FD energy gradient on a cheap k-mesh H₂ toy.
- Regression: offline CP2K ENERGY_FORCE twin (or frozen first-green +
  identities if no pin host) under `regressions/gpw/gpw-kpoint-forces/`.

**Out**

- k-point UKS forces, hybrid/meta-GGA, stress-at-k (separate slug),
  CPMD/XL at k, GPWEvaluator `k_mesh` members (callers use `GpwDft` +
  `gpw_forces_analytic` for this foothold).

## Public surface

- `gpw_forces_analytic(atoms, cfg)` — already public; gains k-path.
- `GpwResult` payload: add `C_occ_k` (complex occupied MOs per k).
- No public `*Config`, no `make_*` / `build_*`.

## Physics

At the SCF fixed point (closed shell, occupancy 2):

\[
F_A = \sum_k w_k\Bigl[
  \mathrm{Tr}(W_k\,\partial S_k/\partial R_A)
  - \mathrm{Tr}(P_k\,\partial H^\mathrm{core}_k/\partial R_A)
\Bigr]
+ \text{grid HF (fixed }V_\mathrm{eff}\text{)}
+ \text{ionic / local-PP HF}.
\]

\(W_k = \sum_i 2\varepsilon_i^k\,c_i c_i^\dagger\). Complex contractions use
\(2\,\mathrm{Re}(\cdots)\) on the bra-side center derivatives (Γ limit
recovers the existing real assembly).

Refs: VandeVondele et al., CPC **167**, 103 (2005); ATV Γ path
`internal/gpw_forces.cpp`; k machinery `internal/kpoints.h`.

## Acceptance (summary)

See `acceptance.md`. Dual gate: unit analytic↔FD + regression pin/identities.
Perf: \(T_\mathrm{atv}/T_\mathrm{cp2k} \le 10\) or `perf_exempt` if no fair
CP2K wall.

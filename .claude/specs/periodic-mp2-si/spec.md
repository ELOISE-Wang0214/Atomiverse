---
slug: periodic-mp2-si
status: approved
created: 2026-07-25
revised: 2026-07-25
domain: gpw
queue: cp2k-method-loop (wave 3)
depends_on: []          # solid Γ GPW energy ✓ (already on dev)
unblocks: [periodic-rpa-si]
branch: wf/periodic-mp2-si
---

# Spec: `periodic-mp2-si` — Γ-point canonical MP2 in GPW (tiny Si cell)

**Slug:** `periodic-mp2-si`
**Domain:** GPW post-SCF correlation (CPU)
**Queue row:** `periodic-mp2-si` | Periodic RI-MP2 (GPW) | gpw | Si (tiny cell)
| MOLOPT+GTH | e_tot Ha | roadmap §10 ○ Periodic MP2; CP2K MP2
**Unblocks:** `periodic-rpa-si` (needs the same Γ-point MO-ERI-via-Poisson
machinery + denominators)

## Summary

Deliver the **periodic MP2 foothold**: second-order Møller–Plesset
correlation energy on top of a converged Γ-point closed-shell GPW SCF, with
the occupied–virtual MO repulsion integrals \((ia|jb)\) computed **the GPW
way** — pair densities collocated on the uniform density grid, periodic FFT
Poisson solve per pair, grid quadrature — exactly CP2K's Quickstep
wavefunction-correlation integral route (Del Ben, Hutter, VandeVondele,
*JCTC* **8**, 4177 (2012); CP2K `&WF_CORRELATION` `MP2 METHOD
DIRECT_CANONICAL`). v1 is **canonical (non-RI)**: for the tiny-cell foothold
the full \((ia|jb)\) tensor is small (`n_occ·n_virt` Poisson solves); the RI
factorization \((ia|P)\,J^{-1}_{PQ}\,(Q|jb)\) is an optimization follow-on
that needs an auxiliary-basis table (cFIT) not yet transcribed in
`basis_data/`. The MO-ERI machinery and denominator handling shipped here
are precisely what `periodic-rpa-si` consumes.

## Domain basis

Γ-point periodic MOs from the converged GPW RKS solve
(`GpwResult::C`, `GpwResult::eps`, S-orthonormal after canonical
orthogonalization):

\[ \psi_p(\mathbf r) = \sum_\mu C_{\mu p}\,\chi_\mu(\mathbf r),\qquad
   \chi_\mu(\mathbf r) = \sum_{\mathbf T} \varphi_\mu(\mathbf r - \mathbf T) \]

(the Γ-point lattice-summed AO values already produced by
`build_compact_basis_grid`). Occupied–virtual pair densities
\(\rho_{ia}(\mathbf r) = \psi_i(\mathbf r)\psi_a(\mathbf r)\) are neutral
(\(\int\rho_{ia} = S_{ia} = 0\)), so the periodic Coulomb with the G=0
component dropped (`solve_poisson`, neutral-background convention) is
well-defined — no divergence at Γ:

\[ (ia|jb) = \int\!\!\int
   \frac{\rho_{ia}(\mathbf r)\,\rho_{jb}(\mathbf r')}{|\mathbf r-\mathbf r'|}
   \,d\mathbf r\, d\mathbf r'
   = \sum_g w\, V_{ia}(\mathbf r_g)\, \rho_{jb}(\mathbf r_g) \]

Closed-shell canonical MP2 (real orbitals):

\[ E^{(2)} = -\sum_{ijab}
   \frac{(ia|jb)\,[\,2(ia|jb) - (ib|ja)\,]}
        {\varepsilon_a + \varepsilon_b - \varepsilon_i - \varepsilon_j},
   \qquad
   E^{(2)} = E_{\rm OS} + E_{\rm SS} \]

with \(E_{\rm OS} = -\sum (ia|jb)^2/\Delta\),
\(E_{\rm SS} = -\sum (ia|jb)[(ia|jb)-(ib|ja)]/\Delta\) (the same
decomposition as molecular `cpu::scf::MP2Result`). MP2 on KS orbitals is the
standard periodic-correlation foothold (identical machinery to the RPA/double
hybrid ingredient; a true MP2@HF reference awaits periodic exact exchange,
which is a separately blocked row `periodic-hybrid-gpw`).

## Public API (target shape — no Config bag, no factory)

Mirrors the molecular convention `RHF::mp2(frame)` → method member on the
concrete GPW method:

```cpp
cpu::gpw::GpwDft m;                 // configured + build() + run() as today
m.run(atoms);                       // converged Γ RKS
m.mp2_n_virt = 0;                   // member: virtual-space cap (0 = all)
const cpu::gpw::Mp2Result& r = m.mp2(atoms);   // Γ-point GPW MP2
// r.e_corr [Ha] < 0, r.e_corr_os, r.e_corr_ss, r.e_total = E_SCF + e_corr
```

- New public header `include/atomiverse/cpu/gpw/mp2.h`: `Mp2Result` value
  type only (like `cpu/scf/mp2.h`).
- New `GpwDft` members: `mp2_n_virt` (0 = all kept virtuals),
  `mp2_os_scale` / `mp2_ss_scale` (Grimme SCS scales, default 1.0 — same
  member names as the molecular methods).
- `GpwDft::mp2(atoms)` throws `std::logic_error` before a converged `run()`,
  `std::invalid_argument` for UKS (foothold is closed-shell RKS only) and
  for a vanishing HOMO–LUMO denominator (< 1e-6 Ha).
- Kernels (MO-on-grid, MO-ERI-via-Poisson, closed-shell doubles sum) are
  **private workspace**: `cpu/gpw/internal/mp2_kernels.h`, unit-testable,
  never on the public path.

## Non-goals (follow-ons)

- RI factorization with cFIT auxiliary basis (needs a `basis_data/` table).
- k-point sampling beyond Γ; UMP2; gradients; frozen-core (GTH valence-only
  makes frozen-core moot).
- `GPWEvaluator` / Frame-meta integration (single-point MP2 driver wiring).
- `periodic-rpa-si` (consumes this machinery; separate queue row).

## Implementation plan

1. `internal/mp2_kernels.h` + `src/cpu/gpw/internal/mp2_kernels.cpp`:
   - `mo_on_grid(cbg, C, p0, np)` — MO values from the compact
     lattice-summed AO grid (CSR expansion).
   - `mo_eri_ovov(grid, psi_occ, n_occ, psi_virt, n_virt, nonperiodic)` —
     `(ia|jb)` tensor via per-pair Poisson (periodic FFT or MT), using
     `(ia|jb) = (jb|ia)` symmetry.
   - `mp2_rhf_energy(v_ovov, eps_occ, eps_virt, &os, &ss)` — canonical
     closed-shell doubles sum.
2. `GpwDft::mp2` in `src/cpu/gpw/mp2.cpp`: rebuild grid + compact basis at
   the given geometry (same `make_uniform_grid` / `build_compact_basis_grid`
   parameters as `run()`), validate, compose kernels, fill `Mp2Result`.
3. Unit tests `tests/cpu/gpw/test_cpu_gpw_mp2.cpp` (single-function, fast):
   spin-orbital cross-check of the closed-shell sum; MO-ERI vs Poisson
   Hartree-energy identity on a toy 2-function cell; guard throws.
4. Regression `regressions/gpw/periodic-mp2-si/` (public API only, ≤ 60 s):
   Si₂ diamond fragment, cubic a = 5.43 Å, DZVP-MOLOPT-GTH + GTH-PBE-q4,
   periodic Γ (settings tuned to budget; see REFERENCE.md), hard internal
   gates + first-green frozen band + CP2K input twin for offline re-pin.

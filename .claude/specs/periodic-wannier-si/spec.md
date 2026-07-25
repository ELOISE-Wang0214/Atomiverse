# periodic-wannier-si — Periodic multi-k MLWF (Marzari–Vanderbilt) on PW

**Loop:** cp2k-method-loop, assign wave 1 2026-07-25. **Domain:** response / pw.
**Backend:** CPU plane-wave. **Branch:** `wf/periodic-wannier-si`.
**Depends on:** `wannier-centers` ✓ (molecular Foster–Boys foothold),
`bands-dos-si` ✓ / `berry-polarization` ✓ (NSCF-on-frozen-Veff, `rebuild_k_tables`,
occupied multi-k overlaps).

## What

**Maximally localized Wannier functions (MLWF)** for a periodic insulator via the
Marzari–Vanderbilt (MV) multi-k construction (Marzari & Vanderbilt,
*Phys. Rev. B* **56**, 12847 (1997); Marzari et al., *Rev. Mod. Phys.* **84**,
1419 (2012)):

1. Restricted fixed-occupation PW SCF (`run_restricted_scf`).
2. NSCF Davidson on the frozen \(V_\mathrm{eff}\) at every Monkhorst–Pack k-point
   for the occupied manifold only.
3. Neighbor shell \(\{\mathbf b\}\) on the mesh (± one step along each reciprocal
   axis); overlap matrices
   \(M_{mn}^{(\mathbf k,\mathbf b)}=\langle u_{m\mathbf k}|u_{n,\mathbf k+\mathbf b}\rangle\)
   with periodic-gauge G-shift when \(\mathbf k+\mathbf b\) wraps the BZ.
4. Unitary gauge \(U(\mathbf k)\): optional sp³ / bond-centered Gaussian
   **projection** initialization (standard for diamond), then MV steepest-descent
   polish minimizing the SMV discretized spread \(\Omega\).
5. Wannier centers from the Berry-phase formula
   \(\mathbf r_n= -\frac1{N_k}\sum_{\mathbf k,\mathbf b} w_b\,\mathbf b\,
   \operatorname{Im}\ln M_{nn}^{(\mathbf k,\mathbf b)}\) [Å].

**Textbook pin (parameter-free geometry):** diamond Si has **four** valence
MLWFs, one per tetrahedral Si–Si bond, with centroids at the four bond
midpoints of the 2-atom primitive cell (Marzari & Vanderbilt 1997 §IV;
Wannier90 Si example). Distances after min-image lattice wrap.

CP2K twin: Quickstep `&LOCALIZE` / Wannier centers on a comparable Si cell
(GTH+Gaussian basis — absolute spreads not 1:1 with ONCV PW; bond-midpoint
geometry is basis-independent). Offline input echo under
`tools/cp2k_oracle/inputs/`.

## Public surface (architecture-conform)

`include/atomiverse/cpu/planewave/periodic_wannier.h` — concrete class, members
as parameters, Frame in / result out; **no** `*Config` bag, **no** factories:

```cpp
cpu::pw::PW method;                  // upf, ecut_ry, k_mesh, xc, …
cpu::pw::PeriodicWannier mlwf;       // members: k_mesh, max_iter, tol, …
auto r = mlwf.compute(method, frame);
// r.centers[i] [Å], r.spreads[i] [Å²], r.total_spread, r.converged
```

Pure free helpers (unit-testable without SCF):
`diamond_bond_midpoints(a_ang)`, `mv_pair_weight_cubic(b_len)`,
`center_from_diagonal_phases(...)` geometry checks.

## Internal work

- Reuse `internal::run_restricted_scf`, `rebuild_k_tables`, Davidson NSCF
  (same pattern as `berry_polarization.cpp` / `bands_dos.cpp`).
- Full complex overlap **matrices** (Berry only needed \(\det S\); MLWF needs
  \(M_{mn}\)).
- MP neighbor table + SMV weights \(w_b=3/(N_b|\mathbf b|^2)\) for the six
  cubic first neighbors.
- Projection init onto four bond-centered Gaussians (Fourier coefficients on
  each k-sphere) + thin MV SD polish (`max_iter`, `tol`, `step` members).
- Centers / spreads in Å / Å² via `units::kBohrInAngstrom`.

## Out of scope

- Molecular Foster–Boys (already `cpu::scf::WannierCenters`).
- GPW-grid MLWF, disentanglement for entangled bands, spin-polarized /
  metallic cases (throw).
- Wannier interpolation of bands, transport, CUDA.
- Fair CP2K wall when no binary on host → `perf_exempt`.

## Case

Bulk Si diamond, \(a=5.43\) Å, 2-atom fcc cell, in-tree ONCVPSP LDA UPF,
\(E_\mathrm{cut}=15\) Ry, SCF/MLWF k-mesh \(2\times2\times2\), LSDA.
Acceptance: four centers match tetrahedral bond midpoints within **0.15 Å**
(min-image); SCF converged; `n_wann == 4`; localization reduces or holds
total spread vs the unrotated gauge (projection/MV).

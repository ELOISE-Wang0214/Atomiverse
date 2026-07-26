# Spec: `gpw-mbd-kpoint-forces`

## Goal

Close the residual that **Γ-point** `assemble_gpw_forces` folds MBD@rsSCS
ionic forces when `GpwDft::mbd = true`, but the **multi-k** assembly
(`assemble_gpw_forces_kpoints`) counts MBD only for mutual-exclusion and
**skips** adding \(F_{\mathrm{MBD}}\). Ship product multi-k analytic forces
with free-atom / periodic MBD@rsSCS matching the Γ path.

## Residual (tree on assign)

In `src/cpu/gpw/internal/gpw_forces.cpp` (`assemble_gpw_forces_kpoints`):

- `n_disp` includes `cfg.mbd`.
- Entry condition is only `d3 || d3bj || d4 || ts || xdm` — **no `mbd`**.
- No MBD branch; d3/d4/ts/xdm fold correctly on the k-path.

Γ path (`assemble_gpw_forces`) already has a full MBD branch (`gpw-mbd-forces` ✓).

## Scope

| In | Out |
| -- | --- |
| Wire MBD forces into multi-k force assembly (isolated + `mbd_k_mesh` periodic) | Density-chain Hirshfeld \(\partial(V/V_{\mathrm{free}})/\partial R\) |
| Same mutual-exclusion rules as Γ | New public API / `*Config` / factories |
| Unit test + `regressions/gpw/gpw-mbd-kpoint-forces` | Live CP2K in ctest |
| Offline libMBD / dual-gate pins | Stock CP2K `TYPE MBD` (2026.1 has none) |

## API (no public Config bags / factories)

- Existing: `GpwDft::mbd`, `GpwDft::mbd_k_mesh`, `GpwDft::k_mesh` / `k_points`.
- Product: `gpw_forces_analytic` / `gpw_forces_from_scf` / evaluator with
  `compute_forces=true` on a multi-k SCF must add geometry-only MBD forces
  exactly as at Γ.
- Internal: reuse `mbd_rsscs_pbe` (forces already product-quality).

## Science

- MBD@rsSCS: Ambrosetti et al., *J. Chem. Phys.* **140**, 18A508 (2014).
- libMBD: Hermann et al., *J. Chem. Phys.* **159**, 174802 (2023).
- PBE: \(\beta=0.83\), Fermi \(a=6\); free-atom tables `basis_data/ts_tables.h`.
- Electronic k-mesh and MBD k-mesh are independent: multi-k **KS** path must
  still fold geometry-only \(F_{\mathrm{MBD}}\) (isolated or Ewald MBD).
- Public Frame energy remains eV; force pins reported Ha/Bohr (internal).

## Oracle pin (offline)

| Quantity | Source | Notes |
| -------- | ------ | ----- |
| E_disp (isolated free-atom) | libMBD Fortran pin (same as `gpw-mbd-ice` / `gpw-mbd-forces`) | −1.4928804171407961×10⁻³ Ha water dimer |
| F_MBD dual-gate | Independent central FD of geometry-only \(E_{\mathrm{MBD}}\) | max\|ΔF\| ≤ 1e-6 Ha/Bohr |
| Multi-k wiring | \(F^{\mathrm{k+MBD}}-F^{\mathrm{k,plain}}\approx F_{\mathrm{MBD}}\) | One multi-k SCF; two force assemblies |
| Full multi-SCF FD of \(E_{\mathrm{tot}}\) | Soft / unit-tier if present | Not required for acceptance |

CP2K 2026.1 stock `PAIR_POTENTIAL TYPE` has no MBD — same policy as sibling
MBD slugs. Twin notes under `tools/cp2k_oracle/inputs/gpw-mbd-kpoint-forces/`.

## Acceptance

1. Multi-k `assemble_gpw_forces_kpoints` with `mbd=true` adds MBD forces
   (isolated and periodic `mbd_k_mesh` switch identical to Γ).
2. Geometry-only dual-gate still holds (product vs independent FD).
3. Product wiring: on multi-k SCF, max\|F_mbd − F_plain − F_geom\| ≤ 1e-8 Ha/Bohr
   (geometry-only; density unchanged by MBD).
4. Regression `reg_gpw-mbd-kpoint-forces` green, ≤ 60 s.
5. No public `*Config` / `make_*` / `build_*`.
6. Perf: `within_10x` **or** `perf_exempt` (geometry-only / unfair vs full QS).

## Files (expected)

- `src/cpu/gpw/internal/gpw_forces.cpp` (k-path MBD branch)
- `tests/cpu/gpw/test_cpu_gpw_mbd.cpp` (multi-k force wiring unit)
- `regressions/gpw/gpw-mbd-kpoint-forces/`
- `tools/cp2k_oracle/inputs/gpw-mbd-kpoint-forces/`
- `regressions/CMakeLists.txt` registration

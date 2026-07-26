# Spec: `gpw-mbd-forces`

## Goal

Ship **ionic forces** for free-atom / volume-scaled **MBD@rsSCS** on GPW so
`gpw_forces_*` / `GPWEvaluator(compute_forces=true)` include the MBD
contribution (today energy + Γ/periodic stress only; forces energy-only v1).

## Scope

| In | Out |
| -- | --- |
| Geometry-only MBD forces \(F_A=-\partial E_{\mathrm{MBD}}/\partial R_A\) | Density-chain Hirshfeld \(\partial(V/V_{\mathrm{free}})/\partial R\) |
| Isolated free-atom MBD (default `mbd_k_mesh={0,0,0}`) | Multi-k electronic GPW force path dispersion (k-path assemble still no disp) |
| Periodic Ewald MBD forces when `mbd_k_mesh` set | D4 forces (`gpw-d4-forces` sibling) |
| Wire into Γ `assemble_gpw_forces` | CP2K live oracle in ctest |
| Unit test + `regressions/gpw/gpw-mbd-forces` | Stock CP2K `TYPE MBD` (CP2K 2026.1 has no stock MBD pair type) |

## API (no public Config bags / factories)

- Internal: `MbdResult` gains `forces` [Ha/Bohr]; `mbd_rsscs_pbe(...)` fills them.
- Product: `GpwDft::mbd = true` + existing force entry points; no new public symbols.
- Forces are the **negative gradient of the same energy** already matched to
  libMBD (isolated pin of `gpw-mbd-ice`). Implementation may use high-order
  central finite differences of the geometry-only energy (product-quality
  gradient; matches libMBD analytic forces to \(O(h^2)\)).

## Science

- MBD@rsSCS: Ambrosetti et al., *J. Chem. Phys.* **140**, 18A508 (2014).
- libMBD: Hermann et al., *J. Chem. Phys.* **159**, 174802 (2023).
- PBE: \(\beta=0.83\), Fermi \(a=6\); free-atom tables `basis_data/ts_tables.h`.
- Public Frame energy still eV; force pins reported Ha/Bohr (internal) with
  conversion documented when Frame eV/Å is compared.

## Oracle pin (offline)

| Quantity | Source | Notes |
| -------- | ------ | ----- |
| E_disp | libMBD Fortran / pure-Python (existing ice pin) | −1.4928804171407961×10⁻³ Ha isolated water dimer |
| F_MBD | Central FD of ATV/libMBD-class energy **or** offline libMBD gradients | Dual-gate: product forces vs independent FD \(h=10^{-5}\) Bohr; max\|ΔF\| ≤ 1e-6 Ha/Bohr |
| Full GPW | Optional: product \(F = F_{\mathrm{elec}}+F_{\mathrm{MBD}}\) vs multi-SCF FD of \(E_{\mathrm{tot}}\) on one O coordinate | Soft gate if present |

CP2K 2026.1 stock `PAIR_POTENTIAL TYPE` has no MBD — same policy as
`gpw-mbd-ice` / `gpw-ts-ice`. Twin notes under
`tools/cp2k_oracle/inputs/gpw-mbd-forces/`.

## Acceptance

1. `mbd_rsscs_pbe` returns forces; isolated water-dimer E_disp still within
   **1e-12 Ha** of the libMBD pin.
2. Unit test: max\|F_product − F_fd\| < **1e-6** Ha/Bohr (independent FD step).
3. Unit test: \(\|\sum_A F_A\| < 1e-8\) Ha/Bohr (translation).
4. `assemble_gpw_forces` with `mbd=true` adds MBD forces (not energy-only).
5. Regression `reg_gpw-mbd-forces`: E_disp pin + force dual-gate; ≤ 60 s.
6. No public `*Config` / `make_*` / `build_*`.
7. Scoped build+test green; `perf_exempt` allowed (geometry-only kernel vs
   full SCF wall unfair, same policy as ice E_disp pins).

## Files (expected)

- `include/.../mbd_dispersion.h`, `src/.../mbd_dispersion.cpp`
- `src/.../gpw_forces.cpp` (MBD branch)
- `tests/cpu/gpw/test_cpu_gpw_mbd.cpp`
- `regressions/gpw/gpw-mbd-forces/`
- `tools/cp2k_oracle/inputs/gpw-mbd-forces/`
- `regressions/CMakeLists.txt` registration

# Acceptance — `gpw-xdm-forces`

## Correctness

| Check | Tol / criterion |
| ----- | --------------- |
| \|E_disp − pin\| (periodic free-atom XDM, rcut 50 Å) | **1×10⁻⁸ Ha** |
| max\|F_product − F_fd\| (h=1e-5 Bohr, isolated + periodic kernels) | **1×10⁻⁶ Ha/Bohr** |
| \|Σ_A F_A\| (translation) | **1×10⁻⁸ Ha/Bohr** |
| max\|F − F_dual_python\| on documented sample components | **1×10⁻¹⁰ Ha/Bohr** |
| `GpwDft`/`xdm=true` short SCF `e_disp` matches pin | **1×10⁻⁸ Ha** |
| Unit: `gpw_forces_analytic` + `xdm` vs multi-SCF FD on H₂ | max\|ΔF\| ≤ **5×10⁻⁵** Ha/Bohr |

## Performance

- **`perf_exempt = true`**
- **Reason:** CP2K 2026.1 has no stock `TYPE XDM`; primary hard gate is dual
  pure-Python E_disp + force FD dual-gate (same policy as `gpw-xdm-ice` /
  `gpw-mbd-forces`). Geometry-only kernel is microseconds; SCF wall is not a
  fair dispersion-force twin.

## Done when

- `regressions/gpw/gpw-xdm-forces/` + `REFERENCE.md` green via
  `run_reg_gpw-xdm-forces`
- `test_cpu_gpw_xdm` + `test_cpu_gpw_forces` (XDM product) green under scoped
  ctest
- Offline note under `tools/cp2k_oracle/inputs/gpw-xdm-forces/`
- Spec acceptance_met when dual gate green

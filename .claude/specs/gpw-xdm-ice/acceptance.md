# Acceptance — gpw-xdm-ice

## Hard gates

1. **Unit** `test_cpu_gpw_xdm`:
   - Isolated free-atom XDM water dimer \(E_{\mathrm{disp}}\) matches offline
     dual-Python pin to \(1\times 10^{-12}\) Ha.
   - Periodic (12 Å, rcut 50 Å) \(E_{\mathrm{disp}}\) matches offline pin to
     \(1\times 10^{-12}\) Ha.
   - Analytic forces match central FD to \(1\times 10^{-7}\) Ha/Bohr (isolated).
   - `GpwDft::xdm` folds `e_disp` into `e_total` (geometry-only; short SCF).

2. **Regression** `reg_gpw-xdm-ice`:
   - Periodic free-atom XDM \(E_{\mathrm{disp}} =
     -2.114815761287903\times 10^{-3}\) Ha within \(1\times 10^{-8}\) Ha
     (see `regressions/gpw/gpw-xdm-ice/REFERENCE.md`).

3. **Architecture**: no public `*Config`; no `make_*` factories; Frame I/O
   unchanged; public energy still eV on Frame when evaluator path is used;
   internal Ha for `e_disp` cited.

## Perf

**perf_exempt = true**: geometry-only \(O(N^2 N_{\mathrm{images}})\)
microseconds; full twin SCF wall is unfair for an E_disp pin (same policy as
`gpw-ts-ice` / `gpw-pbe-d3bj-ice`).

## Oracle

Offline dual pure-Python lattice sum (formula + tables in REFERENCE). CP2K
2026.1 has **no** stock `TYPE XDM` — documented like TS/libMBD.

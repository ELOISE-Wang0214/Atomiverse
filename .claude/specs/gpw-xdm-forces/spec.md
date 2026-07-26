# Spec — `gpw-xdm-forces`

**Slug:** `gpw-xdm-forces`  
**Branch:** `wf/gpw-xdm-forces`  
**Domain:** GPW (periodic DFT / dispersion)  
**Oracle:** dual pure-Python free-atom XDM lattice sum (CP2K has no stock `TYPE XDM`)  
**Dep:** `gpw-xdm-ice` ✓ (E_disp foothold)

## Intent

Ship a **dedicated offline pin + e2e regression** for **analytic ionic forces**
from free-atom XDM dispersion on the same water-dimer geometry as
`gpw-xdm-ice`. Force assembly already folds `xdm_pbe` forces into
`gpw_forces_*` (Γ + multi-k); this slug closes the product dual-gate that the
energy-only foothold left open.

## Scope

**In:**

1. Hard E_disp pin (periodic free-atom XDM, same as `gpw-xdm-ice`).
2. Analytic \(F_A = -\partial E_{\mathrm{disp}}/\partial R_A\) dual-gated vs
   independent central FD of the geometry-only kernel.
3. Offline force-component sample (dual-Python) for isolated + periodic.
4. Short GPW SCF with `xdm=true` reports the same `e_disp`.
5. Unit tier: product `gpw_forces_analytic` with `xdm=true` matches multi-SCF
   energy FD on cheap H₂ (parity with D4 unit gate).

**Out:** density-dependent Hirshfeld / BR-hole XDM; stress-only pin (separate);
CP2K live QS (no stock XDM); molecular post-HF; CUDA.

## Public API

No new public surface. Existing members:

- `cpu::gpw::GpwDft::xdm` / `GPWEvaluator::xdm`
- `gpw_forces_analytic` / `gpw_forces_from_scf` fold geometry-only XDM when
  `xdm=true`

## Acceptance (summary)

See `acceptance.md`. Dual gate: correctness (pin + FD) + `perf_exempt`
(geometry-only primary hard gate; no fair CP2K ENERGY_FORCE twin).

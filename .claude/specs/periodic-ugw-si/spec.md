# Spec: periodic-ugw-si

**Slug:** `periodic-ugw-si`  
**Method:** Unrestricted spectral G0W0 (UGW@UKS) — Γ-point GPW  
**Oracle:** offline CP2K / dual-gate vs closed-shell spectral G0W0  
**Branch:** `wf/periodic-ugw-si`  
**Domain:** `gpw`

## Goal

Ship **Γ-point spectral G0W0 on a collinear UKS reference** (`GpwDft::g0w0`
after UKS `run()`), product twin of closed-shell
`periodic-gw-si`. Screening is **spin-conserving direct RPA** (Coulomb-only
Casida over α→α ⊕ β→β transitions); correlation self-energy is the
frequency-exact spectral sum over RPA poles (molecular spectral form with
spins explicit). Exact-exchange diagonals use **MT Coulomb** per spin;
local \(v_{\mathrm{xc}}^\sigma\) is the spin-resolved SCF functional.

**Multi-k UGW is out of scope** (follow-on after multi-k UMP2). **URPA**
energy product is `periodic-urpa-si` (held behind this slug).

## Public API (no Config, no factories)

```cpp
cpu::gpw::GpwDft m;
m.n_alpha = 4; m.n_beta = 4;   // UKS
m.xc = XcKind::PBE;
m.mp2_ri = true;
m.build(atoms);
m.run(atoms);
const cpu::gpw::G0W0Result& gw = m.g0w0(atoms);
// gw.gap_ev [eV]; gw.qp_energy_ha = α-channel spectrum (primary report)
```

- Parameters remain **members** of `GpwDft` (virtual cap `mp2_n_virt`,
  `mp2_ri`). No `G0W0Config` / `make_ugw`.
- Result type stays `G0W0Result` (extend docs; primary spectrum = α when
  UKS; dual-gate vs RKS when \(n_\alpha=n_\beta\) and \(C^\alpha\approx C^\beta\)).

## Science

| Piece | Choice |
|-------|--------|
| Reference | Γ UKS GPW PBE, artificial sc-Si₂ (same deck as `periodic-gw-si` / `periodic-ump2-si`) |
| Screening | spin-conserving dRPA: \(A=\Delta\varepsilon\delta+V\), \(B=V\) with \(V_{ia\sigma,jb\tau}=(ia\sigma\|jb\tau)\); **no** closed-shell factor 2 (spins explicit) |
| Σ_c | spectral sum over RPA poles; \(W_s=\sum_{ia\tau}(X+Y)_{ia\tau}^s V_{ia\tau}\); same-spin \(G_\sigma W\) |
| Σ_x | MT exchange per spin \(\Sigma_{x,p\sigma}=-\sum_i^{\mathrm{occ}_\sigma}(pi\|pi)\) |
| \(v_{\mathrm{xc}}\) | spin-resolved grid functional (PBE LSD/GGA spin path) matching SCF |
| QP | Newton orbit-by-orbit per spin; **report α** HOMO–LUMO gap as `gap_ev` |

**Dual-gate identity:** when \(n_\alpha=n_\beta\) and orbitals match the
closed-shell reference, UGW QP gap ≡ RKS G0W0 gap (`periodic-gw-si`
first-green **5.719545 eV**).

## Acceptance

1. Unit: `direct_rpa_uhf` E_c dual-gates `direct_rpa_rhf` on identical α=β tensors (1e-10).
2. Unit: UKS multi-k / OT / unconverged guards throw; RKS path unchanged.
3. Regression `reg_periodic-ugw-si`: UKS PBE Si₂, RI, all virtuals; hard gates
   (converged, finite positive QP gap, HOMO_QP < LUMO_QP) + dual-gate band
   \(|\mathrm{gap}-5.719545|\le 0.05\) eV (or tighter after first green).
4. Wall ≤ 60 s; `perf_exempt` if no fair CP2K spectral twin.
5. Offline twin input under `tools/cp2k_oracle/inputs/periodic-ugw-si.inp`.

## Non-goals

- Multi-k UGW, URPA energy product, BSE@UGW, spin-flip RPA, molecular UHF G0W0.
- Live CP2K in ctest.
- Public Config bags / factories.

## Files (expected)

- `include/atomiverse/cpu/gpw/g0w0.h` — docs for UKS
- `include/atomiverse/cpu/gpw/gpw_solver.h` — `g0w0` throws table
- `include/atomiverse/cpu/gpw/internal/mp2_kernels.h` + `src/.../mp2_kernels.cpp` — `direct_rpa_uhf`
- `src/cpu/gpw/g0w0.cpp` — UKS Γ path
- `tests/cpu/gpw/test_cpu_gpw_g0w0.cpp`
- `regressions/gpw/periodic-ugw-si/`
- `tools/cp2k_oracle/inputs/periodic-ugw-si.inp`
- `regressions/CMakeLists.txt`

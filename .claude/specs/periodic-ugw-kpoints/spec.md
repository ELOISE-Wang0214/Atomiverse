---
slug: periodic-ugw-kpoints
status: in_progress
created: 2026-07-26
domain: gpw / response
queue: cp2k-method-loop
depends_on: [periodic-ugw-si, periodic-gw-kpoints, periodic-ump2-kpoints]
---

# Spec: `periodic-ugw-kpoints` — Multi-k unrestricted spectral G0W0

## Summary

Close the **multi-k UGW residual**: product path for spin-conserving spectral
G0W0 on multi-k **GPW UKS** orbitals (finite-\(q\) inter-\(k\) \(W\)), dual-
gated against closed-shell multi-k spectral G0W0 when \(n_\alpha=n_\beta\) and
\(C^\alpha\approx C^\beta\).

Today:
- Γ UGW ✓ (`periodic-ugw-si`)
- multi-k RKS finite-q G0W0 ✓ (`periodic-gw-kpoints` / `periodic-gw-finite-q`)
- multi-k UKS stores `C_k` / `C_k_beta` (UMP2 ✓)
- `GpwDft::g0w0` **throws** on multi-k UKS ("multi-k UGW is a follow-on")

## Science

For each transfer \(q=\mathrm{wrap}(k_a-k_i)\) on the closed MP mesh:

1. Build spin-conserving dRPA over α→α ⊕ β→β transitions \(i_k\to a_{k+q}\)
   (Coulomb-only Casida; spins explicit — **no** closed-shell factor 2):
   \(A=\Delta\varepsilon\,\delta+V\), \(B=V\).
2. Form spectral transition potentials \(W_s^q=\sum_{ia\sigma}(X+Y)_{ia\sigma}^s V_{ia\sigma}^q\)
   (**no** \(\sqrt{2}\) — spins explicit; matches Γ UGW).
3. Same-spin spectral \(\Sigma_c^\alpha\): contract \(W_s\) with α intermediates
   at \(k-q\) with weight \(1/N_k\).
4. Per-\(k\) MT \(\Sigma_x^\alpha\) and spin-resolved
   \(\langle p_\alpha|v_{\mathrm{xc}}^\alpha|p_\alpha\rangle\) on the
   BZ-weighted UKS density.
5. Newton QP orbit-by-orbit on the α spectrum. Report **min direct α QP gap**
   over the mesh (primary report dual-gates RKS when occupations match).

Recovers Γ UGW science at \(N_k=1\) (real-valued complex MOs). When
\(C^\alpha=C^\beta\), multi-k UGW ≡ multi-k RGW (unit dual-gate).

## Public API

```cpp
cpu::gpw::GpwDft m;
m.n_alpha = 8; m.n_beta = 8;  // UKS
m.k_mesh = {2,1,1}; m.k_time_reversal = false;
m.mp2_ri = true; m.mp2_n_virt = 2;
m.build(atoms); m.run(atoms);
const cpu::gpw::G0W0Result& gw = m.g0w0(atoms);  // multi-k UGW, no throw
// gw.n_kpts == 2; gw.gap_ev = min_k direct α QP gap [eV]
```

No public `*Config` bags; no `make_*` factories. Parameters remain members.

## Acceptance

1. `GpwDft::g0w0` no longer throws on multi-k UKS with `C_k`/`C_k_beta`.
2. Unit: `spectral_g0w0_uhf_kpoints` at \(N_k=1\), \(C^\alpha=C^\beta\) recovers
   `spectral_g0w0_rhf_kpoints` QP to 1e-8 Ha.
3. Regression `regressions/gpw/periodic-ugw-kpoints/`: diamond Si bct 4-atom,
   UKS PBE \(n_\alpha=n_\beta=8\), MP 2×1×1, n_virt=2, dual-gate band vs
   `periodic-gw-kpoints` first-green min-direct QP gap (7.529834 eV, tol 0.5 eV).
4. Offline CP2K input twin under `tools/cp2k_oracle/inputs/` (never live in ctest).
5. `perf_exempt` (CP2K RI-G0W0 ≠ spectral + MT; no fair wall twin).
6. Wall ≤ 60 s on product deck.

## Out of scope

URPA energy path; spin-polarized (open-shell) literature pin with \(n_\alpha\neq n_\beta\);
frequency-integrated CP2K RI-GW absolute digit twin; CUDA.

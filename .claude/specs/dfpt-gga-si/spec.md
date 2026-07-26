---
slug: dfpt-gga-si
status: acceptance_met
created: 2026-07-26
domain: response / planewave DFPT
queue: cp2k-method-loop
depends_on: [dfpt-phonon-si]
unblocks: []
---

# Spec: `dfpt-gga-si` — DFPT phonons with GGA XC kernel (PW)

**Slug:** `dfpt-gga-si`  
**Domain:** plane-wave density-functional perturbation theory (CPU)  
**Queue:** CP2K method loop  
**Depends on:** `dfpt-phonon-si` (Γ LSDA DFPT product)

## Summary

Extend `cpu::pw::DfptPhonon` / `internal::dfpt_*_force_constants` so the
self-consistent Sternheimer response uses a **true GGA XC kernel** (PBE and
other pure GGAs), not the local LDA `f_xc = ∂V_xc/∂n` only. Today pure GGA
ground-state SCF works on PW (`XCBuilder` spectral ∇ρ) but DFPT throws
`XCType::LSDA only`. This closes the roadmap residual
(`docs/roadmap/dft.md` §11; queue row `dfpt-gga-si`).

## Science

Baroni et al., *Rev. Mod. Phys.* **73**, 515 (2001). The response potential is

\[
\Delta V_{\mathrm{Hxc}}[\Delta n]
  = \Delta V_{\mathrm{H}}[\Delta n]
  + K_{\mathrm{xc}}[\Delta n],
\]

with \(\Delta V_{\mathrm{H}}(G)=4\pi\Delta n(G)/G^2\). For pure LDA,
\(K_{\mathrm{xc}}[\Delta n]=f_{\mathrm{xc}}(r)\,\Delta n(r)\). For GGA the
kernel is the second functional derivative of \(E_{\mathrm{xc}}[n]\) and
includes density-gradient terms (flux / divergence).

**Implementation choice (discrete consistency):** apply \(K_{\mathrm{xc}}\) by
central finite difference of the **same** `XCBuilder::step` operator the
ground-state SCF uses:

\[
K_{\mathrm{xc}}[\Delta n]
  \approx
  \frac{
    V_{\mathrm{xc}}[n_0+\varepsilon\Delta n]
    - V_{\mathrm{xc}}[n_0-\varepsilon\Delta n]
  }{2\varepsilon},
\]

with sphere filtering identical to the ground-state \(V_{\mathrm{xc}}\) path.
This makes DFPT the exact second derivative of the *discretized* KS energy
(same rationale as the LSDA energy-FD twin in `reg_dfpt-phonon-si`). Pure LDA
keeps the cheap local `lda_fxc_at` path. Hybrids (`hybrid_exchange_fraction ≠ 0`)
remain unsupported (no exact-exchange response).

## Public API

- Unchanged surface: `cpu::pw::DfptPhonon` members + `compute(method, frame)`.
- `method.xc` may be pure LSDA **or** pure GGA (`PBE`, `PBEsol`, `RevPBE`,
  `RPBE`, `BLYP`). Still require UPF, fixed occupations, restricted spin.
- No public `*Config`, no `make_*` / `build_*`.

## System (regression twin)

| Setting | Value |
| ------- | ----- |
| Structure | diamond Si, 2-atom fcc primitive, \(a=5.43\) Å |
| PP | in-tree `tests/data/pw/Si_lda_nc.upf` (NC ONCV; same GGA-smoke PP as `pw-pbe-gga-si`) |
| XC | `XCType::PBE` (true GGA) |
| Cutoff / k | 10 Ry, \(2\times2\times2\) MP (minimal twin; match LSDA DFPT budget) |
| Gate | offline energy second-difference of SCF total energy for \(C_{(0x)(0x)}\); optical window; ASR; acoustic bound |

## Out of scope

- Hybrid / meta-GGA DFPT; spin-polarized DFPT; metallic smearing.
- Live CP2K in ctest (offline pin only).
- Finite-q GGA regression (Γ product path is the acceptance; q-path gets
  the same kernel so LSDA-only does not regress).

## Acceptance

See `acceptance.md`. Green ⇔ unit guards + `reg_dfpt-gga-si` pass; energy-FD
agreement documented in `REFERENCE.md`; `perf_exempt` if no CP2K binary.

---
slug: dfpt-spin-si
status: acceptance_met
created: 2026-07-26
domain: response / planewave DFPT
queue: cp2k-method-loop
depends_on: [dfpt-phonon-si, dfpt-gga-si]
unblocks: [dfpt-metallic-smear]
---

# Spec: `dfpt-spin-si` — DFPT phonons with collinear spin (PW)

**Slug:** `dfpt-spin-si`  
**Domain:** plane-wave density-functional perturbation theory (CPU)  
**Queue:** CP2K method loop  
**Depends on:** `dfpt-phonon-si` (Γ LSDA), `dfpt-gga-si` (pure GGA kernel)

## Summary

Extend `cpu::pw::DfptPhonon` / `internal::dfpt_gamma_force_constants` so
collinear spin (`PW::spin_polarized`) is supported for fixed-occupation
pure LSDA / pure GGA. Today DFPT throws `collinear spin unsupported`.
This closes the wave-7 queue residual and the roadmap §11 note that DFPT
is closed-shell only.

## Science

Baroni et al., *Rev. Mod. Phys.* **73**, 515 (2001) collinear spin:

- Ground-state SCF is collinear UKS (`nspin=2`), fixed occupations
  \(N_\uparrow-N_\downarrow=\mathrm{multiplicity}-1\).
- Response densities \(\Delta n_\uparrow,\Delta n_\downarrow\) with
  \(f=1\) per spin-orbital (not \(f=2\)).
- Hartree from total \(\Delta n=\Delta n_\uparrow+\Delta n_\downarrow\).
- Spin XC kernel \(K_{\mathrm{xc}}^{\sigma\sigma'}[\Delta n]\): local
  \(f_{\mathrm{xc}}\) matrix for pure LDA; central FD of
  `XCBuilder::step_spin` for pure GGA (same discrete operator as SCF).
- Sternheimer per spin channel with frozen \(H_\sigma=T+V_{\mathrm{loc}}+V_H+V_{\mathrm{xc},\sigma}\).

Hybrids and metallic smearing remain unsupported (follow-on).

## Public API

- Unchanged surface: `DfptPhonon` members + `compute(method, frame)`.
- `method.spin_polarized = true` accepted when `smearing == none` and XC is
  pure LSDA/GGA; `method.multiplicity` sets \(2S+1\).
- No public `*Config`, no `make_*` / `build_*`.

## System (regression twin)

| Setting | Value |
| ------- | ----- |
| Structure | diamond Si, 2-atom fcc primitive, \(a=5.43\) Å |
| PP | in-tree `tests/data/pw/Si_lda_nc.upf` |
| XC | `XCType::LSDA` (collinear LSDA kernel) |
| Spin | `spin_polarized=true`, `multiplicity=1` (collinear framework, \(M=0\)) |
| Cutoff / k | 10 Ry, \(2\times2\times2\) MP (match LSDA DFPT budget) |
| Gate | offline energy second-difference of **spin** SCF total energy for \(C_{(0x)(0x)}\); optical window; ASR; acoustic bound |

Closed-shell collinear (\(M=0\)) on Si is the minimal twin that exercises
true \(nspin=2\) machinery while remaining FD-pinnable against the same
Hamiltonian. Net-moment systems are a follow-on, not this slug's acceptance.

## Out of scope

- Metallic smearing DFPT (queued behind this product).
- Noncollinear spin; hybrids; finite-q spin product regression (Γ is the
  product path; commensurate finite-q reuses Γ-supercell and inherits spin).
- Live CP2K in ctest.

## Acceptance

See `acceptance.md`.

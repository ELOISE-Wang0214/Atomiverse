# Spec — dfpt-finite-q-spin-smearing

**Slug:** `dfpt-finite-q-spin-smearing`  
**Method:** Finite-q DFPT collinear spin + metallic smear (PW)  
**Domain:** response / planewave  
**Oracle:** CP2K (offline only; GTH-GPW twin documents capability; absolute pins are NC-UPF same-Hamiltonian)  
**Branch:** `wf/dfpt-finite-q-spin-smearing`  
**Parent:** `dev`

## Gap

On `dev`:

| Path | Fixed occ | Restricted smear | Collinear spin | Spin + smear |
|------|-----------|------------------|----------------|--------------|
| Γ DFPT | ✓ | ✓ | ✓ | ✓ (`dfpt-spin-smearing`) |
| Commensurate finite-q (supercell Γ) | ✓ | ✓ | ✓ | **claimed in docs; unpinned combo** |
| Direct finite-q Sternheimer | ✓ | ✓ | ✓ | **throws** |

`dfpt_q_force_constants` rejects `spin_polarized && smearing != none`:

```
collinear spin + metallic smearing unsupported (restricted smear only)
```

Unit gate: `DfptPhonon.ThrowsOnSpinSmearAtIncommensurateQ`.

## Product target

Enable **collinear spin + metallic occupation smearing** at finite \(\mathbf{q}\):

1. **Commensurate** \(\mathbf{q}\) (e.g. \(X=(1/2,0,0)\)): supercell Γ-DFPT (already spin+smear-capable) + IFC Fourier — ship regression twin.
2. **Direct** finite-\(q\) Sternheimer (`dfpt_q_spin_impl`): fractional spinor \(f_{nk\sigma}\in[0,1]\) from shared \(\mu\), projector on \(f\ge 10^{-8}\) at \(\mathbf{k}+\mathbf{q}\); **no** de Gironcoli \(\delta\mu\) at \(q\neq 0\) (same as restricted finite-q smear).

Public surface: `cpu::pw::DfptPhonon` members + `PW::{spin_polarized,smearing,degauss,...}` — no `*Config`, no factories.

## Acceptance

| Gate | Criterion |
|------|-----------|
| Unit | Direct finite-q spin+smear reaches UPF I/O (no spin+smear gate throw) |
| Regression | `regressions/response/dfpt-finite-q-spin-smearing` — diamond Si at X, collinear M=1, FD degauss 0.02 Ry, LSDA, 10 Ry, k=1 on SC |
| Spectra | \(\lvert\omega\rvert\) windows track fixed-occ / smear-only siblings (gapped FD fill) |
| Perf | `perf_exempt` if no CP2K binary on host; else \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\le 10\) |
| Wall | regression ≤ 60 s |

## Non-goals

- Hybrid XC DFPT; cold / Marzari–Vanderbilt smear enum on this pin path
- GPW DFPT; LO–TO / Born charges; CUDA DFPT
- Full self-consistent Mermin GS density loop inside DFPT (frozen integer-fill GS + smear on frozen spectrum, matching Γ spin+smear product)

## Files (expected)

- `src/cpu/planewave/internal/dfpt.cpp` — `dfpt_q_spin_impl` smear occupations; drop spin+smear throw
- `include/atomiverse/cpu/planewave/dfpt_phonon.h`, `internal/dfpt.h` — docs
- `tests/cpu/planewave/test_dfpt.cpp` — unit gate flip
- `regressions/response/dfpt-finite-q-spin-smearing/` + CMake
- `tools/cp2k_oracle/inputs/dfpt-finite-q-spin-smearing.note`

## Status

- [x] acceptance_met (2026-07-26): unit `test_cpu_pw_dfpt` PASS; reg
  `reg_dfpt-finite-q-spin-smearing` PASS (~11 s); spectra match fixed-occ
  spin incommensurate twin; `perf_exempt` (no CP2K binary).

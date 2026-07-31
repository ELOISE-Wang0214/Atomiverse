---
slug: dfpt-gga-spin-si
status: approved
created: 2026-07-31
revised: 2026-07-31
domain: pw / response (DFPT phonons)
queue: cp2k-method-queue
depends_on: [dfpt-gga-si, dfpt-spin-si]
unblocks: [dfpt-gga-finite-q-si]
grilled: false
---

# Spec: `dfpt-gga-spin-si` — Γ DFPT phonons, joint collinear spin × pure GGA

## Problem

`cpu::pw::DfptPhonon` ships two *separate* pinned products:

- `regressions/response/dfpt-gga-si` — closed-shell restricted response with a
  **true PBE GGA kernel** (FD of `XCBuilder` on the total density).
- `regressions/response/dfpt-spin-si` — **collinear spin** response
  (`spin_polarized = true`, multiplicity = 1) with the local **LSDA**
  \(f_{\mathrm{xc}}^{\sigma\sigma'}\) table.

The joint path — spin-resolved Sternheimer response where
\(K_{\mathrm{xc}}^{\sigma\sigma'}\) comes from central FD of
`XCBuilder::step_spin` with a **pure GGA** functional (PBE) — exists in
`src/cpu/planewave/internal/dfpt.cpp` (`apply_xc_kernel_spin_fd`,
`use_gga` branch of `dfpt_gamma_spin_impl`) but is **unpinned**: no oracle
number gates it, no regression exercises it. Per the loop rule, a foothold
is not a closed capability.

## Product

One regression `regressions/response/dfpt-gga-spin-si/` (family home of all
`dfpt-*` cases) exercising the **public path only**:

```cpp
cpu::pw::PW pw;                 // upf_path, ecut_ry = 10, k_mesh = 2
pw.xc = cpu::scf::XCType::PBE;  // pure GGA kernel
pw.spin_polarized = true;       // collinear, multiplicity = 1 (M = 0)
cpu::pw::DfptPhonon dfpt;       // Γ (q_frac = 0)
auto r = dfpt.compute(pw, frame);
```

Geometry: diamond Si, 2-atom fcc primitive cell, a = 5.43 Å, in-tree
`tests/data/pw/Si_lda_nc.upf` — the exact regression twin of both siblings so
the joint numbers are comparable against each parent.

No new public API. No behavior change expected; if the joint kernel is wrong,
the FD gate below exposes it and the fix lands inside
`internal::` (iron law: no silent debt).

## Oracle / pins (offline; per item detail: CP2K/FD-supercell dual gate)

No CP2K module exists (`module avail cp2k` empty), but the **locally built
CP2K 2026.1** binary (`~/.local/bin/cp2k`, rev `85b8a9b49c` — the same one
the wave-14 `gpw-d3-forces` pin used) runs on this host, so this slug adds a
real CP2K anchor on top of the family's FD dual gate:

1. **Pin 1 (sharp, same-Hamiltonian):** central second difference of the
   **collinear-spin PBE** SCF total energy,
   \(C_{(0x),(0x)} = (E(+h) - 2E_0 + E(-h))/h^2\), h = 0.03 Å,
   SCF tol 1e-9, identical ecut/k/XC/UPF/spin. This is the discrete Hessian
   of the very functional the response differentiates — an independent
   energy-only code path (Davidson SCF total energy), not the DFPT solver.
2. **Pin 2 (physics window, CP2K-anchored):** optical Γ modes inside
   [450, 650] cm⁻¹ (spread ≤ 10 for ATV's degenerate triplet); CP2K 2026.1
   UKS-PBE k-point `VIBRATIONAL_ANALYSIS` on the same cell (GTH/MOLOPT twin)
   prints 463.44 / 552.19 / 552.19 cm⁻¹ inside the same window; acoustic ≈ 0
   with **no** ASR imposition; row-sum ≤ 0.05 eV/Å²; pre-symmetrization
   asymmetry ≤ 1e-5 eV/Å².
3. **Cross-parent sanity:** joint C00 and optical triplet must sit near the
   two parents (GGA closed-shell 17.3456/579.8; spin-LSDA 17.3390/579.5),
   since M = 0 spin doubling is a consistency limit, not new physics.

Pins measured on this host are recorded in `REFERENCE.md` +
`tools/cp2k_oracle/inputs/dfpt-gga-spin-si.note` (re-runnable recipe with
the CP2K input under `inputs/dfpt-gga-spin-si/`).
ctest never launches CP2K; the regression hard-codes the numbers.
Fair wall: ATV 44.7 s vs CP2K 221.8 s (serial) → ratio 0.20, `within_10x`.

## Acceptance

See `acceptance.md`. Green = regression builds under `ATV_BUILD_PW` +
`ATV_BUILD_REGRESSIONS`, passes all five gates, wall ≤ 60 s.

## Out of scope

- Finite-q spin×GGA (`dfpt-gga-finite-q-si` is a held sibling row).
- Hybrid functionals (rejected by `DfptPhonon` by design).
- M ≠ 0 open-shell cells (no in-tree UPF twin; separate row if queued).
- CUDA.

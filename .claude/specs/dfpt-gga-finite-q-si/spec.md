---
slug: dfpt-gga-finite-q-si
status: approved
created: 2026-07-31
revised: 2026-07-31
domain: pw / response (DFPT phonons)
queue: cp2k-method-queue
depends_on: [dfpt-finite-q-si, dfpt-gga-si, dfpt-gga-spin-si]
unblocks: [dfpt-gga-spin-finite-q (census-held joint product)]
grilled: false
---

# Spec: `dfpt-gga-finite-q-si` — finite-q DFPT phonons with a pure GGA (PBE) kernel

## Problem

`cpu::pw::DfptPhonon` ships two *separate* pinned finite-q / GGA products:

- `regressions/response/dfpt-finite-q-si` — finite \(\mathbf{q}=X\) via the
  commensurate supercell Γ-DFPT + IFC Fourier path, **LSDA only**.
- `regressions/response/dfpt-gga-si` — **true PBE GGA kernel** (central FD of
  the `XCBuilder` operator), **Γ only**.

The joint path — finite-\(\mathbf{q}\) Sternheimer response where
\(K_{\mathrm{xc}}\) is the PBE GGA kernel — exists in
`src/cpu/planewave/internal/dfpt.cpp` (`dfpt_gamma_force_constants` accepts
pure GGA and the commensurate path reuses it on the supercell;
`dfpt_q_force_constants` carries the same `use_gga` Re/Im FD branch) but is
**unpinned**: no oracle number gates the finite-q GGA product, no regression
exercises it. Per the loop rule, a foothold is not a closed capability.
GapCensus wave 13 opened this row; the hold ("shares the DFPT GGA kernel with
`dfpt-gga-spin-si`") cleared when that pick merged in shift harness wave 1.

## Product

One regression `regressions/response/dfpt-gga-finite-q-si/` (family home of
all `dfpt-*` cases) exercising the **public path only** — the
`dfpt-finite-q-si` twin with `XCType::PBE`:

```cpp
cpu::pw::PW pw;                 // upf_path, ecut_ry = 10, k_mesh = 1
pw.xc = cpu::scf::XCType::PBE;  // pure GGA kernel
cpu::pw::DfptPhonon dfpt;
dfpt.q_frac = {0.5, 0.0, 0.0};  // X → commensurate 2×1×1 supercell Γ-DFPT
auto r = dfpt.compute(pw, frame);
```

Geometry: diamond Si, 2-atom fcc primitive cell, a = 5.43 Å, in-tree
`tests/data/pw/Si_lda_nc.upf` — the exact regression twin of the LSDA finite-q
parent so the PBE numbers are directly comparable against it.

No new public API. No behavior change expected; if the finite-q GGA kernel is
wrong, the energy-FD gate below exposes it and the fix lands inside
`internal::` (iron law: no silent debt).

## Oracle / pins (offline; per item detail: CP2K supercell twin + sharp FD pin)

CP2K 2026.1 (`~/.local/bin/cp2k`, rev `85b8a9b49c` — the binary the
`dfpt-gga-spin-si` pin used) has no PW-DFPT, so the twin is
`RUN_TYPE VIBRATIONAL_ANALYSIS` on the **2×1×1 supercell** (4 atoms): the
supercell Γ FD spectrum folds primitive Γ ∪ X, so its X-branch values anchor
the finite-q physics window on the same restricted-PBE Hamiltonian class.

1. **Pin 1 (sharp, same-Hamiltonian):** central second difference of the
   restricted **PBE** PW SCF total energy **on the 2×1×1 supercell**,
   \(C_{(0x),(0x)} = (E(+h) - 2E_0 + E(-h))/h^2\), h = 0.03 Å, SCF tol 1e-9,
   identical ecut/k/XC/UPF. This is the discrete Hessian of the very
   functional the supercell response differentiates — an independent
   energy-only code path (Davidson SCF total energy), not the DFPT solver.
   Gate at the `dfpt-gga-spin-si` level (|Δ| ≤ 0.30 eV/Å² class).
2. **Pin 2 (physics window, CP2K-anchored):** ATV |ω(X)| windows anchored by
   the CP2K supercell FD X-branch values (input + output line committed under
   `tools/cp2k_oracle/inputs/dfpt-gga-finite-q-si/`), widened for the known
   minimal-twin bias (10 Ry NC-UPF LDA-generated PP, k=1 supercell — the LSDA
   parent shows the same offset; windows finalized in `REFERENCE.md`).
3. **Cross-parent sanity:** the PBE finite-q spectrum must sit near the LSDA
   parent (`dfpt-finite-q-si`: |ω| = 130, 130, 230, 448, 673, 673 cm⁻¹) —
   XC-kernel change, not new lattice dynamics — and the supercell Hessian
   symmetry/row-sum diagnostics stay at the family level.

Pins measured on this host are recorded in `REFERENCE.md` +
`tools/cp2k_oracle/inputs/dfpt-gga-finite-q-si.note` (re-runnable recipe).
ctest never launches CP2K; the regression hard-codes the numbers.
Fair wall: ATV regression vs the CP2K supercell VIB run (serial), require
ratio ≤ 10 or `perf_exempt` with a stated reason.

## Acceptance

See `acceptance.md`. Green = regression builds under `ATV_BUILD_PW` +
`ATV_BUILD_REGRESSIONS`, passes all gates, wall ≤ 60 s.

## Out of scope

- Joint spin × GGA × finite-q (census-held behind this row).
- Incommensurate q with GGA (direct-Sternheimer GGA path is shared code; a
  dedicated pin is a separate census decision).
- Hybrid functionals (rejected by `DfptPhonon` by design).
- CUDA.

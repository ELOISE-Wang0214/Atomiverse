# Acceptance — `gpw-pbesol-si`

All criteria must be verified (`acceptance_met=true` only then).

## A1 — Public API

- [x] `XCType::PBEsol` exists; `is_gga(PBEsol)==true`;
      `hybrid_exchange_fraction(PBEsol)==0`; `hybrid_pt2_fraction(PBEsol)==0`.
- [x] `XcKind::PBEsol` exists and is accepted by `GpwDft` / `GPWEvaluator`
      (RKS).
- [x] No new public `*Config` bag; no `make_*` / `build_*` factory.
- [x] Si resolves: `GthPotential(Element::Si, "GTH-PBE-q4")` with zion 4;
      `GaussianBasis(Element::Si, "DZVP-MOLOPT-GTH")` non-empty sets.

## A2 — Unit (tests/)

- [x] Pointwise: at \(\sigma=0\), `exc(PBEsol)` matches `exc(PBE)` to 1e-12
      (same LDA base + Fx=1).
- [x] Pointwise: at finite \(\sigma\), `|exc(PBEsol)-exc(PBE)| > 1e-8`.
- [x] GPW grid: constant density → PBEsol energy matches volume × ρ × exc
      (LDA limit); matches PBE@σ=0 on grid.
- [x] Table unit: Si GTH rloc / first h and MOLOPT nexp pinned to CP2K digits.

## A3 — Regression (e2e)

- [x] `regressions/gpw/gpw-pbesol-si/reg_gpw-pbesol-si.cpp` public headers only.
- [x] SCF converges (`last_result().converged`).
- [x] Dual-gate: `|E_PBEsol - E_PBE| > 1e-4` Ha (functional activation;
      measured \(2.011\times10^{-2}\) Ha).
- [x] Offline pin on \(E_{\mathrm{PBEsol}}\) [Ha] within stated tol in
      REFERENCE.md (dual-gate + frozen ATV first-green −13.0542687637 Ha;
      no CP2K binary — never invent CP2K digits from memory).
- [x] ≤ 60 s single-rank Release (~0.7 s OMP=4); `ATV_BUILD_GPW=ON`.
- [x] No live third-party spawn.

## A4 — REFERENCE + oracle

- [x] `REFERENCE.md` cites functional, basis/PP, cell, cutoff, pin strategy,
      tol, perf section.
- [x] Twin input `tools/cp2k_oracle/inputs/gpw-pbesol-si.inp` present for
      offline re-pin.

## A5 — Perf gate

- [x] `perf_exempt=true` with reason (no CP2K binary on host).

## A6 — Build / commit

- [x] Builds with HF+PW+GPW, CUDA OFF.
- [x] Commit on `wf/gpw-pbesol-si` only (never `dev`); no push.

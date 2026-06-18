---
title: Native XC functional library with libxc (and CP2K-native) bit-parity
slug: native-xc-libxc-parity
status: draft
science_required: true
acceptance: native-xc-libxc-parity.acceptance.md
---

## Summary

Build a **native exchange–correlation library** in Atomiverse (`cpu/xc/`) whose
pointwise output (`exc`, `vrho`, `vsigma`, and higher derivatives where needed)
reproduces **libxc bit-for-bit** for every functional CP2K routes through libxc,
and reproduces **CP2K's own native formula** for the few it does not (notably
`PADE`). Validation is the core deliverable: the installed **libxc.so.15** (already
present in the Spack view, headers `xc.h`/`xc_funcs.h`) is used as a test-only
oracle, and each functional is swept over a dense (ρ, ζ, ∇ρ) grid asserting machine
agreement.

This library replaces design item **D7** of [[gpw-gaussian-plane-wave]] (originally
"link libxc at runtime"). It is consumed by *both* electronic-structure paths
(`cpu::scf` and the new `cpu::gpw`) through the existing `evaluate_xc` /
`evaluate_xc_batch` interface, so eliminating the XC mismatch removes one of the
largest sources of divergence between Atomiverse and CP2K.

**End goal:** pass the tutorial examples and their associated regression tests —
which requires `PADE` (CP2K-native LDA) for the first GPW examples and `PBE`
(libxc-routed GGA) for the getting-started tutorial.

## Domain basis (why bit-parity is *achievable* here)

Unlike CP2K total-energy parity (bounded at ~1e-6 Ha by grid/FFT discretization —
see [[gpw-gaussian-plane-wave]]), **XC parity is reachable to machine precision**:
an XC functional is deterministic pointwise math in IEEE float64. If we evaluate the
identical analytic expression libxc evaluates, the only residual is the **last-ulp
behaviour of libm transcendentals** (`pow`, `exp`, `log`, `cbrt`, `atan`). Target:

| Quantity | Acceptance vs libxc oracle |
|---|---|
| `exc`, `vrho`, `vsigma` (and `v2*` where built) | **max rel. error ≤ 1e-12** over the sweep |
| residual cause if any point exceeds | documented last-ulp libm difference, bounded ≤ few ulp |

libxc's per-functional math is Maple-generated from a published energy expression;
the SCF/GPW consumer needs `exc` + first derivatives (the potential). Second
derivatives (`v2rho2`, `v2rhosigma`, `v2sigma2`) are only needed for response
properties (TDDFT/linear-response) — out of scope until a tutorial requires them.

## Design

### Two source tracks (vendored)

- **libxc-routed functionals** (PBE and, later, most GGAs/hybrids): CP2K evaluates
  these via `src/xc/xc_libxc.F`, so *matching libxc == matching CP2K*. **Vendor the
  libxc per-functional source** (the maple-generated `src/*.c` kernels +
  `work_lda.c`/`work_gga.c` machinery) and adapt to standalone compilation — exact
  by construction. The oracle confirms the adaptation is wired correctly.
- **CP2K-native functionals** (`PADE`): CP2K implements these in `src/xc/xc_pade.F`,
  NOT through libxc. **Port `xc_pade.F` Fortran→C++** to match CP2K's `Ar.inp`/
  `H2.inp` exactly (libxc's `LDA_XC_TETER93` is close but differs — do not
  substitute it). Validate against a CP2K run (compare the `E_xc` term).

### Scope is tutorial/test-driven, not "all 600 functionals"

Implement only what the milestone inputs exercise, in order:

- **X0 — Oracle harness.** Wire libxc as a *test-only* oracle (the
  `ATV_BUILD_XC_ORACLE` flag already exists). A sweep generator produces points
  across ρ ∈ [1e-10, 1e3], spin ζ ∈ [0,1] (polarized + unpolarized), reduced
  gradient s ∈ [0, ~5] incl. edge cases (ρ→0, fully polarized, σ=0 and large σ),
  and compares native vs `xc_lda_exc_vxc` / `xc_gga_exc_vxc`.
- **X1 — LDA building blocks.** Vendor libxc `LDA_X`, `LDA_C_PW`, `LDA_C_VWN`
  kernels. (Atomiverse's current *native* LSDA in `cpu/scf/xc_functional.cpp` is
  retired in favour of the vendored versions so values match CP2K's libxc-routed
  LDA exactly.)
- **X2 — PADE (CP2K-native).** Port `src/xc/xc_pade.F` → unlocks **M1** of
  [[gpw-gaussian-plane-wave]] (`Ar.inp` −21.04944231395054 Ha, `H2.inp`
  −1.06345574964568 Ha). Validated against CP2K.
- **X3 — PBE.** Vendor libxc `GGA_X_PBE` + `GGA_C_PBE` (spin-polarized) → unlocks
  **M3**, the getting-started `h2o.inp` tutorial. Exact against libxc.
- **X4+ — on demand.** Add functionals as further tutorials/regtests require
  (e.g. BLYP, SCAN, PBE0/B3LYP exchange fractions already present natively).

### Architecture

- New shared module `cpu/xc/` (promoted out of `cpu/scf/`) so `cpu::scf` and
  `cpu::gpw` both depend on XC without GPW pulling in all of SCF. Keep the existing
  public signature `evaluate_xc(XCType, ρ_α, ρ_β, σ_αα, σ_αβ, σ_ββ) -> {exc, vrho_α,
  vrho_β, vsigma_*}` and the batched `evaluate_xc_batch` (current callers in
  `hamiltonian_builder.cpp` / `xc_builder.cpp` unchanged).
- Internally, a **functional registry**: each functional is a pure pointwise kernel
  with the libxc functional ID attached (for the oracle to look up the counterpart).
  Composite/hybrid types compose X + C kernels + an exact-exchange fraction
  (`hybrid_exchange_fraction` already exists).
- Spin path: implement the **polarized** kernel as the primary; unpolarized is the
  ζ=0 special case (matches how libxc exposes `XC_POLARIZED`/`XC_UNPOLARIZED`).
- Determinism: pointwise, no reductions → bit-reproducible regardless of threads.

### Licensing — DECIDED: vendor source verbatim, Atomiverse relicenses GPL-2.0

Operator decision (2026-06-18): **copy source verbatim, license is not a
constraint.** Consequences, recorded so they are deliberate:

- **libxc kernels are vendored as-is** (MPL-2.0). libxc already ships the
  maple-generated per-functional `.c` files in its `src/` — those are copied
  directly, so agreement with libxc is *exact by construction* (same source, same
  float64 math), not merely ≤ 1e-12. Keep the upstream license headers.
- **PADE is ported from CP2K** `src/xc/xc_pade.F` (GPL-2.0), translated Fortran→C++
  faithfully.
- **Atomiverse's combined binary therefore becomes GPL-2.0** (GPL is the most
  restrictive of the inputs; MPL-2.0 is GPL-compatible via its secondary-license
  clause). Update the project `LICENSE`/`README` accordingly as part of X4. This is
  an accepted, intentional relicense — flag it in the PR that lands the first
  vendored kernel.

Because sources are vendored, the libxc oracle (`ATV_BUILD_XC_ORACLE`) becomes a
*regression guard* (it should pass trivially) rather than the thing that earns
parity. It still runs to catch wiring/units/spin-channel bugs in the adaptation.

## Files

New:
- `include/atomiverse/cpu/xc/{xc_functional.h,xc_registry.h,functionals/*.h}`
- `src/cpu/xc/{xc_functional.cpp, functionals/{lda_x,lda_c_pw,lda_c_vwn,pade_gth,gga_x_pbe,gga_c_pbe}.cpp}`
- `tests/cpu/xc/{test_cpu_xc_oracle_sweep.cpp, test_cpu_xc_<functional>.cpp}`
- `tests/cpu/xc/oracle/` — libxc oracle wrapper (built only under `ATV_BUILD_XC_ORACLE`)

Edit / migrate:
- Move `cpu/scf/xc_functional.*`, `cpu/scf/xc_gga.*` content into `cpu/xc/`
  (keep thin compatibility include in `cpu/scf/` during transition).
- `CMakeLists.txt` — `cpu/xc/` sources compiled whenever RHF or GPW is on; ensure
  `ATV_BUILD_XC_ORACLE` links libxc for tests only (find via Spack view path).
- `.claude/architecture.md` — register `cpu/xc/` as a shared leaf module.

## Tasks

- [ ] **X0:** oracle wrapper + sweep generator (regression guard).
- [ ] **X1:** vendor libxc `LDA_X`, `LDA_C_PW`, `LDA_C_VWN` + adapt; sweep passes.
- [x] **X2:** port `xc_pade.F`→C++ (`cpu/xc/pade`); FD-validated potential, exc
      matches the Padé formula + CP2K's −0.8097 at ρ=1 (`test_cpu_xc_pade`).
      Full bit-parity vs CP2K `Ar.inp`/`H2.inp` deferred to GPW M1.
- [ ] **X3:** vendor libxc `GGA_X_PBE` + `GGA_C_PBE` polarized; sweep passes.
- [ ] **X4:** module at `cpu/xc/`, repoint `cpu::scf` + `cpu::gpw` callers;
      relicense Atomiverse to GPL-2.0 (`LICENSE`/`README`) in this step.
- [ ] **X5:** wire into GPW M1 (PADE) and M3 (PBE) parity gates of
      [[gpw-gaussian-plane-wave]].

## Testing

- Per-functional sweep tests (`tests/cpu/xc/`) vs libxc oracle, both spin channels,
  edge cases enumerated above, ≤ 1e-12 on `exc` + every first derivative.
- PADE: dedicated test comparing the XC energy to a CP2K reference produced by
  running `Ar.inp`/`H2.inp` with the local `cp2k` (`~/.local/bin/cp2k`).
- Optional: vendor a slice of libxc's own regression reference data for a second,
  build-independent check.
- Integration: the GPW/SCF `E_xc` term on the tutorial geometry matches CP2K.

## Out of scope

- Second/third derivatives (response/TDDFT), meta-GGAs beyond what a tutorial needs,
  range-separated kernels, the full 600-functional catalogue, and a CUDA XC port —
  each a later spec. Implement strictly what the tutorials and their regtests touch.

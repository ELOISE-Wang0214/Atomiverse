---
title: f-shell Rys ERI support — raise CPU Rys engine to l=3
status: approved
created: 2026-05-20
---

# f-shell Rys ERI support — raise CPU Rys engine to l=3

## Summary

The CPU Rys-quadrature electron-repulsion-integral (ERI) engine is currently capped at l≤2 (d-shells): `kRysMaxL=2`, `kRysMaxRoots=5`, `kMaxL=2`. This blocks every basis set that carries f functions — including def2-TZVP, whose checked-in table (`def2_tzvp.h`) already contains l=3 shells for C, N, O, F. This spec raises the engine to l=3 (f-shells) so that `build_integrals` succeeds end-to-end on def2-TZVP molecules. The work has two parts: a set of constant/table-sizing bumps (`kRysMaxL`, `kRysMaxRoots`, `kMaxL`, and the Hermite/Boys/E-coefficient buffers in `native_integrals.cpp`), and a genuinely new algorithm port — libcint's iterative Rys root finder for nroots 6 and 7, which the (f,f|f,f) quartet and its first derivative require. After this link, def2-TZVP and heavier-element basis sets become usable for CPU SCF; g-shells (l=4), the CUDA Rys path, and the downstream KS-gradient / denser-grid sub-specs are out of scope.

## Domain basis

Rys quadrature evaluates the ERI by an N-point Gaussian quadrature on the Rys polynomial roots over t∈[0,1]:

- **Number of roots:** N = floor(L_tot / 2) + 1, where L_tot = la + lb + lc + ld is the total angular momentum of the shell quartet. An N-point Rys quadrature is exact for the degree-L_tot polynomial-in-t² ERI integrand (exact to polynomial degree 2N−1). Reference: Dupuis, Rys & King, *J. Chem. Phys.* 65, 111 (1976), DOI 10.1063/1.432807 (King–Dupuis scheme); corroborated by libcint, Sun, *arXiv:1412.0649*.
- **(f,f|f,f):** L_tot = 3+3+3+3 = 12 → N = 7 roots.
- **Derivative of (f,f|f,f):** d/dR of a Gaussian of momentum l yields a momentum-(l+1) term (dominant for the root budget) plus a momentum-(l−1) term. One shell raised +1 gives effective (4,3|3,3), L_tot = 13 → N = 7 roots. Hence `kRysMaxRoots` must be exactly 7.
- The +1 derivative shift means the 2D Rys grid in the block-derivative engine must be extended by one angular unit on the deriv side; the existing `rys_eri_block_deriv` docstring already anticipates this — only the root count and grid dimensions need bumping, not the recurrence structure.
- **Iterative root finder (the new algorithm):** libcint's `CINTrys_roots` (src/rys_roots.c) dispatches nroots 1..5 to closed-form fits `rys_root1..rys_root5`, but switches to an iterative scheme for nroots ≥ 6 — `segment_solve(CINTrys_jacobi, CINTrys_schmidt)`: a Jacobi-polynomial three-term recurrence followed by modified Gram-Schmidt orthogonalization. High-order closed-form Chebyshev fits degrade in the intermediate-T regime, which is exactly why libcint switches schemes; the fitted-polynomial table must **not** be extended to nroots 6–7.

## Design

Three categories of change, in dependency order.

1. **Constant bumps (compile-time propagation).** `kRysMaxL` 2→3 and `kRysMaxRoots` 5→7 in `rys_eri.h`; `kMaxL` 2→3 in `cart.h` (`kMaxNCart` is `n_cart_of(kMaxL)` and follows automatically: 6→10). The stale STO-3G-specific `@note` near `rys_eri.h:96` and the analogous comment in the file header become false for any f-shell basis and must be rewritten to state the new l≤3 / 7-root budget. The `kRysMaxL`/`kRysMaxRoots` docstrings (lines 46–53) likewise need updating.

2. **Buffer sizing in `native_integrals.cpp` (72 KB).** The Hermite / Boys / E-coefficient (`md_e_coeffs`) and Rys-grid tables are already sized from `kMaxL` and `kDerivL = kMaxL + 1` (e.g. `kMaxOverlapCells`, `kMaxKineticCells`, `kERIRDim = 4*kMaxL+2`, the `kMaxNCart` ang-tables). Bumping `kMaxL` propagates these automatically, but each derived constant must be audited: the integral-derivative path raises one shell's ai by +1 (f → g-like temporary), and every stack buffer indexed by `kDerivL`, `4*kMaxL+...`, or `n_cart` must be confirmed large enough for l=3 so no fixed-size array overflows. This is the same buffer-sizing pattern applied for the d-shell derivative path in the prior phase-1 work — extend, do not restructure.

3. **Iterative Rys root finder (new code).** `src/cpu/scf/rys_roots.cpp` + `roots_for_x0.inl` currently retain only `rys_root1..5` and the small-x/large-x polyfit dispatchers; the `rys_roots()` dispatcher `abort()`s for nroots > 5. Port libcint's iterative path: a Jacobi-polynomial recurrence (`CINTrys_jacobi` equivalent) plus modified Gram-Schmidt (`CINTrys_schmidt` equivalent) and the `segment_solve` driver, namespaced into `ATV::cpu::scf::native`, covering nroots 6 and 7. The iterative root finder is a self-contained numerical routine; it may live in a new `rys_roots_iterative.{h,cpp}` unit included by `rys_roots.cpp`, or be added directly into `rys_roots.cpp` — the smaller-files preference favors a separate unit. The `rys_roots()` dispatcher gains `case 6` / `case 7` routing into it. Upstream BSD-2 attribution must be preserved as in the existing file header.

4. **libcint reference harness.** `tests/libcint_reference.cuh`'s `build_eri_libcint_from_pack` computes a per-l Cartesian normalization factor `g_l` only for l∈{0,1,2}. To cross-check an f-quartet it needs an l=3 branch (`g_3 = CINTgto_norm(3, a0) / phys_norm_f(a0)`, or the verified-empirical equivalent if `phys_norm_f` already matches libcint's Cartesian f convention). If `phys_norm_f` does not exist it must be added alongside `phys_norm_{s,p,d}`; if the harness cannot be made l=3-correct within this link, the libcint cross-check criterion degrades to a conditioned/grep-only check (see acceptance `ac-007`).

No public API signatures change: `rys_eri_block` and `rys_eri_block_deriv` keep their declarations; only their reachable l-range widens.

## Files to create or modify

- `include/atomiverse/cpu/scf/rys_eri.h` — `kRysMaxL` 2→3, `kRysMaxRoots` 5→7; rewrite stale STO-3G `@note` (~line 96), file-header comment, and the two constant docstrings.
- `include/atomiverse/cpu/scf/cart.h` — `kMaxL` 2→3 (`kMaxNCart` follows); update the `kMaxL` docstring and the `fill_cart_comps` l=2 example to mention l=3.
- `src/cpu/scf/native_integrals.cpp` — audit and resize every `kMaxL`/`kDerivL`/`4*kMaxL`/`kMaxNCart`-derived stack buffer (Hermite, Boys, `md_e_coeffs`, Rys-grid, R-table) so the l=3 + derivative-bump path does not overflow.
- `src/cpu/scf/rys_roots.cpp` — add `case 6` / `case 7` routing in `rys_roots()`; update the file-header comment that claims "only nroots ≤ 5".
- `src/cpu/scf/rys_roots_iterative.h` (new) — declaration of the iterative nroots-6/7 root finder (Jacobi recurrence + modified Gram-Schmidt), with doxygen `@brief`/`@param`(units)/`@return`.
- `src/cpu/scf/rys_roots_iterative.cpp` (new) — ported implementation, BSD-2 attribution preserved.
- `tests/cpu/scf/test_rys_eri_fshell.cpp` (new) — FD self-consistency tests for f-quartet ERIs and derivatives.
- `tests/cpu/scf/test_rys_roots.cpp` — extend with nroots-6/7 cases for the iterative finder.
- `tests/libcint_reference.cuh` — add the l=3 `g_l` branch so `build_eri_libcint_from_pack` supports f-shells.
- `include/atomiverse/cpu/scf/native_integrals.h` — add `phys_norm_f` declaration if it does not already exist and the harness needs it.
- `tests/CMakeLists.txt` — register `test_rys_eri_fshell.cpp` and the new `rys_roots_iterative` source in the SCF test target.

## Tasks

- [ ] Write failing tests for the iterative Rys root finder, nroots 6 and 7 (tests/cpu/scf/test_rys_roots.cpp) — assert roots/weights reproduce a brute-force or libcint-derived reference across small-x, intermediate-x, large-x regimes
- [ ] Port libcint's iterative Jacobi-recurrence + modified-Gram-Schmidt root finder for nroots 6 and 7 into src/cpu/scf/rys_roots_iterative.{h,cpp} (new), with doxygen and BSD-2 attribution
- [ ] Wire nroots 6/7 into the rys_roots() dispatcher in src/cpu/scf/rys_roots.cpp and update its stale "only nroots ≤ 5" header comment
- [ ] Bump kRysMaxL 2→3 and kRysMaxRoots 5→7 in rys_eri.h, kMaxL 2→3 in cart.h, and rewrite the stale STO-3G comments and constant docstrings in both headers
- [ ] Audit and resize every kMaxL/kDerivL/4*kMaxL/kMaxNCart-derived stack buffer in src/cpu/scf/native_integrals.cpp so the l=3 and derivative-bump paths do not overflow
- [ ] Write failing FD self-consistency tests for f-quartet ERIs and first derivatives — (f,s|s,s), (f,p|s,s), (f,d|s,s), (f,f|f,f) — in tests/cpu/scf/test_rys_eri_fshell.cpp (new), and register it plus the new rys_roots_iterative source in tests/CMakeLists.txt
- [ ] Add the l=3 g_l normalization branch (and phys_norm_f if missing) to tests/libcint_reference.cuh so build_eri_libcint_from_pack handles f-shells; add a libcint cross-check for one f-quartet
- [ ] Verify def2-TZVP build_integrals succeeds end-to-end on a small molecule (e.g. CH4 or H2O) with non-zero, finite ERIs
- [ ] Add docstrings per doxygen style (@brief, @param with units, @return) to every new public symbol in rys_roots_iterative.h
- [ ] Run full check + test suite

## Testing strategy

- **Happy path — root finder:** nroots-6 and nroots-7 iterative roots/weights match an independent reference (brute-force Gauss-quadrature construction or libcint values) to ≤1e-10 across small-x (x<3e-7), intermediate-x (x~1–30), and large-x (x>35+5N) regimes.
- **Happy path — f-quartet ERIs:** `rys_eri_block` produces finite, non-zero blocks for (f,s|s,s), (f,p|s,s), (f,d|s,s), (f,f|f,f).
- **Domain validation — FD self-consistency:** `rys_eri_block_deriv` for each f-quartet agrees with central finite differences of `rys_eri_block` (h=1e-4) to ~1e-7, exercising the +1 derivative angular shift and the 7-root budget.
- **Domain validation — libcint cross-check:** at least one f-quartet ERI from `rys_eri_block` matches `build_eri_libcint_from_pack` to ≤1e-9, confirming absolute correctness against the ground-truth engine (conditioned on the harness supporting l=3 after this link's harness change).
- **Integration:** def2-TZVP `build_integrals` runs to completion on a small molecule and yields a finite ERI tensor.
- **Regression:** the full `ctest` suite stays green — existing s/p/d ERI, gradient, and SCF tests show no change (the kMaxL bump must not alter any l≤2 result).
- **Edge cases:** intermediate-T regime where closed-form fits degrade (the motivation for the iterative scheme); shell-swap / non-canonical bra ordering with an f-shell.

## Out of scope

- g-shells (l=4) and beyond — `kRysMaxL` stays at 3.
- The CUDA Rys ERI path — this link is CPU-only.
- The downstream chain links: 02-ks-dshell-gradient and 03-denser-grid (both independent, neither depends on this one).
- Generating new def2-TZVP basis data — the table already exists for H, C, N, O, F with f-shells; no `gen_basis_data.py` run is required.
- Performance tuning of the iterative root finder — correctness only; optimization is deferred.

# Architecture notes

## [2026-06-10] cpu::pw migrated to PWSolver/PWEvaluator + MPI k-point parallelism; SCF audit fixes

**PW now mirrors the gaussian-SCF shape exactly.** `cpu::pw::PWSolver`
(public low-level inner solver, NOT a Driver subclass) owns the
facade-private `internal::PWWorkspace` and six step helpers in
`cpu::pw::internal::` (`PoissonSolver`, `XCBuilder`, `VeffBuilder`,
`DavidsonSolver`, `DensityBuilder`, `ConvergenceCheck`). The
`cpu::pw::PWEvaluator` Evaluator component (`Capability::Energy`)
plugs into `cpu::Evaluate`; `pw_scf()` survives as a thin shim over a
transient `PWSolver`. Forward declarations + `unique_ptr` keep
internal headers out of the public surface (same pattern as
`scf_solver.h`).

**Deleted:** `ATV::scf::Phase` (facade/stage.h), the dead
`Driver<scf::Phase, cpu::device_tag>` instantiation (src/cpu/driver.cpp),
`BuildContext::pw_workspace` (`void*` leak in component.h), the seven
Plan-shaped `PW_*` components (pw_poisson/pw_lda/pw_hamiltonian/
pw_veff/pw_davidson/pw_density/pw_convergence), and the dead
`PW_Density` duplicate — which carried a real bug: its Pulay history
was never trimmed, so after `pulay_depth` iterations it mixed against
the oldest-5 residuals forever. The consolidated `DensityBuilder` uses
a correct sliding window (`while (size > pulay_depth_) erase front`).

**Grep gate widened:** `tests/cuda/scf_phase_scrub.sh` now scans
include/, src/, examples/, tests/ and runs unconditionally as
`test_scf_phase_scrub` (was CUDA-only).

**Convention fixes folded in:** `kAngstrom2Bohr`/`kRy2Hartree` from
core/units.h (no inline magic numbers), lookup tables throw on
unsupported input (`valence_from_z` beyond Z=18, `sph_bessel_for_l`
beyond l=2 — no-fallback rule), `ATV_INFO`/`ATV_DEBUG` logging
replaces std::cout (rank-0 only), PW test targets renamed
`test_cpu_pw_*`, and a latent OOB fixed: `ws.T_G` is resized per
k-point (k-sphere G-count can exceed the Γ-sphere count).

**MPI k-point parallelism** (gated `ATV_HAS_MPI` + runtime
`MPI_Initialized`): `PWSolver::run()` distributes k-points round-robin
over `MPI_COMM_WORLD` (`ws.owned_k` / per-k `ws.gvec_k`), Allreduces
the density accumulation (`DensityBuilder::accumulation()`, MPI_SUM,
IN_PLACE) plus the band energy each SCF iteration, then finishes the
Pulay mix identically on every rank. Davidson's initial-guess noise
uses a deterministic per-subspace-index LCG (`TinyRng`, replaces
`std::rand`) so the SCF trajectory is rank-layout invariant — energy
parity across n=1/2/4 is exact, gated by `test_mpi_pw_kpoint_n{1,2,4}`
(tests/mpi/test_pw_kpoint.cpp, own MPI_Init_thread main).

**Gaussian-SCF audit fixes:** toy solver wrapped in
`cpu::scf::toy::` (ODR clash with production `SCFResult`);
`scf_solver.h` switched to forward decls + `unique_ptr` and the
isolation gate revived via `ATV_REPO_ROOT` compile definition
(`test_cpu_scf_solver_internal_isolation` actually scans the tree
again); `SCFEvaluator` pre-builds its ASPC corrector `SCFSolver` in
`build()` (no allocation in `launch()`); `unwrap_for_scf` uses
`Box::min_image` with a zero-box gas-phase guard; build() doc
corrected (only n_basis-sized prealloc is per-build; X, shell pairs,
KS grid are per-run because atoms move).

**Four pre-existing issues surfaced by the verification run:**

1. *OpenBLAS 0.3.32 / LAPACK 3.12.1 xSYGVD workspace-query regression*
   (Arrhenius `buildenv-gcc/2026.03` module — the only OpenBLAS on the
   system, all threading variants affected): the LWORK=-1 query returns
   `info=0` but writes `WORK(1)=0`, so the follow-up call aborts with
   "parameter number 11 had an illegal value". Every gaussian-SCF test
   (anything through `diagonalize_fock_generalized`) aborted. Fixed in
   `scf_utils.h` by calling `dsygvd_` directly with the documented
   workspace sizes (LWORK = 1+6N+2N², LIWORK = 3+5N) — same routine,
   same numerics, no query. A 12-line standalone repro confirmed the
   bug is environmental, not Atomiverse code.

2. *Davidson subspace-expansion heap overflow* (verbatim in the
   historical `pw_davidson.cpp`, faithfully ported): the expansion
   wrote the new column at slot `n_sub` BEFORE growing
   `sub_vectors`/`H_sub_vectors` (resize happened after the band
   loop) — one full ng-sized column OOB per added residual. Latent
   under the old allocator layout; surfaced as "double free or
   corruption" in `test_cpu_pw_si`/`test_mpi_pw_kpoint_*` after the
   refactor shifted heap layout. Located with ASan (`build-asan`,
   davidson.cpp:503); fixed by resizing before the write. NB:
   valgrind cannot run these binaries — `-march=native` emits AVX-512
   the installed valgrind doesn't decode (SIGILL at the first Mat3
   copy); use ASan on this machine.

3. *Becke weight-derivative catastrophic cancellation*
   (`src/cpu/scf/grid.cpp`, `get_smooth_function_derivative` — verbatim
   from PySCF's `grids_response_becke`): the iterated smoothing
   polynomial p3 was used as `s_uab = 0.5*(1 - p3 + 1e-200)`. For
   far-tail grid points (|mu| → 1) p3 rounds to exactly ±1, `1-p3`
   collapses to 0, and the log-derivative ratio `pt = -t/s` — whose
   mathematical divergence is cancelled downstream by multiplying with
   `pbecke ∝ s` — blows up against the 1e-200 floor instead, leaving
   ~1e150 garbage in dE_xc/dR. Reproduced 6/6 at -O0 (xc_terms
   max_abs 6.2e150, net translational force 4e149 Ha/Bohr in
   `test_cpu_scf_ks_grid_response_gradient`); -O2's FMA kept p3 at
   1-1ulp by luck, hiding the same defect. ASan/valgrind useless here
   — it is pure FP, not memory. Fixed with the cancellation-free
   recurrences `1-f(x) = (1-x)²(x+2)/2`, `1+f(x) = (1+x)²(2-x)/2`
   (s exact down to denormals, t/s a faithful log-derivative, the
   pbecke cancellation restored). After the fix: -O0 3/3 pass with
   xc_terms bitwise equal to the -O2 value.

4. *Unoptimized default build* — the documented configure line had no
   `CMAKE_BUILD_TYPE` and the top-level CMakeLists set no default, so a
   bare `cmake -B build ...` compiled everything at **-O0** (with
   `-march=native` but no -O level). Measured: the same PBE/def2-SVP
   SCF takes 97 s at -O0 vs 9.3 s in the *ASan-instrumented* -O2 tree.
   This, stacked on login-node contention, is why "second-level" tests
   ran for minutes-hours. Fixed: CMakeLists now defaults
   `CMAKE_BUILD_TYPE` to Release when unset (single-config gens), and
   the CLAUDE.md configure lines pass `-DCMAKE_BUILD_TYPE=Release`
   explicitly.

**Test-contract update — `test_cpu_scf_evaluator_no_alloc` (ac-014):**
the old gate asserted *exact* launch-to-launch alloc-count equality.
That held only while `SCFEvaluator::launch()` constructed a transient
`SCFSolver` per step (identical cold-start allocation sequence). With
the pre-built persistent solver (the no-alloc-in-launch rule), buffer
capacities and DIIS history carry across launches and OpenMP workers
grow per-thread buffers lazily, so the count wobbles ~1% (e.g.
647–657) without any leak — measured: iteration count constant at 3,
wobble persists even single-threaded (409/413). The gate now asserts
what ac-014 actually guards: bounded spread (≤5%) + no upward trend
over a 12-launch window.

**Deferred:** Phase 5b IntegralPack/ShellList hot-path split, PW
augmentation charges (E_aug = 0 for the norm-conserving UPF v1 path),
PW forces (PWEvaluator is Energy-only).

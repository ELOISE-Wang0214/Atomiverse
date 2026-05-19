# Architecture Notes

Accumulated architectural decisions and corrections from development iterations.
Agents read this alongside CLAUDE.md. When a note stabilizes, promote it into CLAUDE.md.

---

## [2026-03-25] DeviceSystem removed — metadata lives in cuda::Frame

DeviceSystem class was deleted. Host-cached metadata (n_atoms, n_bonds, n_angles, n_dihedrals) now lives directly in `cuda::Frame` as plain Int fields with accessor methods. The device-side atom count pointer is `d_n_atoms` (not `n_atoms`). (Note: `DeviceFrame` was later renamed to `cuda::Frame`.)

**Rule**: Contexts (BuildContext, LaunchContext, DestroyContext) hold `cuda::Frame& frame` directly. Components access `ctx.frame.n_atoms()` for host-cached count and `ctx.frame.fx` etc. for device pointers. No wrapper class around `cuda::Frame`.

## [2026-05-19] bonded → potentials reorg

`cpu/bonded/` retired. The three bonded evaluators
(`HarmonicBond`, `HarmonicAngle`, `MultiHarmonicDihedral`) now live in
`cpu/potentials/` alongside `LJ` and the future pair/many-body
evaluators. Reason: under the backend-root layering rule (CLAUDE.md
*Backend directory layering*), force-field evaluators are a single
domain. The split into `bonded/` (bonded interactions) vs (implicit)
`pair/` (pair potentials) reproduced the old LAMMPS pair/bond divide
without giving us anything — both produce Force/Energy/Stress on the
same Frame.

**Rule**: any new force-field evaluator (Tersoff, EAM, Morse, …)
goes into `cpu/potentials/`. Never reintroduce `cpu/bonded/`.

## [2026-05-19] PW SCF productized — solve_si retired

`cpu/pw/solve_si.{h,cpp}` and `cpu/pw/internal/crystal.{h,cpp}` were
Si-specific scaffolding from the spec-driven PW-DFT bring-up
(pw-dft-si chain). They're replaced by a generic
`pw::pw_scf(frame, upf_path, ecut, k_mesh, tol, max_iter)` in
`cpu/pw/pw_scf.{h,cpp}`, with the Python surface (`atomiverse/pw.py`
+ `bind_pw.cpp`) as the user-facing entry point. Si-specific frame
construction now belongs in user code or test fixtures, not in the
PW driver.

**Rule**: PW-DFT additions parameterize on `(frame, upf_path, ecut,
k_mesh)` — no element-specific helpers in `cpu/pw/`. Reference cases
(Si, etc.) live in tests/examples/benchmarks.

## [2026-05-19] minimizer/init.h and observer/thermo_component.h folded into set_*

Both files are gone. Their callers now use the LAMMPS-style
categorical setter API on the CPU drivers (`set_optimizer`,
`set_thermo`, etc., see CLAUDE.md *Example*). The free-function /
explicit-component idioms they exposed were duplicated by setter
auto-injection from spec
`.claude/specs/driver-infra-auto-injection.md` (now closed).

**Rule**: no new free-function "init_*" helpers or
auto-registered components in the observer/minimizer roots — they
either go through the categorical setter API or through
`driver.add()` / `driver.on(phase).add()`.

## [2026-05-19] libxc is test-only scaffolding, never a runtime dependency

libxc lives entirely under `tests/cpu/scf/oracle/` in namespace
`ATV::test::xc_oracle`. `libatomiverse.{a,so}` contains **zero** libxc
symbols regardless of build flags — verified by `nm` in the parity
test's acceptance criterion. CMake gate: `ATV_BUILD_XC_ORACLE`
(default OFF), enabling only test targets. The oracle is intentional
scaffolding: once the native XC primitive-parity gate (1e-12 vs
libxc) has stayed green long enough in CI, a "retire libxc oracle"
spec can delete `tests/cpu/scf/oracle/` and the CMake option in a
single commit, leaving production unchanged.

**Rule**: never link libxc into `libatomiverse`, never add a public
`ATV_HAS_LIBXC` macro on the production target, never put a
libxc-using function in `ATV::cpu::scf` (or any namespace under
`ATV::cpu::*`). libxc is an *oracle*, not a backend.

## [2026-05-19] PBE uses PW92 (not VWN5) as its LSDA correlation base

Native `cpu::scf::evaluate_xc(XCType::PBE, …)` calls
`pw92_correlation()` for its LDA correlation contribution. The PBE
paper (Perdew-Burke-Ernzerhof, PRL 77, 3865 (1996)) specifies PW92
(Perdew-Wang 1992), and libxc's `XC_GGA_C_PBE` follows that.
VWN5 and PW92 differ by ~mHa at typical densities; mixing them
produced a uniform ~0.4 mHa bias on every PBE energy until the spec
3 fix. **LSDA mode** continues to use VWN5 (matches libxc's
`XC_LDA_C_VWN` id 7).

**Rule**: any GGA functional that names PW92 as its LDA base
(PBE, BLYP, PW91, BP86) must call `pw92_correlation()`, not
`vwn5_correlation()`. Cross-check by reading the corresponding
libxc Maple source in `libxc/maple/`.

## [2026-05-19] Default radial grid is Treutler-Ahlrichs (M4), not Mura-Knowles

`build_molecular_grid` calls `radial_grid_treutler_ahlrichs(n, Z)` by
default. This is PySCF's `dft.gen_grid.Grids.radi_method` default and
the reason our UKS-LSDA total energies agree with PySCF references
to numerical precision (~1e-10 Ha at default grid, well within the
spec's 1e-6..1e-4 Ha targets). The earlier Mura-Knowles + Bragg-Slater
radii combination differed from PySCF's grid by O(mHa) per atomic
center and was a constant systematic bias against PySCF references.

Treutler-Ahlrichs uses atom-specific xi scaling
(`treutler_ahlrichs_xi(Z)` in grid.cpp, matches PySCF
`dft.radi._treutler_ahlrichs_xi` for Z=1..36).

**Rule**: when validating any grid-based DFT quantity against PySCF,
default to `radial_grid_treutler_ahlrichs`. Mura-Knowles is retained
in the public API for other DFT codes that use it (NWChem, ORCA), but
it is NOT the right default for PySCF parity work.

## [2026-05-19] Grid radii match PySCF's BRAGG_RADII for Z=1..36

`bragg_slater_radius(Z)` for Z=1..36 reproduces
`pyscf.dft.radi.BRAGG_RADII` exactly. Notably H uses 0.35 Å
(Slater's non-bonded radius, not the 0.25 Å covalent value), and
noble gases (He, Ne, Ar, Kr) use PySCF's Bondi-style values
(1.40, 1.50, 1.80, 1.90 Å). Z=37..86 still uses the original
Slater 1964 table. Matching PySCF radii is what makes our UKS
total energies comparable to PySCF references at the few-mHa level;
using the bare Slater 1964 values produces an extra ~3 mHa
systematic offset.

**Rule**: when validating any grid-based quantity against PySCF,
make sure the `BRAGG_RADII` table matches. Source of truth:
`pyscf/dft/radi.py`.

## [2026-04-04] molrs integration complete

molrs CXX FFI integration is done. Key decisions that diverged from original spec:
- No `MolFrame` type in C++ -- all frame data lives in `molrs::MolRec` on the Rust side
- `cpu::MolRec` uses PIMPL (not direct `rust::Box`), keeping `rust/cxx.h` out of public headers
- Bridge is precision-templated (`bridge.rs` generated from template via `build.rs`)
- `ZarrWriter` delegates to `molrs::MolRec::write_zarr` (no separate `zarrs_ffi` dependency)
- Mulliken analysis stays in C++ (xtensor-based), not migrated to molrs
- `ATV_HAS_ZARR` / `ATV_BUILD_ZARR` flags removed -- Zarr is always available via molrs

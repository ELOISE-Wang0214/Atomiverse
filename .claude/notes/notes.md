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

## [2026-04-04] molrs integration complete

molrs CXX FFI integration is done. Key decisions that diverged from original spec:
- No `MolFrame` type in C++ -- all frame data lives in `molrs::MolRec` on the Rust side
- `cpu::MolRec` uses PIMPL (not direct `rust::Box`), keeping `rust/cxx.h` out of public headers
- Bridge is precision-templated (`bridge.rs` generated from template via `build.rs`)
- `ZarrWriter` delegates to `molrs::MolRec::write_zarr` (no separate `zarrs_ffi` dependency)
- Mulliken analysis stays in C++ (xtensor-based), not migrated to molrs
- `ATV_HAS_ZARR` / `ATV_BUILD_ZARR` flags removed -- Zarr is always available via molrs

# Spec: `pw-pbe-gga-si` — plane-wave true PBE GGA (Si)

**Slug:** `pw-pbe-gga-si`  
**Domain:** `pw` (plane-wave DFT)  
**Queue:** CP2K capability loop — incomplete `[~]`  
**Branch:** `wf/pw-pbe-gga-si`  
**Status:** open

## Problem

`cpu::pw::internal::XCBuilder` already calls `evaluate_xc(XCType::PBE, …)` but
**never passes density gradients** (`sigma = 0` default). That is the
**LDA limit of PBE**, not true GGA. GPW already has spectral GGA
(`xc_grid_pbe`); plane-wave must match that physics.

Roadmap / queue: “PW XC still LDA-limit (`sigma=0`); need true GGA twin”.

## Goal

Ship **true GGA PBE** on the plane-wave FFT grid for spin-unpolarized KS:

- Spectral \(\nabla\rho\) from \(\rho(G)\) (or r-space + FFT)
- Pointwise PBE via shared `cpu::scf::evaluate_xc` (libxc-parity)
- Potential \(v_{\mathrm{xc}} = v_\rho - \nabla\cdot\big[(v_{\sigma_{\alpha\alpha}}+\tfrac12 v_{\sigma_{\alpha\beta}})\nabla\rho\big]\) (closed shell)
- Public selector on concrete `cpu::pw::PW` (member, **no** `*Config` bag)

## Out of scope

- Spin-polarized PW UKS (`pw-spin-polarized-o2`)
- Metallic smearing (`pw-smearing-al`)
- Analytic Nielsen–Martin stress (`pw-analytic-stress-si`)
- New UPF download at ctest time (ship/use in-tree NC-UPF only)
- Live CP2K/QE from ctest
- CUDA PW GGA (CUDA path may still be LDA-limit; not this slug)

## Public API (minimal)

```cpp
cpu::pw::PW pw;
pw.upf_path = "…/Si_lda_nc.upf";   // in-tree NC-UPF
pw.ecut_ry = …;
pw.k_mesh = …;
pw.xc = cpu::scf::XCType::PBE;     // true GGA when is_gga(xc)
pw.build(frame);
const auto& r = pw.run(frame);     // r.energy [Ry], energy on Frame via evaluator [eV]
auto F = pw_forces_analytic(frame, pw); // [eV/Å]
```

Rules:

- Parameters are **members** of `cpu::pw::PW` (iron law: no public `*Config`).
- Reuse `cpu::scf::XCType` / `is_gga` — do **not** invent a parallel public enum bag.
- `XCBuilder` stays internal (not a `Component`).

## Tasks

1. **Spec** (this file) + INDEX Live row.
2. **`PW::xc` member** (`XCType`, default `PBE`) wired through SCF → `XCBuilder::step`.
3. **True GGA path in `XCBuilder`:** when `is_gga(xc)`, compute spectral density
   gradients and the divergence term; when LDA, keep pointwise-only path.
4. **Unit tests** under `tests/cpu/planewave/` (no e2e):
   - Uniform \(\rho\): PBE GGA \(E_{\mathrm{xc}}\) matches `sigma=0` PBE (LDA limit).
   - Non-uniform \(\rho\): GGA \(E_{\mathrm{xc}}\) **differs** from `sigma=0` by a
     meaningful threshold (proves \(\sigma\neq 0\)).
5. **Offline pin + regression** `regressions/pw/pw-pbe-gga-si/`:
   - Si diamond primitive (2 atoms), in-tree `tests/data/pw/Si_lda_nc.upf`,
     modest `ecut_ry` + small k-mesh (≤ 60 s).
   - Observables: total energy finite/converged; ionic forces near zero by
     crystal symmetry; GGA vs LDA-limit total-energy gap documents active GGA.
6. **Perf gate:** fair wall vs CP2K/QE if binary present; else `perf_exempt`
   with reason recorded in `REFERENCE.md`.
7. Commit on `wf/pw-pbe-gga-si` only (never `dev`).

## Acceptance criteria

| # | Criterion | How verified |
| - | --------- | ------------ |
| A1 | `PW` exposes `xc` as a public member (`XCType`); no public `*Config` / `make_*` | Code review + compile |
| A2 | For `XCType::PBE`, XC uses non-zero \(\sigma\) on nonuniform density (not LDA-limit only) | Unit test A2 |
| A3 | Uniform density: GGA PBE \(E_{\mathrm{xc}}\) equals PBE@\(\sigma=0\) within \(10^{-10}\) Ha (relative to grid integral) | Unit test A3 |
| A4 | Regression `reg_pw-pbe-gga-si` builds with `ATV_BUILD_PW` + `ATV_BUILD_REGRESSIONS` | CMake + link |
| A5 | Si diamond SCF with `xc=PBE` converges (or reaches `max_iter` with finite energy); energy finite | Regression exit 0 |
| A6 | Ionic \(|F_i| < 0.5\) eV/Å on both atoms (ideal diamond sites; symmetry pin) | Regression |
| A7 | \(\|E_{\mathrm{PBE}}-E_{\mathrm{LSDA}}\|\) or \(\|E_{\mathrm{GGA}}-E_{\mathrm{PBE@}\sigma=0}\|\) documents GGA activity above stated floor in REFERENCE | Regression / REFERENCE |
| A8 | `REFERENCE.md` cites offline oracle (CP2K/QE version) **or** explicit unavailability + literature/symmetry basis for pins | File present |
| A9 | Dual gate: correctness (A5–A7) **and** (`within_10x` **or** `perf_exempt=true` with reason) | REFERENCE + report |
| A10 | No live third-party process spawned from ctest / regression binary | Grep / design |

## Offline pin strategy

1. Prefer CP2K Quickstep **or** Quantum ESPRESSO twin (same cell, NC-UPF or
   equivalent, PBE, comparable cutoff/k-mesh) run **offline** via
   `tools/cp2k_oracle/` notes; record version + inputs in
   `regressions/pw/pw-pbe-gga-si/REFERENCE.md`.
2. If no CP2K/QE binary on the worker node:  
   - **Force pin:** crystal symmetry \(\mathbf{F}=\mathbf{0}\) at ideal diamond
     sites (literature lattice \(a=5.43\) Å, FCC primitive 2-atom basis).  
   - **GGA activation pin:** internal dual run LSDA vs PBE **or** unit-test
     \(\sigma\neq 0\) with a floor on \(\Delta E_{\mathrm{xc}}\).  
   - **perf_exempt=true**, reason: no offline twin binary on host.  
   - Do **not** invent absolute total-energy digits from memory.

## Files (expected touch list)

- `.claude/specs/pw-pbe-gga-si/spec.md` (this)
- `.claude/specs/INDEX.md` (Live row)
- `include/atomiverse/cpu/planewave/pw_solver.h` (`xc` member)
- `include/atomiverse/cpu/planewave/internal/xc_builder.h`
- `src/cpu/planewave/internal/xc_builder.cpp`
- `src/cpu/planewave/pw_solver.cpp` (pass `xc`)
- `tests/cpu/planewave/test_pw_components.cpp` (or dedicated GGA unit tests)
- `regressions/pw/pw-pbe-gga-si/{reg_pw-pbe-gga-si.cpp,REFERENCE.md}`
- `regressions/CMakeLists.txt`

## Done when

All acceptance rows A1–A10 hold; `acceptance_met=true`; branch committed;
report JSON with `status=green` (or `blocked` only for speed after one optimize
pass / missing dep that cannot be fixed in-slug).

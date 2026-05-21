---
title: B3LYP and PBE0 global-hybrid single-point energies on CPU KS SCF
status: approved
created: 2026-05-21
---

# B3LYP and PBE0 global-hybrid single-point energies on CPU KS SCF

## Summary

Add the B3LYP and PBE0 global-hybrid density functionals to the Atomiverse CPU Kohn-Sham SCF path so that closed-shell single-point energies match PySCF to sub-microHartree accuracy. A global hybrid mixes a functional-specific fraction `alpha_hf` of exact Hartree-Fock exchange into the KS Fock matrix on top of a semi-local exchange-correlation part. This is delivered without any new driver, solver, or "Hybrid" class: hybrid is a property of the chosen functional, carried through the existing `SCFConfig`, the existing `cpu::scf::SCFSolver`, and its internal `HamiltonianBuilder`. Users select a hybrid through new `b3lyp_config` / `pbe0_config` convenience helpers, exactly parallel to `default_config` / `baseline_config`. Hybrid analytic gradients and range-separated hybrids are explicitly excluded (link 3 of this chain handles gradients).

## Domain basis

Global hybrids combine a semi-local XC functional with a fixed scalar fraction of exact HF exchange. The two functionals added here, with coefficients verified against PySCF + libxc source:

- **B3LYP** — libxc native id 402 (`XC_HYB_GGA_XC_B3LYP`). Semi-local components with ABSOLUTE coefficients: `0.08*Slater-X + 0.72*B88-X + 0.19*VWN_RPA-C + 0.81*LYP-C`, plus `alpha_hf = 0.20` exact HF exchange. `VWN_RPA` is libxc id 8 (the RPA parameterization), NOT VWN5 and NOT VWN3 — supplied by link 1 (`phase3-hybrid-dft-01-b88-lyp-functionals`).
- **PBE0** — libxc native id 406 (`XC_HYB_GGA_XC_PBEH`). Semi-local components: `0.75*PBE-X + 1.00*PBE-C`, plus `alpha_hf = 0.25`. PBE0 reuses the EXISTING native PBE engine; only the PBE exchange coefficient is scaled to 0.75 while PBE correlation stays at 1.0. No new semi-local functional is introduced for PBE0.
- **Hybrid restricted KS Fock**:
  `F = h_core + J[P] - alpha_hf*K[P] + V_xc^sl[rho]`
  where `V_xc^sl` is the SEMI-LOCAL XC potential only, evaluated with the absolute component coefficients above. `alpha_hf` appears in the Fock build at exactly one place — the `alpha_hf*K` term. The semi-local coefficients are already absolute; there is no `(1-alpha_hf)` rescaling of them. No double counting.
- **Hybrid energy decomposition**:
  `E = E_1e + E_J + alpha_hf*E_K + E_xc^sl + E_nuc`
  with `E_J = 0.5*Tr[P*J]`, `E_K = -0.5*Tr[P*K]` (exact HF exchange energy at `alpha_hf=1`), and `E_xc^sl` the semi-local XC energy from the component mix.
- **Convention pin (resolved here).** The existing Atomiverse K builders normalize K per the current HF Fock conventions: the RHF stored-ERI path produces `G = J - 0.5*K` (`build_coulomb_exchange_fused`, P = total double-occupied density), and the UHF path produces per-spin `K_alpha` / `K_beta` such that `F^sigma = Hcore + J - K^sigma` (`build_jk_uhf_direct`). For the closed-shell restricted hybrid path delivered here, the implementer MUST first read the K builder used and confirm its factor convention, then apply `alpha_hf` at ONE consistent location so the Fock term `-alpha_hf*K` and the energy term `alpha_hf*E_K` use the identically-scaled K. The recommended approach: obtain a "bare" K (the matrix that, used as `F = Hcore + J - K`, reproduces pure HF) and multiply that K by `alpha_hf` once. Do not split `alpha_hf` across builder and accumulator.
- B3LYP and PBE0 are GLOBAL hybrids: a single scalar `alpha_hf`, range-separation parameter `omega = 0`. Range-separated hybrids (CAM-B3LYP, wB97X) are out of scope.

References:
- B3LYP: Becke, J. Chem. Phys. 98, 5648 (1993), DOI 10.1063/1.464913; Stephens, Devlin, Chabalowski, Frisch, J. Phys. Chem. 98, 11623 (1994), DOI 10.1021/j100096a001.
- PBE0: Adamo, Barone, J. Chem. Phys. 110, 6158 (1999), DOI 10.1063/1.478522; Ernzerhof, Scuseria, J. Chem. Phys. 110, 5029 (1999), DOI 10.1063/1.478401.
- B88 exchange: Becke, Phys. Rev. A 38, 3098 (1988), DOI 10.1103/PhysRevA.38.3098.
- LYP correlation: Lee, Yang, Parr, Phys. Rev. B 37, 785 (1988), DOI 10.1103/PhysRevB.37.785.
- VWN: Vosko, Wilk, Nusair, Can. J. Phys. 58, 1200 (1980), DOI 10.1139/p80-159 (RPA parameterization).

## Design

**Functional selection.** Extend `enum class XCType` in `xc_functional.h` with two new entries: `XCType::B3LYP` and `XCType::PBE0`. This is the minimal reconciliation consistent with the real `SCFConfig`, which already carries `XCType xc_type` and is KS-only when `method == KS`. No `XCSpec` struct is introduced — a separate struct would fragment the existing `xc_type` field and force every call site to migrate; the stale monolithic draft's `XCSpec` proposal is rejected on that ground.

**Hybrid metadata.** Add a single free function in `xc_functional.h` / `.cpp`:
`double hybrid_exchange_fraction(XCType type)` returning `alpha_hf` — `0.20` for `B3LYP`, `0.25` for `PBE0`, and `0.0` for every non-hybrid type (`LSDA`, `PBE`). A `0.0` return is the signal "pure functional, no K build". This keeps the hybrid flag derivable from `xc_type` alone, with no new `SCFConfig` field. (An `is_gga(XCType)` predicate — true for `PBE`, `B3LYP`, `PBE0` — is also added here for link 3's reuse.)

**Semi-local evaluation.** `evaluate_xc` / `evaluate_xc_batch` gain two new cases. `B3LYP` evaluates the absolute-coefficient mix `0.08*Slater-X + 0.72*B88-X + 0.19*VWN_RPA-C + 0.81*LYP-C` using the link-1 native B88 / LYP / VWN_RPA component evaluators. `PBE0` evaluates `0.75*PBE-X + 1.00*PBE-C` by reusing the existing native PBE engine with the exchange contribution scaled by 0.75. The `XCOutput` struct is unchanged — hybrid affects only the semi-local energy density / potential values, not the data shape.

**Fock build.** In `HamiltonianBuilder::step`, the KS branch currently computes `F = Hcore + J + Vxc` with no K matrix. Extend it: when `hybrid_exchange_fraction(config.xc_type) > 0`, additionally build the exchange matrix K via the existing HF-path K builder (`build_jk_uhf_direct` reused with the restricted-density convention already used for the J build, or the stored-ERI fused builder, matching whatever path the J build uses for ERI availability), then form `F_alpha = Hcore + J - alpha_hf*K + Vxc_alpha`. `alpha_hf` multiplies K at this single site. The non-hybrid KS path is byte-for-byte unchanged (no K build, no extra cost).

**Energy decomposition.** The hybrid HF-exchange energy `alpha_hf*E_K` must be added to the total in `scf_solver.cpp`. The KS energy branch currently computes `E_total = E_one + E_J + fock.E_xc_grid + V_nuc` and derives `exchange_energy` by subtraction. `FockResult` gains a field `double E_exact_exchange = 0` carrying `alpha_hf*E_K` (zero for pure KS); the KS energy branch adds it to `E_total`, and `result_.exchange_energy` is set directly to that field rather than by subtraction. This keeps the decomposition exact and the `alpha_hf=0` regression bit-identical.

**User entry points.** Add `b3lyp_config(const IntegralPack&)` and `pbe0_config(const IntegralPack&)` inline helpers in `scf_solve.h`, cloned from `default_config` (Restricted, KS, MINAO guess, `GridLevel::MatchPySCF`, `tol = 1e-9`, DIIS on) but with `xc_type` set to `B3LYP` / `PBE0` respectively.

**Reference data.** The current per-molecule reference JSON carries a single `xc` key. Hybrid references are added under new keyed slots in each molecule JSON (`b3lyp` and `pbe0` energy + orbital-energy blocks) so the existing PBE references are untouched. `generate_references.py` gains a single-point hybrid-reference mode that, for a given molecule and `xc`, runs `pyscf.dft.RKS(mol, xc=...).kernel()` at 6-31G* on the starting geometry and records total energy plus HOMO/LUMO orbital energies, alongside the existing `pyscf_version` provenance field.

Lifecycle: no new objects with non-trivial lifetime. `hybrid_exchange_fraction` is a pure function; `XCType` enum values are compile-time constants; the K matrix lives for the duration of one `HamiltonianBuilder::step` call, same as the existing J matrix.

Performance note (non-binding): a hybrid SCF builds K every iteration, so a B3LYP / PBE0 single point is roughly 2x the cost of a pure-PBE run and comparable to PySCF's hybrid path. This is expected and is NOT a correctness gate; no wall-time acceptance criterion is imposed.

## Files to create or modify

- `include/atomiverse/cpu/scf/xc_functional.h` — add `XCType::B3LYP`, `XCType::PBE0`; declare `hybrid_exchange_fraction(XCType)` and `is_gga(XCType)`.
- `src/cpu/scf/xc_functional.cpp` — implement `hybrid_exchange_fraction` / `is_gga`; add B3LYP and PBE0 cases to `evaluate_xc` / `evaluate_xc_batch`.
- `include/atomiverse/cpu/scf/scf_solve.h` — add `b3lyp_config` and `pbe0_config` inline helpers.
- `include/atomiverse/cpu/scf/internal/hamiltonian_builder.h` — add `double E_exact_exchange` to `FockResult`.
- `src/cpu/scf/internal/hamiltonian_builder.cpp` — KS branch: build K and apply `-alpha_hf*K` when hybrid; populate `E_exact_exchange`.
- `src/cpu/scf/scf_solver.cpp` — KS energy branch: add `fock.E_exact_exchange` to `E_total`; set `result_.exchange_energy` from the field.
- `benchmarks/qc/reference/generate_references.py` — add hybrid single-point reference mode emitting energy + HOMO/LUMO into keyed JSON slots.
- `benchmarks/qc/reference/h2o.json` — add `b3lyp` / `pbe0` reference blocks.
- `benchmarks/qc/reference/ch4.json` — add `b3lyp` / `pbe0` reference blocks.
- `benchmarks/qc/reference/nh3.json` — add `b3lyp` / `pbe0` reference blocks.
- `benchmarks/qc/reference/benzene.json` — add `b3lyp` / `pbe0` reference blocks.
- `benchmarks/qc/reference/glycine.json` — add `b3lyp` / `pbe0` reference blocks.
- `tests/cpu/scf/test_hybrid_singlepoint.cpp` (new) — B3LYP/PBE0 energy + HOMO/LUMO parity vs PySCF, `alpha_hf=0` regression.
- `tests/CMakeLists.txt` — register `test_hybrid_singlepoint`.

## Tasks

- [ ] Write failing tests for hybrid single-point energy and HOMO/LUMO parity (tests/cpu/scf/test_hybrid_singlepoint.cpp)
- [ ] Add XCType::B3LYP / XCType::PBE0, hybrid_exchange_fraction and is_gga in include/atomiverse/cpu/scf/xc_functional.h
- [ ] Implement hybrid_exchange_fraction / is_gga and B3LYP/PBE0 semi-local cases in src/cpu/scf/xc_functional.cpp
- [ ] Add E_exact_exchange to FockResult and build alpha_hf*K in the KS branch of src/cpu/scf/internal/hamiltonian_builder.cpp
- [ ] Add fock.E_exact_exchange into the KS total energy in src/cpu/scf/scf_solver.cpp
- [ ] Add b3lyp_config and pbe0_config helpers in include/atomiverse/cpu/scf/scf_solve.h
- [ ] Add hybrid single-point reference mode to benchmarks/qc/reference/generate_references.py and regenerate b3lyp/pbe0 blocks in the h2o/ch4/nh3/benzene/glycine JSONs
- [ ] Register test_hybrid_singlepoint in tests/CMakeLists.txt
- [ ] Add docstrings per doxygen with units for hybrid_exchange_fraction, is_gga, the new XCType cases, FockResult::E_exact_exchange, and the config helpers
- [ ] Run full check + test suite

## Testing strategy

- Happy path: B3LYP single-point total energy on H2O / CH4 / NH3 / benzene / glycine at 6-31G* matches the checked-in PySCF `RKS(mol, xc='b3lyp')` reference to <= 1e-7 Ha.
- Happy path: PBE0 single-point total energy on the same molecule set matches the checked-in PySCF `xc='pbe0'` reference to <= 1e-7 Ha.
- Domain validation: HOMO and LUMO orbital energies for B3LYP and PBE0 on the molecule set match the PySCF references to <= 1e-5 eV (correct mixing of exact exchange into the eigenspectrum).
- Edge case / regression: running the hybrid Fock+energy code path with `alpha_hf` forced to 0 reproduces the pure-PBE result bit-for-bit (the hybrid path must collapse exactly onto the existing pure-DFT path; no `(1-alpha_hf)` rescaling artifacts).
- Domain validation: energy-decomposition self-consistency — `E_total == E_1e + E_J + alpha_hf*E_K + E_xc^sl + E_nuc` to machine precision, confirming `alpha_hf` is applied at a single consistent site in both Fock and energy.
- Config check: `b3lyp_config` / `pbe0_config` produce a `SCFConfig` with `xc_type` equal to `B3LYP` / `PBE0` and `hybrid_exchange_fraction` returning 0.20 / 0.25 for those types.
- Full `ctest` suite stays green (existing pure-HF/PBE/LSDA tests unaffected).

## Out of scope

- Hybrid analytic gradient and KS hybrid forces — link 3 of this chain (`phase3-hybrid-dft-03-*`).
- Geometry optimization with hybrid functionals — link 3.
- Range-separated / long-range-corrected hybrids (CAM-B3LYP, wB97X, `omega != 0`).
- Unrestricted (open-shell) hybrid SCF — the closed-shell restricted path is the deliverable; UKS hybrid is deferred (no open-shell hybrid reference is added here).
- Wall-time / performance acceptance — hybrid is K-heavy (~2x pure PBE); performance is a non-binding note in the Design section, not a correctness gate.
- Double-hybrid functionals (B2PLYP and similar) requiring an MP2 correlation term.

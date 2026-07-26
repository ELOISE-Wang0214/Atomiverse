---
slug: esmm-qmmm-01-embed-forces
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: multiscale / qmmm (SCF embedding physics)
queue: esmm-campaign
deferred: false
depends_on: []
unblocks: [qmmm-mm-site-forces]
grilled: true
---

# Spec: `esmm-qmmm-01-embed-forces` — analytic embedding gradient + forces on the MM charges

**Wave:** esmm Wave 2 (QM/MM redesign), stage 1 of 5.
**Status:** approved — grilled via plan-mode design review + operator Q&A, 2026-07-26.

## Summary

The SCF point-charge embedding (`MeanFieldMethod::point_charges`) computes energy and
QM-nuclei gradients, but (a) the electron–MM leg is central finite difference —
`O(6N × n_charges)` full `core_matrix` rebuilds per force call
(`src/cpu/scf/internal/gradient.cpp:663-722`; FD block `:690-703`) — and (b) **no force
is ever computed on the MM sites**: Newton's third law is violated and any QM/MM MD
leaks momentum every step. This spec replaces the FD leg with the analytic
nuclear-attraction derivative already present in the engine, adds the site-force legs
(electron + nucleus), routes them out through `GradientResult`, and unifies the three
redundant site structs (`MeanFieldMethod<D>::PointCharge` — a **template-nested type**,
so `RHF::PointCharge` ≠ `UKS::PointCharge` and no erased interface can traffic in it —
plus `internal::PointChargeSite`) into one public `cpu::EmbeddingSite`.
**No new integral math is written** — both primitives exist as file-statics in
`src/cpu/scf/internal/integral_derivatives.cpp`.

## Domain basis

Embedding energy (Hartree, internal a.u.):

E_embed = Σ_C q_C Tr[P · v_C] + Σ_A Σ_C Z_A q_C / |R_A − R_C|  (+ Σ_{C<C'} q_C q_C' / |R_C − R_C'| when self-energy is enabled),

with v_C the nuclear-attraction operator at the site. Gradients:

- **QM-center leg** (∂/∂R_A of Tr[P v_C]): bra/ket-center derivative of the
  nuclear-attraction integral — the standard raising-formula path already used by
  `IntegralDerivPack::dV` (`nuclear_deriv_prim_center_general`,
  `integral_derivatives.cpp:147`).
- **Site leg** (∂/∂R_C): translational invariance gives
  `dV_dC[d] = −(dV_dA[d] + dV_dB[d])` (`integral_derivatives.cpp:306-308`) — that
  **is** the electric-field integral of the QM density at R_C; plus the analytic
  nucleus–site Coulomb term `+Σ_A Z_A q_C (R_A−R_C)/|R_A−R_C|³` (sign flip of the
  existing code at `gradient.cpp:704-718`).
- With site self-energy **off**, Σ_QM F + Σ_MM F = 0 to SCF convergence.

References: Warshel & Levitt, J. Mol. Biol. 103, 227 (1976); Senn & Thiel, Angew.
Chem. Int. Ed. 48, 1198 (2009); PySCF `examples/qmmm/30-force_on_mm_particles.py`
(the MM-side force recipe used for the pin). Units: internal Ha/Bohr; public eV/Å, e, Å.

## Design

| Piece | Location | Change |
| ----- | -------- | ------ |
| `cpu::EmbeddingSite{x,y,z,q}` [Å, e], `cpu::EmbeddingForce{fx,fy,fz}` [eV/Å] | **new** `include/atomiverse/cpu/component/qm_engine.h` (structs only; the protocol classes land in `esmm-qmmm-02`) | replaces `MeanFieldMethod<D>::PointCharge` + `internal::PointChargeSite` |
| derivative primitives | `internal/integral_derivatives.{h,cpp}` | promote `nuclear_deriv_prim_center_general` (`:147`) and `nuclear_deriv_prim_nuc_general` (`:221`) from file-static to internal-header declarations (unit-testable) |
| `PointChargeEmbedding` | `internal/point_charges.{h,cpp}` | ctor takes `std::vector<cpu::EmbeddingSite>`; **three new accumulators**: `accumulate_qm_center_gradient(grad, shells, P_cart, shell_to_atom)`, `accumulate_site_gradient(site_grad, shells, P_cart, atoms)`, `accumulate_self_gradient(site_grad)` — all [Hartree/Bohr], contracted on the fly (no [n_atoms,3,n,n] tensor) |
| `GradientResult` | `scf/gradient_result.h` | new `xt::xtensor<double,2> embedding_forces; ///< [n_charges,3] eV/Å` (empty when no charges — honest absence, not a fallback) |
| gradient assembly | `internal/gradient.cpp:663-722` | FD block deleted; analytic QM-center leg + site-gradient assembly into `embedding_forces`; self-gradient only when enabled |
| `MeanFieldMethod<D>` | `scf/mean_field_method.{h,cpp}` | `point_charges` becomes `std::vector<EmbeddingSite>`; new public member `bool embed_self_energy = true` (gates `emb.self_energy()` at `mean_field_method.cpp:~585` **and** the self-gradient); new private `std::vector<EmbeddingForce> embed_site_forces_` republished from `last_grad_.embedding_forces` + accessor `const std::vector<EmbeddingForce>& embedding_site_forces() const` |

`embed_self_energy = true` default **preserves** the standalone-embedding energy
convention pinned by `regressions/scf/h2o_embedding_field` (E = −74.98418931920534 Ha
includes MM–MM). Additive QM/MM couplings (Wave-2 stage 3) set it `false` so MM–MM
Coulomb is owned in exactly one place — the MM force field.

**EmbeddingSite rename touch list (8 sites, no aliases, no back-compat):**
`mean_field_method.h` (delete nested struct `:308-316`), `mean_field_method.cpp`
(drop the `PointChargeSite` conversion loop `:573-589`), `internal/point_charges.h`,
`internal/point_charges.cpp`, `internal/integrals.h` (`:92` `embed_charges` element
type), `internal/gradient.cpp`, `src/cpu/multiscale/qmmm_electrostatic.cpp`
(`rhf.point_charges` fill site), tests + regressions
(`tests/cpu/scf/test_embedding.cpp`, `tests/cpu/semiempirical/test_am1_scaffold.cpp`,
`regressions/scf/embed_h2o_grad/`, `regressions/scf/h2o_embedding_field/`).

## Files to create or modify

- `include/atomiverse/cpu/component/qm_engine.h` (new — structs only this spec)
- `include/atomiverse/cpu/scf/internal/integral_derivatives.h` + `src/.../integral_derivatives.cpp`
- `include/atomiverse/cpu/scf/internal/point_charges.h` + `src/.../point_charges.cpp`
- `src/cpu/scf/internal/gradient.cpp`, `include/atomiverse/cpu/scf/gradient_result.h`
- `include/atomiverse/cpu/scf/mean_field_method.h`, `src/cpu/scf/mean_field_method.cpp`
- `include/atomiverse/cpu/scf/internal/integrals.h`
- rename fallout: `src/cpu/multiscale/qmmm_electrostatic.cpp`, 2 tests, 2 regressions
- `tests/cpu/scf/test_embedding.cpp` (extend), `tests/CMakeLists.txt` if new files
- `regressions/scf/embed_h2o_mm_forces/{reg_embed_h2o_mm_forces.cpp,REFERENCE.md}` (new) + `regressions/CMakeLists.txt`

## Tasks

- [ ] RED: unit test — analytic QM-center leg vs central FD of `core_matrix` (h = 1e-4 Bohr) ≤ 1e-6 Ha/Bohr on a 2-charge H₂O toy (`tests/cpu/scf/test_embedding.cpp`)
- [ ] RED: unit test — analytic site gradient vs central FD of the total energy w.r.t. one charge position ≤ 1e-6 Ha/Bohr; `Σ_QM F + Σ_MM F = 0` ≤ 1e-9 eV/Å with `embed_self_energy = false`; `embedding_forces` empty when `point_charges` empty
- [ ] Promote the two derivative primitives; add the three `PointChargeEmbedding` accumulators
- [ ] Introduce `cpu::EmbeddingSite`/`EmbeddingForce`; rename across the 8 sites; add `embed_self_energy` + `embedding_site_forces()`
- [ ] Rewrite `gradient.cpp:663-722` (analytic legs; fill `GradientResult::embedding_forces`; self-gradient gated)
- [ ] Verify ECP (`ecp != ""`) + `point_charges` reachability: if reachable, typed `std::runtime_error` on the analytic path (no silent FD fallback) + guard test
- [ ] Offline PySCF 2.14 pin (tools/pyscf_oracle, x86_64 venv; `examples/qmmm/30-force_on_mm_particles.py` recipe) on the `embed_h2o_grad` geometry; verify ATV reproduces BEFORE pinning; write `REFERENCE.md`
- [ ] Regression `regressions/scf/embed_h2o_mm_forces` (public API only; ≤ 60 s; no live PySCF)
- [ ] GREEN: full unit suite + existing `embed_h2o_grad` / `h2o_embedding_field` regressions unchanged; Doxygen (units) on every new public symbol

## Testing strategy

Unit: parity FD↔analytic both legs; third-law sum; empty-absence; self-energy gate
(energy with flag on/off differs by exactly the closed-form MM–MM Coulomb).
Regression: `embed_h2o_mm_forces` — MM-side gradient element-wise ≤ 1e-5 Ha/Bohr vs
the offline PySCF pin, plus Σ F = 0 with self-energy off. Perf: FD block removal is a
strict win (`O(6N·n_charges)` sweeps → one contraction); `perf_exempt` with reason if
no twin wall exists.

## Out of scope

- The `QmEngine`/`EmbeddingEngine` protocol classes (`esmm-qmmm-02`)
- Coupling rewrite / mask column / cutoffs (`esmm-qmmm-03`)
- GPW/PW/AM1 embedding hooks, GEEP, PME-QM/MM, periodic embedding sums
- CUDA; public Config bags / factories

## Notes for implementers

- **HARD GATE:** before touching `src/cpu/scf/*` (`mean_field_method.*`,
  `internal/gradient.cpp`, `internal/integral_derivatives.*`, `internal/point_charges.*`)
  or any CMakeLists, verify `git status` shows them clean. As of 2026-07-26 the
  operator's libcint-removal/C++20 work (~760 files) is uncommitted and includes them.
  If dirty: STOP and report. (`cpu/multiscale/` itself is clean.)
- The republish into `embed_site_forces_` happens in the `drive_scf_core_` epilogue
  (`mean_field_method.cpp:958-985` region) — keep it out of the stage components.
- Keep `h2o_embedding_field`'s include-self-energy pin standing: default `true`.

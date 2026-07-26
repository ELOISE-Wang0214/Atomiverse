---
slug: esmm-qmmm-03-electrostatic-rewrite
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: multiscale / qmmm (coupling rewrite)
queue: esmm-campaign
deferred: false
depends_on: [esmm-qmmm-02-qm-engine-protocol, esmm-contract-04-forces-requested, esmm-contract-05-accumulate-switch]
unblocks: [qmmm-engine-slot]
grilled: true
---

# Spec: `esmm-qmmm-03-electrostatic-rewrite` — mask-column QM region + engine slot + real MM coupling

**Wave:** esmm Wave 2 (QM/MM redesign), stage 3 of 5 — the connection layer itself.
**Status:** approved — grilled via plan-mode design review + operator Q&A, 2026-07-26.

## Summary

Rewrite `QmmmElectrostatic` from the shipped toy (hard-coded `RHF{STO3G}` rebuilt cold
every launch; uniform `mm_charge` scalar; single global LJ; cutoff-free embedding sum
with raw positions; per-launch sub-Frame allocation; `type` column silently read as Z)
into the real coupling: **any embedding-capable engine** via the `QmEngine` slot
(persistent across steps ⇒ ASPC warm start works inside QM/MM), **QM region declared
by a per-atom mask column** (operator decision, 2026-07-26), per-atom MM charges from
`kAtom_charge`, an embedding cutoff with minimum image, per-type QM–MM LJ, third-law
site-force scatter, and a persistent private QM frame. The coupling owns **only** the
QM-involving cross terms; MM–MM LJ/Coulomb/bonded/PME belong to the user's ordinary
MM styles — that is what makes QM/MM compose with the force-field library.

## Operator decision — the QM mask column (binding)

**QM region = per-atom mask column on the Frame `atoms` block** (operator Q&A
2026-07-26; overrides the earlier ctor-tag-list design):

- Canonical key **`kAtomQmZ`** = `ColKey<int>{"atoms", "qm_z", kDynamicFieldId, …}`,
  published in new header `cpu/multiscale/qm_mask.h`.
- Value **0 = MM atom; > 0 = QM atom whose atomic number Z = the value** — one column
  answers membership AND element, killing the type-column-as-Z hazard (`type` keeps
  its force-field meaning and keys the LJ table).
- v1 = **runtime in-memory column** (`Block::insert`, `wire_id = kDynamicFieldId`;
  `gather.cpp:396` precedent; never crosses the molrs boundary). The user creates and
  fills it before `build()`; the coupling resolves it at `build()` and **throws with a
  how-to message** when the column is missing or all-zero.
- Runtime columns do not travel in Transaction patches or halo packs ⇒ QM membership
  is **per-run static**: `QmSubsystem::refresh` throws on QM-count change (MC
  insert/delete of a QM atom); `build()` throws when `mpi_size > 1` (v1 serial).
- `LinkBond{qm_tag, mm_tag}` stays a member input (topology cuts, few in number).
- molrs schema round-trip of the column is a named follow-on seam, not v1.

## Domain basis

Additive electrostatic embedding (Senn & Thiel 2009 §2; Warshel & Levitt 1976):
E = E_QM[ρ; sites] + E_MM(MM styles) + E_QM–MM^vdW, with the MM sites inside the QM
Hamiltonian (esmm-qmmm-01 physics) and `set_site_self_energy(false)` so MM–MM Coulomb
is owned exactly once — by the force field. Exclusion story requires **zero new
machinery**: the coupling validates at `build()` that every QM atom carries
`kAtom_charge == 0` (throw listing offending tags) — then any Coulomb style yields
exactly 0 for every QM-involving pair; QM types simply have no diagonal in the user's
`PairLJCut` table ("a type with no diagonal has no LJ site", `pair_lj_cut.h:63-66`),
and the coupling supplies QM–MM LJ itself. Units: eV, Å, e, eV/Å.

## Design

| Piece | Location | Role |
| ----- | -------- | ---- |
| `kAtomQmZ` + Doxygen how-to | new `cpu/multiscale/qm_mask.h` | canonical mask key |
| `internal::QmSubsystem` | new `cpu/multiscale/internal/qm_subsystem.{h,cpp}` | mask resolve (build-time), persistent private QM `Frame` (allocated once in `build`, refreshed per launch), link caps, **the engine slot**, site collection, force scatter. Internal member class, not a Component and not a base (`internal::Pcm`/`SOSCFStepper` precedent) |
| `internal::QmMmVdw` | new `cpu/multiscale/internal/qm_mm_vdw.{h,cpp}` | per-(type,type) 12-6 LJ, Lorentz–Berthelot mixing from diagonals, shifted at cutoff, minimum image, QM–MM pairs only, O(N_qm × N_local) |
| `QmmmElectrostatic : Potential` | rewrite `cpu/multiscale/qmmm_electrostatic.{h,cpp}` | composition of the two internals + embedding wiring |

### Public API (binding shape)

```cpp
/// ctor: (cross_lj type-pair rows {type_i, type_j, eps [eV], sigma [Å]},
///        lj_cutoff [Å], embed_cutoff [Å]) — NO qm-atom argument: the QM
/// region comes from the kAtomQmZ mask column. Reuse PairLJCut's coefficient
/// row shape if public; else define the row struct in qm_mm_vdw.h.
QmmmElectrostatic(std::vector<LJTypeCoeff> cross_lj, double lj_cutoff,
                  double embed_cutoff);
std::vector<LinkBond>& links();          // covalent cuts (empty = none)
double& link_scale();                    // g ∈ (0,1), default 0.709
bool& reset_history_each_launch();       // MC posture: engine cold start
template <typename E> requires(!std::is_lvalue_reference_v<E>)
QmmmElectrostatic& set_engine(E&& e);    // owned, heap-stable address
template <typename E> QmmmElectrostatic& set_engine_ref(E& e);  // borrowed
template <typename E> QmmmElectrostatic& set_engine(E&) = delete;
double last_qm_energy() const;           // [eV]
double last_vdw_energy() const;          // [eV]
int last_n_links() const;
std::array<double,3> last_force_residual() const;  // Σ_all F [eV/Å] witness
```

Per-launch flow: `refresh` (re-resolve mask→rows, copy x into the private frame with
`z` from the mask value, place link H, zero private energy/forces) → `collect_sites`
(MM atoms with |ΔR_min-image| ≤ embed_cutoff of any QM atom; charges from
`kAtom_charge`; link MM partners omitted — Z1 charge deletion) → `evaluate` (engine on
the private frame; honors `ctx.forces_requested` from esmm-contract-04) → `scatter`
(`+=` QM forces with link chain rule, `+=` site forces onto MM rows, `+=` energy into
`kEnergy` **and the QM share into `kEnergyElectronic`** per the Wave-1 contract).
`build()` throws when: no engine; engine not embedding-capable
(`as_embedding()==nullptr`, message names `qm_engine_name()`); mask column
missing/all-zero; a QM atom has nonzero `kAtom_charge`; a link tag is absent;
`mpi_size > 1`. **Not** a `NeighborConsumer` in v1 (declaring it would make a 5-atom
`Evaluate` illegal — drivers throw on a consumer with no neighbor structure; the cost
left after MM–MM removal is O(N_qm×N_mm), noise next to one SCF; seam recorded).

### Usage (end to end — the 自由组合 witness)

```cpp
cpu::scf::RKS ks{cpu::scf::XCType::PBE0};
ks.basis = cpu::scf::BasisSet::Def2SVP;  ks.predictor(3, 2);
frame.blocks.at("atoms").insert("qm_z", /* int column spec */);
auto qm = frame.span(cpu::kAtomQmZ);
qm[0] = 8; qm[1] = 1; qm[2] = 1;              // QM water; all others stay 0
cpu::QmmmElectrostatic qmmm{{{2,2,0.0065,3.15}}, 10.0, 12.0};
qmmm.set_engine(std::move(ks));
md.add_potential(cpu::PairLJCut{...});         // MM–MM only (no QM diagonals)
md.add_potential(cpu::CoulLong{10.0});         // QM q=0 ⇒ no double count
md.add_potential(std::move(qmmm));
```

## Files to create or modify

- new: `cpu/multiscale/qm_mask.h`; `cpu/multiscale/internal/qm_subsystem.{h,cpp}`;
  `cpu/multiscale/internal/qm_mm_vdw.{h,cpp}`
- rewrite: `cpu/multiscale/qmmm_electrostatic.{h,cpp}`
- CMake re-gating: root `CMakeLists.txt` (multiscale sources out of the
  `ATV_BUILD_SEMIEMPIRICAL` block — Electrostatic/Polarizable need `ATV_BUILD_HF`
  only; Mechanical stays wherever its default-free engine story lands in stage 4);
  `tests/CMakeLists.txt`; `regressions/CMakeLists.txt` (per-case gates)
- tests: new `tests/cpu/multiscale/test_qm_subsystem.cpp`,
  `tests/cpu/multiscale/test_qm_mm_vdw.cpp`; rewrite
  `tests/cpu/multiscale/test_qmmm_electrostatic.cpp`
- regression: new `regressions/multiscale/qmmm_es_h2o_pyscf/` (supersedes and deletes
  `qmmm_es_h2` — superseding intent recorded here)

## Tasks

- [ ] RED: `test_qm_subsystem.cpp` — mask resolve (values→rows+Z); missing column /
  all-zero column throw with how-to text; private frame allocated once (stable address
  across launches); QM-count change throws; nonzero QM `kAtom_charge` throws listing
  tags; link-H placement + (1−g)/g force split; site collection honors embed_cutoff +
  minimum image + per-atom charges; link MM partner omitted; no engine ⇒ build throws
- [ ] RED: `test_qm_mm_vdw.cpp` — LB mixing; no-diagonal ⇒ ε = 0; shift continuity at
  cutoff; QM–QM and MM–MM pairs untouched
- [ ] RED: `test_qmmm_electrostatic.cpp` — non-embeddable engine (AM1) throws naming
  the engine; RHF↔RKS engine swap both run; `forces_requested=false` skips force
  assembly; `+=` semantics (pre-seeded `kEnergy` preserved); `kEnergyElectronic`
  carries the QM share
- [ ] Implement `qm_mask.h`, `QmSubsystem`, `QmMmVdw`, `QmmmElectrostatic`
- [ ] CMake re-gating (build multiscale with `ATV_BUILD_HF` alone)
- [ ] Offline PySCF pin reuse + `REFERENCE.md`; regression
  `qmmm_es_h2o_pyscf` (all-zero cross-LJ table so the coupling reproduces the bare
  PySCF `qmmm.mm_charge` numbers): energy, QM forces, MM-site forces, |ΣF| ≈ 0
- [ ] Delete `regressions/multiscale/qmmm_es_h2/`; GREEN: full suite + Doxygen

## Testing strategy

Unit: the three RED batches above (mirror paths, seconds). Regression: one —
`qmmm_es_h2o_pyscf`, proving the coupling adds no error over the esmm-qmmm-01 SCF pin
(energy + both force sides ≤ 1e-5 Ha/Bohr-equivalent; ΣF ≤ 1e-9 eV/Å). Engine-swap
has no independent literature value ⇒ unit-test territory, per the
one-regression-per-feature rule.

## Out of scope

- Mechanical/Polarizable port (`esmm-qmmm-04`); MD conservation closure (`esmm-qmmm-05`)
- NeighborConsumer path; PME-QM/MM; switching functions at embed_cutoff; GEEP;
  ONIOM; adaptive QM/MM; MPI replicated-frame mode; molrs schema column; CUDA

## Notes for implementers

- **HARD GATE:** before touching `src/cpu/scf/*` or any CMakeLists, verify
  `git status` shows them clean. As of 2026-07-26 the operator's libcint-removal/
  C++20 work (~760 files) is uncommitted and includes them. If dirty: STOP and
  report. (`cpu/multiscale/` itself is clean.)
- Depends on `esmm-contract-04` (`forces_requested` plumbing) and `esmm-contract-05`
  (`kEnergyElectronic` key + accumulate contract) — verify both landed.
- Engine storage must be heap-stable (`unique_ptr` for owned) so warm-start state
  survives the coupling being moved into a Driver's `owned_components_`.
- Hoist every `span()` at the top of `launch()` (resolution budget); use
  `owned_count(ctx)`-consistent loops even though v1 is serial.
- `reset_history_each_launch` calls `engine->reset_history()` before evaluate —
  document (MC: rejected geometry must not poison the next trial's guess).
- Sharp embed_cutoff ⇒ energy discontinuity when an MM atom crosses it: document; the
  regression uses embed_cutoff ≥ half the box diagonal for exactness.

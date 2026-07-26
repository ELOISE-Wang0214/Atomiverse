---
slug: esmm-qmmm-02-qm-engine-protocol
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: multiscale / qmmm (component protocol)
queue: esmm-campaign
deferred: false
depends_on: [esmm-qmmm-01-embed-forces]
unblocks: [qm-engine-protocol]
grilled: true
---

# Spec: `esmm-qmmm-02-qm-engine-protocol` — the erased QM-engine slot

**Wave:** esmm Wave 2 (QM/MM redesign), stage 2 of 5.
**Status:** approved — grilled via plan-mode design review + operator Q&A, 2026-07-26.

## Summary

Define the protocol by which a multiscale coupling drives **any** concrete
electronic-structure method — the connection layer the operator asked for. Engines are
three different shapes today: `MeanFieldMethod<D>` (a `Driver<scf::Phase>` + `HF` tag,
×7 concrete methods), `AM1` (a plain method class), `DFTB2` (a Component). This spec
adds `cpu::QmEngine` (lifecycle = the Component lifecycle:
`prepare(BuildContext)/evaluate(LaunchContext)`) and the opt-in capability
`cpu::EmbeddingEngine` (probed via `as_embedding()`, default `nullptr` — the
`as_neighbor_consumer()` precedent), then implements them on `MeanFieldMethod` (once,
in the CRTP base — all 7 Gaussian methods), `AM1`, and `DFTB2`. **Pure interface work;
no behavior change anywhere.** The structs (`EmbeddingSite`/`EmbeddingForce`) and the
site-force physics landed in `esmm-qmmm-01`.

## Domain basis

Not a physics spec. The one physical contract: `evaluate` on the coupling's private QM
frame writes total energy [eV] to `kEnergy` and (when `qm_ctx.forces_requested`)
forces [eV/Å] to `kAtom_f*`; engines that can compute per-atom charges write them [e]
into the QM frame's `kAtom_charge` (Frame as sole currency — the same channel a future
GPW engine uses; `AM1Evaluator` already does this, `am1_evaluator.cpp:89`; the
Gaussian path reuses `mulliken_charge` from `observe/mulliken.h`).

## Design

### `include/atomiverse/cpu/component/qm_engine.h` (extend the structs-only header)

```cpp
/// Two protocols, one header (the neighbor_protocol.h precedent). Neither is a
/// Component: implementors are Drivers (MeanFieldMethod<D>), plain method
/// classes (AM1), or Components (DFTB2). Contract on the frame handed to
/// prepare/evaluate: it is the coupling's PRIVATE QM frame (QM atoms + link
/// caps only), pre-zeroed on energy and forces — an engine may overwrite or
/// accumulate, both are correct. The engine MUST NOT consult
/// ctx.runtime.pair_list or halo state (whole-system replicated, no neighbors).
class QmEngine {
  public:
    virtual ~QmEngine() = default;
    /// @brief Engine label for diagnostics / exception messages.
    virtual std::string qm_engine_name() const = 0;
    /// @brief One-time setup against the QM reference geometry.
    virtual void prepare(const BuildContext& qm_ctx) = 0;
    /// @brief One QM evaluation on the private QM frame (energy [eV] to
    ///        kEnergy; forces [eV/Å] iff qm_ctx.forces_requested; per-atom
    ///        charges [e] to kAtom_charge when the engine can).
    virtual void evaluate(const LaunchContext& qm_ctx) = 0;
    /// @brief Electrostatic-embedding capability, or nullptr.
    virtual EmbeddingEngine* as_embedding() { return nullptr; }
    /// @brief Drop warm-start / density-extrapolation history (MC posture).
    virtual void reset_history() {}
};

class EmbeddingEngine {
  public:
    virtual ~EmbeddingEngine() = default;
    /// @brief Replace the external charge set (positions [Å], charges [e]).
    virtual void set_embedding_sites(std::vector<EmbeddingSite> sites) = 0;
    /// @brief Site forces [eV/Å] of the last forces-enabled evaluate, ordered
    ///        as given; with QM forces these satisfy Σ F = 0 to SCF tol.
    virtual const std::vector<EmbeddingForce>& embedding_site_forces() const = 0;
    /// @brief false ⇒ MM–MM Coulomb among the sites is EXCLUDED (additive
    ///        QM/MM; MM–MM belongs to the force field). Default true.
    virtual void set_site_self_energy(bool enabled) = 0;
};
```

Rejected alternative (recorded): `evaluate(Frame&, bool want_forces)` — a
Component-shaped engine (`DFTB2::launch`) cannot be reached without a `Context&`, so
every such engine would need an adapter holding a fabricated context. Expressing the
protocol in Build/LaunchContext makes Component-shaped engines one-line forwarders and
"energy-only" rides the existing `LaunchContext::forces_requested` — no new concept.

### Implementors

| Engine | Bases added | Mapping |
| ------ | ----------- | ------- |
| `MeanFieldMethod<D>` (×7) | `QmEngine, EmbeddingEngine` | `qm_engine_name()` → method label; `prepare` → no-op (integrals need geometry); `evaluate` → `DriveOpts{write_energy=true, forces=qm_ctx.forces_requested, unwrap=true}` → `drive_scf_` + Mulliken writeback to `kAtom_charge` (reads the **parameter**, never the `compute_force` member — that stays the nested-driver knob, `mean_field_method.cpp:288`); `as_embedding()` → `this`; `reset_history()` → clear `aspc_alpha_/beta_`, `has_last_density_`; `set_embedding_sites` → `point_charges = std::move(sites)`; `embedding_site_forces()` → `embed_site_forces_` (from esmm-qmmm-01); `set_site_self_energy` → `embed_self_energy` |
| `AM1` | `QmEngine` | `prepare` → `build(frame)`; `evaluate` → `run` + (iff forces_requested) `compute_forces` + energy/forces/charges writeback onto the QM frame; `as_embedding()` → nullptr (NDDO has no external-charge path yet); `reset_history()` → drop persisted density |
| `DFTB2` | `QmEngine` (already `Semiempirical : Component`) | `prepare(c)` → `build(c)`; `evaluate(c)` → `launch(c)` — one-line forwards |

`AM1Evaluator::launch` then **delegates** to the solver's `evaluate` (dedup — the
energy/forces/charges writeback exists once).

## Files to create or modify

- `include/atomiverse/cpu/component/qm_engine.h` (extend)
- `include/atomiverse/cpu/scf/mean_field_method.h`, `src/cpu/scf/mean_field_method.cpp`
- `include/atomiverse/cpu/semiempirical/am1_solver.h`, `src/cpu/semiempirical/am1_solver.cpp`
- `src/cpu/semiempirical/am1_evaluator.cpp` (delegate)
- `include/atomiverse/cpu/semiempirical/dftb2.h`, `src/cpu/semiempirical/dftb2.cpp`
- `tests/cpu/component/test_qm_engine.cpp` (new), `tests/CMakeLists.txt`

## Tasks

- [ ] RED: `tests/cpu/component/test_qm_engine.cpp` — `as_embedding()` non-null for RHF/UHF/RKS, null for AM1/DFTB2; erased `evaluate` == direct `run` energy **bit-for-bit** on H₂/STO-3G; `forces_requested=false` leaves the QM frame's forces untouched; `reset_history` clears the warm start (second evaluate after reset reconverges from scratch — iteration count equals the cold count)
- [ ] Extend `qm_engine.h` with the two protocol classes (Doxygen, units, the private-frame contract)
- [ ] Implement on `MeanFieldMethod` (CRTP base, once) incl. Mulliken `kAtom_charge` writeback
- [ ] Implement on `AM1`; forward on `DFTB2`; delegate `AM1Evaluator::launch`
- [ ] GREEN: full unit suite bit-identical (protocol-only change); Doxygen on every new public symbol

## Testing strategy

Unit only — this spec ships **no** regression of its own (pure interface; the
engine-slot capability's e2e witness is `esmm-qmmm-03`'s coupling regression, per the
one-regression-per-delivered-feature rule the engine swap belongs to the coupling).
Bit-for-bit parity between erased and direct evaluation is the load-bearing gate.

## Out of scope

- Coupling rewrite / mask column (`esmm-qmmm-03`)
- GPW/PW engines (seam: `GPWEvaluator` later derives `QmEngine` + `EmbeddingEngine`
  via `integrate_potential(cbg, nbasis, V_r)`, `collocate.h:101`)
- Any new physics; CUDA; Config bags / factories

## Notes for implementers

- **HARD GATE:** before touching `src/cpu/scf/*` (`mean_field_method.*`,
  `internal/gradient.cpp`, `internal/integral_derivatives.*`, `internal/point_charges.*`)
  or any CMakeLists, verify `git status` shows them clean. As of 2026-07-26 the
  operator's libcint-removal/C++20 work (~760 files) is uncommitted and includes them.
  If dirty: STOP and report. (`cpu/multiscale/` itself is clean.)
- `HF` stays a pure tag; `QmEngine` is NOT a Component — do not register engines into
  any Plan. Couplings reach them by direct call only (nested DriverNodes are
  second-class in three replay paths; see the Wave-1 specs).
- Multiple inheritance order: append the protocol bases after the existing ones;
  `MeanFieldMethod` is move-only — protocol bases must not break that.

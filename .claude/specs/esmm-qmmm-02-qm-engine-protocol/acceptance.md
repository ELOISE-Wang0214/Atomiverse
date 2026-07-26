---
spec: esmm-qmmm-02-qm-engine-protocol
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Protocol header complete and documented
    type: code
    pass_when: |
      cpu/component/qm_engine.h declares QmEngine (qm_engine_name/prepare/
      evaluate/as_embedding/reset_history) and EmbeddingEngine
      (set_embedding_sites/embedding_site_forces/set_site_self_energy) with
      Doxygen carrying units and the private-QM-frame contract; neither
      derives Component; no *Config, no make_* factory.
    status: pending
  - id: ac-002
    summary: All 7 Gaussian methods are engines via one CRTP implementation
    type: runtime
    pass_when: |
      tests/cpu/component/test_qm_engine.cpp probes as_embedding() non-null
      through QmEngine* for RHF, UHF, RKS (representative of the 7); the
      implementation lives once on MeanFieldMethod<D>, not per method.
    status: pending
  - id: ac-003
    summary: AM1 and DFTB2 are engines; AM1 embedding probe is null
    type: runtime
    pass_when: |
      AM1 and DFTB2 are reachable through QmEngine*; as_embedding() returns
      nullptr for both; DFTB2 forwards are one-line (evaluate→launch,
      prepare→build); AM1Evaluator::launch delegates to the solver evaluate
      (writeback code exists once — grep shows no duplicated energy/forces
      writeback between am1_evaluator.cpp and am1_solver.cpp).
    status: pending
  - id: ac-004
    summary: Erased evaluation is bit-identical to direct
    type: runtime
    pass_when: |
      For RHF on H2/STO-3G, energy via QmEngine::evaluate equals
      MeanFieldMethod::run bit-for-bit; forces_requested=false leaves the
      frame's force spans untouched and skips gradient work.
    status: pending
  - id: ac-005
    summary: reset_history restores cold-start behavior
    type: runtime
    pass_when: |
      After a converged evaluate, reset_history() followed by evaluate at the
      same geometry uses the cold iteration count (no ASPC/last-P reuse);
      unit test pins the iteration counts.
    status: pending
  - id: ac-006
    summary: No behavior change anywhere
    type: runtime
    pass_when: |
      Full unit suite and existing regressions pass unchanged; no existing
      public signature removed except as specified in esmm-qmmm-01.
    status: pending
out_of_scope:
  - Coupling rewrite / mask column (esmm-qmmm-03)
  - GPW/PW engine implementations
  - Any physics change; CUDA
---

# Acceptance — `esmm-qmmm-02-qm-engine-protocol`

**Done** means any of {RHF, UHF, ROHF, GHF, RKS, UKS, ROKS, AM1, DFTB2} can be
held and driven through one erased `QmEngine*`, with embedding capability
probeable, evaluation bit-identical to the direct path, and warm-start history
resettable — the interface the Wave-2 couplings consume.

## AC-001 — the protocol
One header, two protocols, neighbor_protocol.h precedent.

## AC-002 — CRTP once
Seven methods, one implementation site.

## AC-003 — heterogeneous engines
The protocol's reason to exist: three engine shapes, one slot.

## AC-004 — erasure adds nothing
Bit-for-bit parity is the guard against silent divergence between paths.

## AC-005 — MC posture hook
`reset_history` is what `reset_history_each_launch` (esmm-qmmm-03) calls.

## AC-006 — pure interface
This stage must be invisible to every existing user.

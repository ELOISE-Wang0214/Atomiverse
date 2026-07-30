# Iron law — no freestanding integral bag (2026-07-30)

**Rule (binding):** electronic-structure methods own integral state as **private
members of the concrete method** (`MeanFieldMethod` / `RHF` / …), rebuilt from
**Frame** when geometry moves. Engines borrow those members via a non-owning
`internal::IntegralView` (refs only — not storage).

## Public surface

```cpp
cpu::scf::RHF rhf{cpu::scf::BasisSet::STO3G /*, charge, tol, max_iter */};
rhf.run(frame);
auto g = rhf.forces(frame);
```

**Forbidden:**

- Freestanding integral bags constructed from Frame in user code, regressions,
  examples, apps, or docs
- Teaching any pack/workspace type as a user step
- `rhf.basis = …` after default construction (ctor params only)

## Internal

- Host matrices (S, Hcore, ERI, shells, …) are **expanded private fields** on
  the method (`ints_S_`, `ints_Hcore_`, …).
- Engines take `IntegralView` bound to those fields (Fock, gradient, post-HF,
  Mulliken, …).
- A **temporary** fill object may exist only inside `rebuild_integrals_` to
  move data into members; it is never a long-lived method member and never a
  public type.
- CUDA mirrors CPU: `cuda::scf::RHF{basis}` + Frame; Fock upload takes
  discrete matrices.

## Harness, not unit-test gates

Documented here and in CLAUDE.md. Do **not** add ctest/unit greps as the
primary gate for this style rule.

## Relation

- **no-public-config-bags** — same family: no bag parallel to Frame/members
- **frame-only-token** — Frame is the only public I/O
- **no public factories** — no `build_integrals` / pack factories

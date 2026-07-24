---
slug: scf-review-02-eri-hoists
criteria:
  - id: ac-001
    summary: mp2.h no longer declares ao2mo_*/mo_hcore; internal/ao2mo.h does
    type: code
    evaluator_hint: "grep gate on public vs internal header"
    pass_when: |
      grep -E 'ao2mo_(ovov|full)|mo_hcore' include/atomiverse/cpu/scf/mp2.h
      prints nothing, AND include/atomiverse/cpu/scf/internal/ao2mo.h exists
      and declares all three (ao2mo_ovov, ao2mo_full, mo_hcore) with their
      Doxygen @brief retained. MP2Result still declared in mp2.h.
    status: pending
  - id: ac-002
    summary: no source locally redefines the spin-orbital ERI helpers
    type: code
    evaluator_hint: "grep gate under src/cpu/scf"
    pass_when: |
      grep -rE 'struct SO\b|so_eri|inline (int|double) +(spat_|spin_|phys_|antisym_)\b'
      src/cpu/scf returns no matches; include/atomiverse/cpu/scf/internal/spin_orbital_eri.h
      holds the sole inline SpinOrbital + antisym_eri (struct and bit-packed
      overloads) in namespace ATV::cpu::scf.
    status: pending
  - id: ac-003
    summary: aufbau MP denominator guard defined once as mp_denominator
    type: code
    evaluator_hint: "grep gate: single throw site"
    pass_when: |
      the literal "non-negative orbital-energy denominator" appears exactly
      once under include/ + src/cpu/scf (inside mp_denominator in
      internal/spin_orbital_eri.h); mp2.cpp, mp3.cpp, mp4.cpp, ccsd.cpp call
      mp_denominator(...) at their same-spin/doubles denominators with no
      inline `denom >= 0.0` throw block remaining. cisd.cpp/adc.cpp guards
      left untouched (different semantics).
    status: pending
  - id: ac-004
    summary: spin-orbital ERI helper unit test passes
    type: code
    evaluator_hint: "ctest test_cpu_scf_spin_orbital_eri"
    pass_when: |
      test_cpu_scf_spin_orbital_eri builds and passes: antisym_eri struct and
      bit-packed overloads agree bitwise on the fixture ERI across all spin
      branches; mp_denominator returns the hard-coded negative value within
      1e-12 and throws std::invalid_argument on a non-negative denominator.
    status: pending
  - id: ac-005
    summary: ao2mo-leak architecture scrub registered and green
    type: code
    evaluator_hint: "ctest no_ao2mo_public_scrub"
    pass_when: |
      tests/architecture/no_ao2mo_public_scrub.sh is registered in
      tests/CMakeLists.txt and passes; both new tests appear in the ctest list.
    status: pending
  - id: ac-006
    summary: existing post-HF regressions reproduce cited values unchanged
    type: runtime
    evaluator_hint: "regressions target; post-HF scf families"
    pass_when: |
      building the `regressions` target and running the regressions/scf/
      post-HF cases each report PASS against their REFERENCE.md hard-coded
      literals with no tolerance relaxation, covering at least:
      h2o_mp2_sto3g, ump2-o2-sto3g, romp2-o2-sto3g, frozen-core-mp2-h2o,
      scs-mp2-h2o-sto3g, sos-mp2-h2o-sto3g, h2o_ri_mp2, mp2-grad-h2o-sto3g,
      ump2-grad-o2-sto3g, mp3-h2-sto3g, mp4-h2-sto3g, ccsd-h2-sto3g,
      ccsd-grad-h2o-sto3g, ccsd-t-h2o-sto3g, cc2-h2-sto3g, eom-ip-ccsd,
      eom-ea-ccsd, eom-ee-ccsd, eom-cc2-h2, eom-sf-ccsd-o2, cisd-h2co-sto3g,
      adc2-h2co-sto3g, adc2x-h2co-sto3g, ip-adc2-h2o-sto3g, ea-adc2-h2-sto3g,
      gf2-h2-sto3g, g0w0-h2-sto3g, bse-h2-sto3g, casscf-h2-sto3g,
      casci-h2-sto3g, nevpt2-h2-sto3g. No live third-party oracle invoked.
    status: pending
  - id: ac-007
    summary: full check + unit-test suite green after the refactor
    type: runtime
    evaluator_hint: "ctest --output-on-failure"
    pass_when: |
      a clean CPU configure+build (ATV_BUILD_HF=ON) plus full ctest run
      completes with every test passing and no new warnings/skips introduced.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 / ac-003 (code, grep gates).** Encode the two headline
  review findings: the public `mp2.h` no longer leaks the AO→MO transform
  (Item 12) and the four duplicated spin-orbital ERI helper copies plus the
  copy-pasted aufbau guard collapse to one definition each (Items 5 and 8).
  All three are statically greppable without building.
- **ac-004 (code).** The new single-function unit test proves the two
  `antisym_eri` overloads are numerically identical (justifying the "one
  canonical encoding, ordering preserved" decision) and that `mp_denominator`
  keeps the guard's throw behavior.
- **ac-005 (code).** The architecture scrub is wired into ctest so the
  layering/dedup invariants stay enforced against future drift.
- **ac-006 (runtime, the binding regression gate).** No new regression is
  added — a behavior-preserving refactor ships none; the delivery bar is that
  the existing post-HF families under `regressions/scf/` still reproduce their
  cited literature values, verified by `/mol:impl` at delivery.
- **ac-007 (runtime).** Full check + test suite green confirms no numeric or
  build regression across the whole tree.

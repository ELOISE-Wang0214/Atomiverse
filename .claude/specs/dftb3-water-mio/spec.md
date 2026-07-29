---
slug: dftb3-water-mio
status: approved
created: 2026-07-27
updated: 2026-07-29
domain: semiempirical / tight-binding
queue: tb-method-loop
method: dftb3
oracle: dftb+
depends_on: []
unblocks: [dftb3-analytic-forces]
---

# Spec: `dftb3-water-mio` — SCC-DFTB3 energy (third-order) on H2O

**Slug:** `dftb3-water-mio`
**Domain:** `cpu/semiempirical` tight-binding
**Method:** `dftb3`
**Oracle:** DFTB+ (family-locked; never OpenMOPAC / xtb)
**Status:** wave-3 resume dual-gate (product already on `dev` after
`68df18a8`..`148a17ed`; this tip re-verifies acceptance)

## Summary

**SCC-DFTB3** total energy product foothold: Frame in, energy out on
`frame.meta[kEnergy]` [eV], via Evaluate + `Semiempirical` Component.
Third-order charge fluctuation (Hubbard derivatives \(U^{\mathrm{d}}\),
\(\Gamma_{AB}\)) on static mio-style H/C/N/O SK tables.

## Domain basis

Gaus et al., *J. Chem. Theory Comput.* **7**, 931 (2011), DOI 10.1021/ct100684s
on top of Elstner DFTB2 (PRB 58, 7260 (1998)):

\[
E^{(2)}=\tfrac12\sum_{AB}\gamma_{AB}\,\Delta q_A\Delta q_B,\qquad
E^{(3)}=\tfrac13\sum_{AB}\Delta q_A^2\,\Delta q_B\,\Gamma_{AB},
\]

\[
V^{(3)}_C=\tfrac23\,\Delta q_C\sum_B\Gamma_{CB}\Delta q_B
         +\tfrac13\sum_B\Delta q_B^2\,\Gamma_{BC},
\]

\[
E=\mathrm{Tr}(\rho H^0)+E^{(2)}+E^{(3)}+E_{\mathrm{rep}}.
\]

SCC uses Anderson charge mixing (`mix_beta`, default 0.3); exhausting
`max_scc` throws. H–H reads transcribed mio-1-1 via
`internal::MioHHTable` (same as DFTB0/DFTB2).

## Public API (no Config, no factories)

```cpp
cpu::semiempirical::DFTB3 dftb3;
dftb3.scc_tol = 1e-10;
dftb3.max_scc = 50;
dftb3.mix_beta = 0.3;
dftb3.with_forces = false;
eval.on(evaluate::Phase::evaluate).add(std::move(dftb3));
```

## Geometry pin

Water off-eq (Å), cubic 20 Å: O(0,0,0), H(0.96,0,0), H(−0.24,0.93,0).

## Acceptance

1. Type `cpu::semiempirical::DFTB3` : `Semiempirical`.
2. Unit `test_cpu_dftb3` + `test_cpu_dftb3_scc` PASS (H2O finite; H2 reduces
   to DFTB2; Euler/FD on \(V^{(3)}\); SCC budget-independent; third-order
   gap commensurate with \(E^{(3)}\)).
3. Regression `reg_dftb3-water-mio`: pin
   \(E=-94.273801189131\) eV (tol \(10^{-6}\)); SCC residual \(<10^{-10}\);
   \(|E_3-E_2|\) commensurate with \(|E^{(3)}|\).
4. `REFERENCE.md`: `oracle: dftb+`, `oracle_version`, eV; Hamiltonian twin note.
5. SE-only cmake; scoped build only.
6. `perf_exempt: true` — no host DFTB+ binary for wall ratio.

## Non-goals

- Analytic forces → `dftb3-analytic-forces`.
- Absolute DFTB.org mio-1-1 identity on multi-element (SK fidelity).
- Spin / PBC / TD-DFTB / OpenMOPAC / xtb pins.

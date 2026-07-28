---
slug: dftb3-water-mio
status: approved
created: 2026-07-27
domain: semiempirical / tight-binding
queue: tb-method-loop
method: dftb3
oracle: dftb+
depends_on: []
unblocks: [dftb3-water-mio]
---

# Spec: `dftb3-water-mio` — SCC-DFTB3 energy (third-order) on H2O

**Slug:** `dftb3-water-mio`  
**Domain:** `cpu/semiempirical` tight-binding  
**Method:** `dftb3`  
**Oracle:** DFTB+ (family-locked; never OpenMOPAC / xtb)  
**Status:** approved for wave-1 implement

## Summary

Deliver **SCC-DFTB3** total energy as a product foothold: Frame in, energy
out on `frame.meta[kEnergy]` [eV], via Evaluate + `Semiempirical` Component.
DFTB2 (second-order SCC) already ships; this adds the **third-order charge
fluctuation** term (Hubbard derivatives \(U^\mathrm{d}\), \(\Gamma_{AB}\))
on the same static mio-style H/C/N/O SK tables.

## Domain basis

Elstner SCC-DFTB2 + Gaus DFTB3 third-order expansion
(Gaus et al., *J. Chem. Theory Comput.* **7**, 931 (2011),
DOI [10.1021/ct100684s](https://doi.org/10.1021/ct100684s)):

\[
E^{(2)} = \tfrac12\sum_{A,B}\gamma_{AB}(R_{AB})\,\Delta q_A\Delta q_B,
\]

\[
\Gamma_{AB} = \frac{\partial\gamma_{AB}}{\partial U_A}\,U_A^{\mathrm{d}},
\quad
E^{(3)} = \tfrac13\sum_{A,B}\Gamma_{AB}\,\Delta q_A^{2}\,\Delta q_B.
\]

Hamiltonian shift (atomic potential on Mulliken site \(A\)):

\[
V_A = \sum_B \gamma_{AB}\Delta q_B
  + \tfrac13\sum_B(\Gamma_{AB}+\Gamma_{BA})\Delta q_B^{2}
  + \tfrac13\sum_B(\Gamma_{BA}-\Gamma_{AB})\Delta q_A\Delta q_B.
\]

With eigenvalues of \(H_0 + V\):

\[
E = E_{\mathrm{band}} - E^{(2)} - E^{(3)} + E_{\mathrm{rep}}.
\]

\(\gamma_{AB}\) is the same Ohno–Klopman-like form as DFTB2; on-site
\(\gamma_{AA}=U_A\), \(\Gamma_{AA}=U_A^{\mathrm{d}}\).

**SK tables:** static `basis_data/dftb_mio_tables.h` (mio-**style** sketch,
not bit-identical DFTB.org mio-1-1). Absolute numerical identity with a
DFTB+ run on production mio-1-1 is **out of scope** for this slug
(owned by `dftb2-3ob-organic` / energy-pin SK fidelity). The formula and
module boundary follow DFTB+.

## Public API (no Config, no factories)

```cpp
cpu::semiempirical::DFTB3 dftb3;
dftb3.scc_tol = 1e-10;
dftb3.max_scc = 50;
dftb3.with_forces = false;  // energy foothold; forces optional FD like DFTB2
eval.on(evaluate::Phase::evaluate).add(std::move(dftb3));
```

Parameters are **public members**. No `DFTB3Config`, no `make_dftb3`.

## Geometry pin

Water, slightly off-equilibrium (Å), cubic 20 Å box (isolated):

| Atom | x    | y    | z   |
|------|------|------|-----|
| O    | 0.00 | 0.00 | 0.0 |
| H    | 0.96 | 0.00 | 0.0 |
| H    |−0.24 | 0.93 | 0.0 |

## Acceptance

1. Type `cpu::semiempirical::DFTB3` : `Semiempirical` under
   `include/atomiverse/cpu/semiempirical/dftb3.h`.
2. Unit test `test_cpu_dftb3`: H2O energy finite; rejects unsupported Z;
   third-order path differs from DFTB2 on same geometry when \(|\Delta q|>0\).
3. Regression `regressions/semiempirical/dftb3-water-mio/`:
   - `REFERENCE.md` with `oracle: dftb+`, `oracle_version`, units eV,
     Hamiltonian note (SK twin + Gaus formula).
   - Hard-coded energy gate (finite, non-zero; pin ATV Hamiltonian twin
     once measured; document that DFTB.org mio-1-1 absolute pin is deferred).
   - Optional: \(|E_{\mathrm{DFTB3}}-E_{\mathrm{DFTB2}}| > \varepsilon\) so
     third-order is active on water.
4. SE-only cmake; scoped build of `atomiverse`, `test_cpu_dftb3`,
   `reg_dftb3-water-mio` only.
5. Perf: optional wall \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\le 10\);
   `perf_exempt` if no offline DFTB+ binary on host.

## Non-goals

- Analytic forces (FD ok, same as current DFTB2).
- Spin, PBC, TD-DFTB, 3ob production tables.
- OpenMOPAC / xtb pins.
- `reg_am1_water_bomd` / minimize (MD/MINIMIZER off).

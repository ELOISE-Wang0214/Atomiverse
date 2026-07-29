---
slug: dftb2-kpoints
status: approved
created: 2026-07-29
revised: 2026-07-29
domain: semiempirical / dftb2
queue: tb-method-queue
oracle: dftb+
depends_on: [dftb2-periodic-gamma]
unblocks: []
---

# Spec: `dftb2-kpoints` — SCC-DFTB2 Monkhorst–Pack k-mesh beyond Γ

**Slug:** `dftb2-kpoints`  
**Method:** SCC-DFTB2 (MioSketch)  
**Family oracle:** **DFTB+ only** (never OpenMOPAC / xtb)  
**Delivery:** Evaluate + `cpu::semiempirical::DFTB2` Semiempirical component

## Problem

Γ-point periodic SCC-DFTB2 lattice sums ship (`reg_dftb2-periodic-gamma`).
Product cells need **k-mesh sampling beyond Γ**: Bloch-phased complex
\(H^0(\mathbf{k})\), \(S(\mathbf{k})\), k-weighted density / Mulliken Δq, and
band energy. Without this, 1-D / 3-D solids only have the Γ foothold.

## Goal

Public **member** k-mesh on `DFTB2` (no `*Config`, no `make_*`):

```cpp
cpu::semiempirical::DFTB2 dftb;
dftb.k_mesh = {4, 1, 1};          // Monkhorst–Pack per reciprocal axis
dftb.k_time_reversal = true;      // fold ±k (default)
dftb.with_forces = false;         // multi-k energy path (this slug)
eval.on(evaluate::Phase::evaluate).add(std::move(dftb));
```

Default `k_mesh = {1,1,1}` keeps the existing **real Γ** path bit-compatible
with `reg_dftb2-periodic-gamma` and molecular pins.

## Physics (binding)

- Bloch AO lattice sum (ket image \(T = n_0\mathbf{a}_0+\cdots\)):
  \[
  H^0_{\mu\nu}(\mathbf{k})=\sum_T e^{2\pi i\,\kappa\cdot n}\,H^0_{\mu\nu}(0,T),\quad
  S_{\mu\nu}(\mathbf{k})\ \text{likewise}
  \]
  with fractional \(\kappa\) from CP2K / DFTB+ Monkhorst–Pack
  \(\kappa_d=(2r-N_d-1)/(2N_d)\), \(r=1..N_d\).
- SCC shift is **k-independent**: \(\varepsilon_A=\sum_B\gamma_{AB}\Delta q_B\),
  \(H_{\mu\nu}(k)=H^0_{\mu\nu}(k)+\tfrac12(\varepsilon_A+\varepsilon_B)S_{\mu\nu}(k)\).
- Closed-shell occupancy: \(n_\mathrm{occ}=N_e/2\) filled bands at every k
  (no fractional metallic smearing in this foothold).
- Density / charges: \(P^k_{\mu\nu}=2\sum_i^{\mathrm{occ}}C_{\mu i}C^*_{\nu i}\),
  \(\Delta q_A=\sum_k w_k\sum_{\mu\in A}\mathrm{Re}[(P^k S^k)_{\mu\mu}]-q^0_A\).
- Energy: \(E=\sum_k w_k\,\mathrm{Re}\,\mathrm{Tr}(P^k H^{0k})+E^{(2)}+E_\mathrm{rep}\)
  (same double-counting convention as the Γ path).
- Real-space \(\gamma\) and repulsive remain **lattice-summed** (Γ-style),
  independent of k-mesh (screened, cutoff-truncated — no Ewald).

## Scope IN

| Item | Detail |
| ---- | ------ |
| API | `std::array<int,3> k_mesh{{1,1,1}}`, `bool k_time_reversal=true` |
| SK | MioSketch only (ThreeOb multi-k throws) |
| Spin | restricted closed shell only (open shell multi-k throws) |
| Forces / stress | multi-k throws if `with_forces` or `with_stress` |
| Γ identity | `{1,1,1}` uses the existing real Γ solver (pin-safe) |
| Regression | `regressions/semiempirical/dftb2-kpoints/` + offline DFTB+ twin |
| Unit tests | mesh generation / Γ vs multi-k separation / `{1,1,1}` ≡ Γ |

## Scope OUT

- Analytic nuclear forces / stress at general k
- Metallic smearing / fractional occupation
- Spin-polarized multi-k
- ThreeOb / D3 under multi-k
- Ewald long-range Coulomb for \(\gamma\)
- Live DFTB+ in ctest
- CUDA

## Oracle pin policy

```text
oracle: dftb+
oracle_version: DFTB+ family (offline twin; ATV mio-style SK Hamiltonian pin)
```

No DFTB+ binary on the build host → absolute energy is the **ATV Hamiltonian
twin** at fixed geometry / mesh (same hygiene as `dftb2-periodic-gamma`).
Family numerical oracle remains DFTB+. Input twin under
`tools/tb_oracle/inputs/dftbplus/dftb2-kpoints/`.

## Geometry (regression)

Same 1-D C dimer cell as Γ: two C at (0,0,0) and (1.4,0,0) Å;
orthorhombic \(L_x=2.8\), \(L_y=L_z=12.0\) Å. Mesh `{4,1,1}` with time-reversal.

## Acceptance summary

See `acceptance.md`. Dual gate: unit tests in `test_cpu_dftb2` +
`run_reg_dftb2-kpoints`. SE-only cmake. Perf: `perf_exempt` until offline
DFTB+ wall is available.

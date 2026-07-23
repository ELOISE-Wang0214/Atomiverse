---
slug: pw-spin-polarized-o2
status: approved
created: 2026-07-23
track: cp2k-method-loop
domain: pw
method: PW collinear spin-polarized KS (LSD)
geom: O2
basis: NC-UPF (SG15 ONCV PBE oxygen)
observable: e_tot / magnetization
grilled: false
---

# pw-spin-polarized-o2

## Summary

Ship **collinear spin-polarized plane-wave KS** (nspin=2, fixed occupations) on
`cpu::pw::PW`, with a textbook **O₂ triplet** (³Σg⁻, multiplicity 3) regression
against an offline CP2K/literature pin for total energy and magnetization.
Today’s PW path is spin-unpolarized only (factor-of-2 occupations). k-point
periodic UKS metals remain a later slug (`pw-uks-kmesh-fe`).

## Domain basis

Collinear spin-DFT (von Barth–Hedin / LSD):

- Separate α/β Kohn–Sham orbitals and densities ρ↑, ρ↓.
- Hartree from total density ρ = ρ↑ + ρ↓.
- Spin-dependent XC potentials Vxc↑(ρ↑,ρ↓), Vxc↓(ρ↑,ρ↓) via existing
  `cpu::scf::evaluate_xc` (LSDA / PBE-σ=0 path already used by PW).
- Fixed occupations for non-metallic O₂: N↑ − N↓ = 2S = multiplicity − 1.
- Magnetization M = N↑ − N↓ (electrons); for O₂ triplet M = 2.

O₂ ground state is the triplet ³Σg⁻ (Herzberg; standard computational test).
Valence with NC oxygen PP: Z_val = 6 per O → 12 valence electrons;
N↑ = 7, N↓ = 5 at multiplicity 3.

## Design

### Public API (no `*Config` bag, no factories)

On `cpu::pw::PW` members (defaults preserve today’s restricted path):

| Member | Default | Meaning |
| ------ | ------- | ------- |
| `spin_polarized` | `false` | When true, collinear nspin=2 KS |
| `multiplicity` | `1` | 2S+1; used only if `spin_polarized` |

On `PWResult` (after `run`):

| Field | Units | Meaning |
| ----- | ----- | ------- |
| `energy` / `energy_ha` | Ry / Ha | Total energy (existing) |
| `magnetization` | e⁻ | N↑ − N↓ from fixed occupations (and density integral cross-check) |
| `converged` / `iterations` | — | Existing |

Usage:

```cpp
cpu::pw::PW pw;
pw.upf_path = "tests/data/pw/O_ONCV_PBE-1.2.upf";
pw.ecut_ry = 30.0;
pw.k_mesh = 1;                 // Γ for molecular cell
pw.spin_polarized = true;
pw.multiplicity = 3;           // O2 triplet
pw.build(frame);
const auto& r = pw.run(frame); // r.magnetization ≈ 2
```

`PWEvaluator` reuses the owned `PW` method (no API bag). Forces/stress for
spin-polarized are **not** required for this slug (energy + mag only).

### Internal

- `setup_pw_problem` (or `PW::run`) computes `n_alpha`, `n_beta` from
  `n_elec` and `multiplicity` when spin-polarized; throws if parity invalid.
- Dual densities ρ↑/ρ↓ (G-space), dual Vxc / Veff, dual Davidson solves per
  k-point with occupations 1 (no spin-degeneracy factor of 2).
- Restricted path (`spin_polarized == false`) must remain bit-for-bit equivalent
  to current behaviour (existing unit tests / regressions still pass).
- Ship oxygen NC-UPF under `tests/data/pw/` for offline use (no network at
  ctest).

### Reuse decision

- `reuse` Poisson, Davidson, FFT, UPF, Pulay, `evaluate_xc`.
- `pattern` molecular UKS occupation split (Nα, Nβ from multiplicity).
- `new` collinear spin channel on PW workspace + SCF loop branch.

## Files to create or modify

- `.claude/specs/pw-spin-polarized-o2/*` (this)
- `include/atomiverse/cpu/planewave/pw_solver.h` — members + result fields
- `include/atomiverse/cpu/planewave/internal/pw_workspace.h` — spin buffers
- `include/atomiverse/cpu/planewave/internal/{density,xc,veff}_builder.h` as needed
- `src/cpu/planewave/pw_solver.cpp` — spin SCF branch
- `src/cpu/planewave/internal/*` — spin density/XC/energy
- `tests/data/pw/O_ONCV_PBE-1.2.upf` — NC oxygen PP
- `tests/cpu/planewave/test_pw_spin.cpp` — unit (API + mag identity)
- `regressions/pw/pw-spin-polarized-o2/{reg_*.cpp,REFERENCE.md}`
- `regressions/CMakeLists.txt` — register case
- `tools/cp2k_oracle/inputs/` — optional twin input for offline re-pin

## Tasks

- [ ] Spec + acceptance (this directory)
- [ ] Unit test: `spin_polarized` + `multiplicity` members; O₂ toy/UPF mag ≈ 2
- [ ] Implement collinear spin path; keep restricted default unchanged
- [ ] Offline pin in `REFERENCE.md` (CP2K if available, else cited literature /
      fixed-occupation identity + energy gate with versioned method string)
- [ ] E2E regression `regressions/pw/pw-spin-polarized-o2`
- [ ] Build + run unit + regression; record perf gate or `perf_exempt`
- [ ] Commit on `wf/pw-spin-polarized-o2`

## Testing strategy

- **Unit** (`tests/`): construct API; Γ O₂ with small ecut converges with
  `|M − 2| < 1e-6` (fixed occupations); restricted path still finite energy.
- **Regression**: public API only; assert energy vs offline pin and
  magnetization = 2 within tol. No live CP2K in ctest.

## Out of scope

- k-point UKS metals / smearing (`pw-smearing-al`, `pw-uks-kmesh-fe`)
- Noncollinear / SOC
- Spin-polarized analytic forces / stress
- GPW UKS (`gpw-uks-o2-gamma`)
- CUDA PW spin
- Hybrid XC / meta-GGA

## Dual gate

1. **Correctness:** `REFERENCE.md` offline pin + regression PASS.
2. **Perf:** fair wall vs CP2K if binary available (≤10×); else
   `perf_exempt=true` with reason (no CP2K PW/NC-UPF twin on host).

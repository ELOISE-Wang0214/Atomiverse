---
slug: pw-uks-kmesh-fe
status: approved
created: 2026-07-23
track: cp2k-method-loop
domain: pw
method: k-point collinear UKS plane-wave
geom: O2 cubic cell (Fe metal deferred — needs spin+smearing)
basis: NC-UPF (SG15 ONCV PBE oxygen)
observable: e_tot / mag
grilled: false
---

# pw-uks-kmesh-fe

## Summary

Close the **k-point collinear UKS** gap on `cpu::pw::PW`: spin-polarized KS
(`nspin=2`, fixed occupations) must run with **Monkhorst–Pack `k_mesh ≥ 2`**,
not only Γ. Wave-1 shipped Γ-only O₂ spin (`pw-spin-polarized-o2`) and
restricted metallic smearing (`pw-smearing-al`); this slug is the periodic UKS
ladder step with a multi-k mesh.

**Geometry:** Queue allows **Fe or O₂ cell**. Bulk Fe metal needs collinear
spin **plus** metallic smearing (currently rejected:
`spin + smearing not supported together`) and a Fe NC-UPF with Z>Ar toy
valence. This slug uses the **O₂ cubic cell** with `k_mesh = 2` (fixed
occupations, mag = 2) so the dual gate is achievable offline without inventing
a new public API. Fe metal + spin+smearing is a follow-on (not this acceptance).

## Domain basis

Collinear spin-DFT (same as `pw-spin-polarized-o2`) with Brillouin-zone
sampling:

- Dual α/β orbitals and densities at **every** k-point of the MP mesh.
- Hartree from total ρ = ρ↑ + ρ↓; spin XC via `step_spin`.
- Fixed occupations: N↑ − N↓ = multiplicity − 1 (O₂ triplet → M = 2).
- k-point weights from `generate_mp_mesh` (existing); spin path already loops
  `ik` over owned k-points — product work is pin + regression + unit proof that
  `k_mesh ≥ 2` converges with correct mag, not a new Config bag.

## Design

### Public API (no `*Config`, no factories)

**No new members required** if the existing collinear path already accepts
`k_mesh > 1` (it does). Public surface remains:

| Member / field | Role |
| -------------- | ---- |
| `PW::spin_polarized`, `multiplicity` | collinear UKS (existing) |
| `PW::k_mesh` | MP mesh size per axis (existing; **must be ≥ 2** in regression) |
| `PWResult::magnetization`, `energy` / `energy_ha` | observables |

Usage:

```cpp
cpu::pw::PW pw;
pw.upf_path = "tests/data/pw/O_ONCV_PBE-1.2.upf";
pw.ecut_ry = 10.0;
pw.k_mesh = 2;                 // multi-k UKS (this slug)
pw.spin_polarized = true;
pw.multiplicity = 3;
pw.build(frame);
const auto& r = pw.run(frame); // r.magnetization == 2; e_tot vs pin
```

### Internal

- Reuse spin SCF branch in `PW::run` (dual Davidson × k).
- Fix only if multi-k spin fails correctness/convergence under the pin settings.
- Restricted path and Γ-spin regression must keep working.
- Do **not** enable spin+smearing in this slug.

### Reuse decision

- `reuse` spin channel, MP mesh, O NC-UPF, Poisson/FFT/Davidson.
- `pattern` Γ-UKS occupations + multi-k restricted loop structure.
- `new` offline multi-k UKS pin + unit/e2e gates only (API already present).

## Files to create or modify

- `.claude/specs/pw-uks-kmesh-fe/*` (this + acceptance)
- `tests/cpu/planewave/test_pw_spin.cpp` — multi-k mag identity unit
- `regressions/pw/pw-uks-kmesh-fe/{reg_*.cpp,REFERENCE.md}`
- `regressions/CMakeLists.txt` — register case
- `src/cpu/planewave/pw_solver.cpp` — only if multi-k spin bug found
- `tools/cp2k_oracle/inputs/` — optional CP2K twin note (not ctest)

## Tasks

- [x] Spec + acceptance (this directory)
- [x] Unit: `k_mesh≥2` + spin_polarized + multiplicity 3 → |M−2| unit tol
- [x] Offline pin in REFERENCE (method string; CP2K if available else frozen ATV)
- [x] E2E regression public API only
- [x] Build + run unit + regression; perf gate or `perf_exempt`
- [x] Commit on `wf/pw-uks-kmesh-fe`

## Testing strategy

- **Unit** (`tests/`): multi-k O₂ (tiny ecut / short budget) still yields
  |M − 2| ≈ 0 from fixed occupations; no full e2e in unit suite.
- **Regression:** SCF converges; |mag − 2| ≤ 1e−9; |E − pin| ≤ tol; ≤ 60 s.
  No live CP2K in ctest.

## Out of scope

- Spin + metallic smearing (Fe bulk metal occupations)
- Fe NC-UPF product pin / valence_from_z(Z=26) without UPF override
- Noncollinear / SOC
- Spin-polarized analytic forces / stress
- GPW k-points (`gpw-kpoints-si`)
- CUDA PW spin

## Dual gate

1. **Correctness:** offline pin + regression PASS (mag + e_tot + converged).
2. **Perf:** ≤10× vs CP2K twin wall if available; else `perf_exempt` with reason.

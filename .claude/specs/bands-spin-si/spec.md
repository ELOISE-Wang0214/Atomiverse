# bands-spin-si — collinear spin-polarized PW band structure

## Goal

Lift the `Bands::compute` spin-polarized foothold throw and ship product
**collinear spin band structure** (NSCF on frozen \(V_{\mathrm{eff},\sigma}\))
for plane-wave DFT, pinned on diamond Si with multiplicity 1 (\(M=0\)).

Sibling product: closed-shell `reg_bands-dos-si`. Residual after that
foothold: spin path threw. Metallic DOS remains a separate slug
(`dos-metallic-si`).

## Scope

| In | Out |
| -- | --- |
| `cpu::pw::Bands::compute` with `method.spin_polarized == true` | Metallic smear bands (hold `dos-metallic-si` / thin sibling) |
| Fixed occupations, collinear UKS/LSDA | Spin-polarized DOS (`Dos::compute`) |
| NSCF eigenvalues per spin on diamond path | Berry / Wannier spin (other slugs) |
| Regression `regressions/pw/bands-spin-si` | Live CP2K in ctest |
| Offline CP2K oracle note under `tools/cp2k_oracle/inputs/` | CUDA ports |

## Public API (no `*Config`, no factories)

```cpp
cpu::pw::PW method;
method.spin_polarized = true;
method.multiplicity = 1;   // M = N↑ − N↓ = 0 for bulk Si
method.xc = cpu::scf::XCType::LSDA;
// ecut, k_mesh, upf, …

cpu::pw::Bands bands;
bands.n_bands = 8;         // per spin channel
bands.set_diamond_path(4);
auto br = bands.compute(method, frame);
// br.nspin == 2
// br.eigenvalues_ry       → spin-up [ik][ib] [Ry]
// br.eigenvalues_beta_ry  → spin-down [ik][ib] [Ry]
// br.magnetization        → N↑ − N↓
// br.gap_ry from joint fixed-occ VBM/CBM (α/β)
```

Parameters remain members of `PW` / `Bands`. Frame in; owned result out.

## Algorithm

1. **SCF (collinear fixed-occ):** same shape as product UKS/DFPT spin:
   \(N_\uparrow=(N_e+M)/2\), \(N_\downarrow=(N_e-M)/2\); per iteration
   Poisson on \(\rho=\rho_\uparrow+\rho_\downarrow\), `XCBuilder::step_spin`,
   \(V_{\mathrm{eff},\sigma}=V_{\mathrm{loc}}+V_H+V_{\mathrm{xc},\sigma}\),
   Davidson per (spin, k), linear mix on spin densities.
2. Freeze \(\rho_\sigma\) and recompute final \(V_{\mathrm{eff},\sigma}\).
3. `rebuild_k_tables` for the high-symmetry path (density frozen).
4. **NSCF:** for each spin \(\sigma\), set `ws.V_eff_r ← V_{\mathrm{eff},\sigma}\)
   and diagonalize `n_bands` states at every path \(\mathbf{k}\).
5. Gap: VBM = max occupied \(\varepsilon_{n\sigma}(\mathbf{k})\), CBM = min
   empty (fixed-occ cutoffs per spin). Report Ry + path indices.

Reject `smearing != none` (metallic residual is out of scope). Reject
invalid multiplicity / electron parity.

## Acceptance (summary)

See `acceptance.md`. Dual gate: correctness pins in
`regressions/pw/bands-spin-si/REFERENCE.md`; perf ≤ 10× CP2K twin or
`perf_exempt` with reason.

## Files

| Path | Change |
| ---- | ------ |
| `include/.../bands_dos.h` | Spin fields on `BandsResult`; docs |
| `src/.../bands_dos.cpp` | Collinear SCF + per-spin NSCF path |
| `tests/cpu/planewave/test_pw_bands_dos.cpp` | Unit: throw policy / path unchanged |
| `regressions/pw/bands-spin-si/*` | E2E + REFERENCE |
| `tools/cp2k_oracle/inputs/bands-spin-si.note` | Offline twin note |
| `regressions/CMakeLists.txt` | Register `reg_bands-spin-si` |

## Non-goals

- Spin DOS productization
- Metallic occupation on the band path
- Absolute eigenvalue match to CP2K GTH (PP family differs; character +
  literature Γ→X window + \(M=0\) + α≈β pins)

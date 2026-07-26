# Spec: bands-metallic-si

**Slug:** `bands-metallic-si`  
**Method:** PW band structure with metallic occupation smearing  
**Oracle:** CP2K (offline pin only)  
**Branch:** `wf/bands-metallic-si`  
**Domain:** `regressions/pw/bands-metallic-si`

## Problem

`cpu::pw::Bands::compute` rejects metallic occupations:

```text
Bands::compute: metallic smearing not supported for bands foothold
```

Product already ships:

| Capability | Status |
| ---------- | ------ |
| Insulator bands / DOS | `reg_bands-dos-si` ✓ |
| Collinear spin bands | `reg_bands-spin-si` ✓ |
| Metallic DOS (FD on denser MP) | `reg_dos-metallic-si` ✓ |
| Metallic SCF (Al energy) | `reg_pw-smearing-al` ✓ |

Residual: **high-symmetry-path bands after a metallic SCF density**.

## Goal

Enable restricted (spin-unpolarized) `Bands::compute` when
`method.smearing ∈ {fermi_dirac, methfessel_paxton, gaussian}`:

1. SCF with fractional occupations (`run_scf_smear_restricted`) on
   `method.k_mesh`.
2. Freeze \(V_\mathrm{eff}\); NSCF diagonalize on the k-path
   (`set_diamond_path` / `kfrac`).
3. Report eigenvalues [Ry], integer-fill spectrum gap (when empty bands
   exist), and \(E_F\) [Ry] from the SCF free-energy result (path is not
   a BZ mesh — do not re-bin \(E_F\) from path weights).

**Out of scope this slug:** collinear spin + metallic bands (keep throw);
spin DOS (`dos-spin-si`).

## Public API (no `*Config`, no factories)

```cpp
cpu::pw::PW method;
method.smearing = Smearing::fermi_dirac;
method.degauss = 0.02;   // Ry
method.extra_bands = 4;
// ecut, k_mesh, upf, xc=LSDA …

cpu::pw::Bands bands;
bands.n_bands = 8;
bands.set_diamond_path(4);
const BandsResult br = bands.compute(method, frame);
// br.eigenvalues_ry, br.e_fermi_ry, br.gap_ry, br.scf_converged
```

`BandsResult` gains `e_fermi_ry` (0 for fixed occupations).

## Twin

Diamond Si (gapped; FD occupations near-integer) — same geometry class as
`dos-metallic-si` / `bands-dos-si`. Product gate is **path enabled under
smear** + open gap + finite \(E_F\) near the gap, not a vanishing metal gap.

## Acceptance (see acceptance.md)

Dual gate: unit helpers stay green; regression converges + pins; offline
CP2K note under `tools/cp2k_oracle/inputs/`. `perf_exempt` if no CP2K
binary (same host policy as siblings).

## Files

| Path | Change |
| ---- | ------ |
| `include/.../bands_dos.h` | Document metallic path; `e_fermi_ry` on result |
| `src/cpu/planewave/bands_dos.cpp` | Smear SCF branch for restricted bands |
| `tests/cpu/planewave/test_pw_bands_dos.cpp` | Unit: default `e_fermi_ry`, reject spin+smear |
| `regressions/pw/bands-metallic-si/` | E2E + REFERENCE.md |
| `tools/cp2k_oracle/inputs/bands-metallic-si.note` | Offline twin note |
| `regressions/CMakeLists.txt` | Register `reg_bands-metallic-si` |

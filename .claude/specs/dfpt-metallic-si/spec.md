# Spec — `dfpt-metallic-si`

## Goal

Productize **Γ-point DFPT phonons** for the plane-wave method when the
ground state uses **metallic occupation smearing** (Fermi–Dirac /
Methfessel–Paxton / Gaussian). Today `dfpt_gamma_force_constants` and
`dfpt_q_force_constants` throw on `method.smearing != none`.

Scope: **restricted** (closed-shell) PW + pure LSDA or pure GGA + UPF.
Collinear spin + smear and finite-q + smear stay rejected (follow-on:
`dfpt-spin-si` covers spin fixed-occ; true-metal Al DFPT can follow).

## Public API

No new public types. Callers set members on existing objects:

```cpp
cpu::pw::PW pw;
pw.upf_path = "Si.upf";
pw.ecut_ry = 10.0;
pw.k_mesh = 2;
pw.xc = cpu::scf::XCType::LSDA;
pw.smearing = cpu::pw::Smearing::fermi_dirac;
pw.degauss = 0.02;       // Ry
pw.extra_bands = 4;

cpu::pw::DfptPhonon dfpt;
auto r = dfpt.compute(pw, frame);  // ω [cm⁻¹], C [eV/Å²]
```

Parameters remain members — no `*Config`, no factories.

## Algorithm

1. Ground-state SCF with the same metallic occupations as `PW::run`
   (shared μ, fractional \(f_{nk}\in[0,2]\), `extra_bands` empty states).
2. Freeze \(V_{\mathrm{eff}}[n]\); re-diagonalize every k with
   \(n_b = \lceil N_e/2\rceil + \mathrm{extra\_bands}\).
3. Fill \(f_{nk}\) from the frozen spectrum (same `fill_occupations`).
4. Sternheimer response only for bands with \(f > f_{\min}\)
   (\(f_{\min}=10^{-8}\)); projector \(P_c\) removes those bands.
5. Weight density response and force-constant terms by \(f_{nk}\)
   (replace fixed \(f=2\)).
6. No explicit \(\delta\mu\) Fermi-shift term in v1 (gapped Si twin
   negligible; true metals may need it later for ASR).

Cite: Baroni et al., *Rev. Mod. Phys.* **73**, 515 (2001);
de Gironcoli, *Phys. Rev. B* **51**, 6773 (1995).

## Acceptance

| Gate | Criterion |
| ---- | --------- |
| Unit | `DfptPhonon` accepts FD smear (no throw on UPF guard path after XC); still throws on spin+smear, hybrid, empty UPF, finite-q+smear |
| Regression | `regressions/response/dfpt-metallic-si` Si Γ twin with FD: acoustic ≈ 0, optical window, DFPT \(C_{00}\) vs offline energy-FD pin |
| Dual gate | Correctness pins + `within_10x` or `perf_exempt` |
| Arch | No public Config/factory; Frame I/O only; public energy still eV on Frame (spectra on result) |

## Out of scope

- Collinear spin + metallic smear DFPT
- Finite-q / incommensurate Sternheimer + smear
- Explicit \(\delta\mu\) occupation response for metals at Γ
- Live CP2K in ctest (offline pins only)

## Deps

- `dfpt-spin-si` ✓ (spin fixed-occ product)
- `pw-smearing-al` ✓ (PW metallic occupations)
- `dfpt-phonon-si` ✓ (Γ LSDA insulator baseline)

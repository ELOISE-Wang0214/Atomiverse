# Spec — pw-uks-fe-smearing

**Slug:** `pw-uks-fe-smearing`  
**Loop:** CP2K capability (`cp2k-method-loop`)  
**Method domain:** plane-wave (CPU)  
**Branch:** `wf/pw-uks-fe-smearing`  
**Parent:** `dev`

## Problem

Roadmap / queue: bulk **Fe** needs **collinear UKS and metallic smearing
together**. Today the pieces ship separately:

| Capability | Regression | Gap |
| ---------- | ---------- | --- |
| Restricted metallic smearing | `pw-smearing-al` | no spin |
| Collinear UKS fixed occupations | `pw-spin-polarized-o2`, `pw-uks-kmesh-fe` | no smear |
| Combined | — | `PW::run` throws if `spin_polarized && smearing != none` |

`pw-uks-kmesh-fe` deliberately used O₂ multi-k fixed fillings because Fe metal
is the spin+smear product case (REFERENCE there names this follow-on).

## Goal

Enable **collinear spin-polarized plane-wave KS with metallic occupations** on
`cpu::pw::PW`, pin **BCC Fe** (NC-UPF) for free energy + self-consistent
magnetization + electron count, and dual-gate offline vs a CP2K algorithm twin
(documentation only — GTH-GPW ≠ NC-UPF absolute energy).

## Public API (no `*Config`, no factories)

Existing members on `cpu::pw::PW` already cover the knobs:

```cpp
cpu::pw::PW pw;
pw.spin_polarized = true;
pw.multiplicity   = /* seeds initial mag split only when smearing is on */;
pw.smearing       = Smearing::fermi_dirac; // or mp / gaussian
pw.degauss        = 0.02; // Ry
pw.extra_bands    = 4;
pw.run(frame);
// PWResult::{energy, energy_ha, magnetization, e_fermi, n_elec_occ, e_entropy}
```

**Semantics when both spin and smearing are on:**

1. **Shared chemical potential** \(\mu\) over spin-up and spin-down spin-orbitals
   (QE `nspin=2` metallic default).
2. Occupations \(f_{\sigma kb}\in[0,1]\); \(N_e = \sum_{\sigma kb} w_k f\).
3. Magnetization \(M = N_\uparrow - N_\downarrow\) is **self-consistent** (not
   forced to `multiplicity-1`). Multiplicity only seeds the initial
   \(\rho_\uparrow/\rho_\downarrow\) split (same as fixed-occ path).
4. Total energy is the **Mermin free energy** \(F = E_{\mathrm{KS}} - TS\)
   (entropy summed over spin-orbitals).
5. `extra_bands` adds empty states **per spin** beyond \(\lceil N_e/2\rceil\).

No new public types. No `make_*` / `*_config`. Frame remains sole I/O.

## System (regression pin)

| Quantity | Value |
| -------- | ----- |
| Material | BCC Fe, 1 atom primitive cell |
| \(a_{\mathrm{cubic}}\) | 2.87 Å (primitive BCC vectors from cubic) |
| PP | SG15 ONCV PBE Fe `tests/data/pw/Fe_ONCV_PBE-1.2.upf` (\(Z_{\mathrm{val}}=16\)) |
| XC | PBE (default) |
| \(E_{\mathrm{cut}}\) | modest (algorithm pin; ≤ 60 s wall) — record exact in REFERENCE |
| k-mesh | Monkhorst–Pack (record exact; prefer ≥ Γ with multi-k if wall allows) |
| Smearing | Fermi–Dirac, `degauss` [Ry] recorded in REFERENCE |
| Multiplicity seed | ferromagnetic seed (e.g. 5 → target 2S=4 for initial split) |

## Offline pins (REFERENCE.md)

1. **Hard gates (ATV fixed settings):**
   - SCF converges.
   - \(|N_{\mathrm{occ}} - 16| \le 10^{-4}\).
   - Finite \(E_F\) (not absurd).
   - Free energy \(F\) [Ry or Ha] within stated tol of frozen ATV offline pin.
   - Magnetization \(M > 0\) (ferro seed) with loose band (Fe bulk ~2 μB order;
     low-ecut pin may undershoot — document).
2. **CP2K twin (algorithm / documentation, not absolute E):**
   - Quickstep UKS + `&SMEAR METHOD FERMI_DIRAC` + Fe GTH (or note if binary
     unavailable) under `tools/cp2k_oracle/inputs/`.
   - Absolute energy **not** compared (same class as `pw-smearing-al` /
     `pw-uks-kmesh-fe`).

## Implementation sketch

1. Remove the hard error in `PW::run` for spin+smear.
2. Extend internal smearing helpers for **spin-orbital** occupations
   \(f\in[0,1]\) (shared \(\mu\), entropy without closed-shell ×2).
3. New SCF branch: dual-spin Davidson → stash evals/wfc → Fermi search over
   both spins → fractional \(\rho_\uparrow,\rho_\downarrow\) → mix → Mermin energy.
4. Keep restricted-smear and fixed-occ spin paths bit-identical for existing regs.
5. Unit tests: spin-orbital Fermi recovers \(N_e\); public members still no Config.
6. Regression `regressions/pw/pw-uks-fe-smearing/` + CMake target.
7. In-tree Fe UPF (SG15, same license class as O_ONCV).

## Out of scope

- Noncollinear / SOC magnetism.
- Constrained magnetization (two chemical potentials).
- GPW spin+smear productization.
- CUDA PW spin or smear.
- Production-quality Fe \(E_{\mathrm{cut}}\) convergence study.
- Absolute energy twin vs CP2K/QE.

## Acceptance (see `acceptance.md`)

Green when dual gate passes and commit lands on `wf/pw-uks-fe-smearing`.

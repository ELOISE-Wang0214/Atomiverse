---
slug: gpw-metallic-smearing
status: approved
created: 2026-07-26
domain: gpw
queue: cp2k-method-loop
depends_on: []
unblocks: [gpw-metallic-smearing]
grilled: false
---

# Spec: `gpw-metallic-smearing` — Metallic occupation smearing on GPW

**Slug:** `gpw-metallic-smearing`  
**Domain:** GPW DFT (CPU), multi-k closed-shell RKS  
**Queue:** CP2K method loop  
**Status:** approved for impl

## Summary

Deliver **Fermi–Dirac / Methfessel–Paxton metallic occupation smearing** on the
product GPW path, mirroring the already-shipped PW surface
(`cpu::pw::Smearing` + `degauss` + `extra_bands` + Mermin free energy). Today
GPW multi-k and Γ SCF fill **fixed integer occupations** (`n_occ` double-filled
or `n_alpha`/`n_beta`). Metals require fractional occupations and a free-energy
functional \(F = E_{\mathrm{KS}} - TS\).

CP2K counterpart: `&SCF &SMEAR METHOD FERMI_DIRAC|METHFESSEL_PAXTON` with
`ELECTRONIC_TEMPERATURE` / `ELECTRONIC_TEMPERATURE [K]` and `ADDED_MOS`
(`tools/cp2k_oracle/inputs/al_fermi_dirac.inp` documents the Al twin; ATV pin
uses **BCC Na** with in-tree `GTH-PBE-q9` + `DZVP-MOLOPT-SR-GTH` because Al is
not yet in `gpw_tables.h` and GPW cells are orthorhombic-only).

## Domain basis

**Spin-unpolarized (RKS) band occupation** with shared chemical potential \(\mu\):

\[
f_b(\varepsilon) = 2\,f_{\mathrm{spinless}}\!\left(\frac{\varepsilon-\mu}{\sigma}\right),
\quad f_{\mathrm{spinless}}\in[0,1],
\]

with \(\sigma =\) `degauss` converted **Ry → Ha** (same convention as
`cpu::pw::PW`). Target electron count \(N_e = 2\,n_{\mathrm{occ}}\).

**Fermi level** \(\mu\) solves \(\sum_{k,b} w_k f_b(\varepsilon_{kb}) = N_e\)
(reusing `cpu::pw::internal::find_fermi_level` / `fill_occupations`).

**Density matrix** per k:

\[
P^k_{\mu\nu} = \sum_b f_{kb}\, C^k_{\mu b}\,(C^k_{\nu b})^\ast.
\]

**Mermin free energy** (reported as `GpwResult::e_total`):

\[
F = E_{\mathrm{KS}} - \sigma\sum_{k,b} w_k s(\varepsilon_{kb}),
\]

with entropy weight \(s\) from `smear_entropy_band` (FD / MP / Gaussian).

Reuse **only** the pure occupation math in
`cpu/planewave/internal/smearing.h` (shared physics, no PW solver dependency
beyond that internal header — already unit-tested by `test_pw_smearing`).

## Public API (no Config bags, no factories)

On `cpu::gpw::GpwDft` (and optionally mirrored on `GPWEvaluator`):

| Member | Type | Meaning |
| --- | --- | --- |
| `smearing` | `gpw::Smearing` | `none` (default), `fermi_dirac`, `methfessel_paxton` |
| `degauss` | `double` | width [Ry], required `> 0` when smear on |
| `mp_order` | `int` | Methfessel–Paxton order \(N\ge 1\) |
| `extra_bands` | `int` | API parity with CP2K `ADDED_MOS` / PW; GPW already diagonalizes the full AO spectrum so this is documentation + future band-cap |

On `GpwResult`:

| Field | Unit | Meaning |
| --- | --- | --- |
| `e_fermi` | Ha | \(\mu\) (0 when fixed occupations) |
| `e_entropy` | Ha | \(-\sigma S\) contribution folded into free energy (0 when fixed) |
| `n_elec_occ` | e | \(\sum w f\) after fill (should recover \(N_e\)) |

`e_total` remains the **free energy** \(F\) [Ha] when smearing is on (CP2K /
QE convention; PW reports free energy too).

## Scope / non-goals (this slug)

**In:** multi-k RKS GPW SCF with FD and MP-1 occupations + free energy; Γ RKS
may share the same fill path when `smearing != none`.  
**Out:** UKS + smear (throw fail-fast); OT + smear (throw); forces/stress with
smear (leave integer-occ force path; may throw or ignore for v1 — energy
regression only); adding Al basis tables; non-orthorhombic cells.

## Acceptance (see `acceptance.md`)

1. Unit: public members + FD fill recovers \(N_e\) on a toy spectrum
   (`test_cpu_gpw_smearing`).
2. Regression `reg_gpw-metallic-smearing`: BCC Na, multi-k, FD smear —
   converges; \(|n_{\mathrm{elec}}-N_e|<10^{-5}\); free-energy pin vs frozen
   first-green ATV; CP2K twin input documented offline.
3. Dual gate: within_10x or `perf_exempt` if no fair CP2K wall.

## References

- Mermin, *Phys. Rev.* **137**, A1441 (1965) — free-energy functional.
- Methfessel & Paxton, *Phys. Rev. B* **40**, 3616 (1989).
- CP2K Quickstep `&SMEAR`; ATV PW path `reg_pw-smearing-al`.

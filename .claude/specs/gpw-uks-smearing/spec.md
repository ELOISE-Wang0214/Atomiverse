---
slug: gpw-uks-smearing
status: approved
created: 2026-07-26
domain: gpw
queue: cp2k-method-loop
depends_on: [gpw-metallic-smearing]
unblocks: [gpw-uks-smearing]
grilled: false
---

# Spec: `gpw-uks-smearing` — Metallic occupation smearing on GPW UKS

**Slug:** `gpw-uks-smearing`  
**Domain:** GPW DFT (CPU), multi-k collinear UKS  
**Queue:** CP2K method loop  
**Status:** approved for impl

## Summary

Productize **Fermi–Dirac / Methfessel–Paxton metallic occupation smearing** on
the **collinear UKS** multi-k GPW path. Wave-6 closed multi-k **RKS** smear
(`gpw-metallic-smearing`); `GpwDft::build` still **rejects**
`n_alpha`/`n_beta` + `smearing != none`. PW already ships UKS+smear
(`reg_pw-uks-fe-smearing`). This slug removes the GPW UKS gate and wires
shared-\(\mu\) spin-orbital occupations into `run_gpw_scf_kpoints_uks`.

CP2K counterpart: Quickstep `LSD` + `&SCF &SMEAR METHOD FERMI_DIRAC` (or
`METHFESSEL_PAXTON`) with `ELECTRONIC_TEMPERATURE` and `ADDED_MOS`.

## Domain basis

**Collinear spin-orbital occupation** (no ×2), shared chemical potential \(\mu\)
over \(\alpha\) and \(\beta\) at all k:

\[
f_{\sigma kb} = f_{\mathrm{spinless}}\!\left(\frac{\varepsilon_{\sigma kb}-\mu}{\sigma}\right)\in[0,1],
\quad
\sum_{\sigma kb} w_k f_{\sigma kb} = N_e = n_\alpha^{\mathrm{seed}} + n_\beta^{\mathrm{seed}}.
\]

Reuse `cpu::pw::internal::{find_fermi_level,fill_occupations,smear_f_spinor,smear_entropy_spinor}`
with **`spin_degen = 1`** (same convention as PW collinear smear).

**Density** per spin/k:

\[
P^{\sigma,k}_{\mu\nu} = \sum_b f_{\sigma kb}\, C^{\sigma k}_{\mu b}\,(C^{\sigma k}_{\nu b})^\ast.
\]

**Magnetization** \(M = N_\uparrow - N_\downarrow\) is **self-consistent** under
smear; `n_alpha`/`n_beta` only seed the cold-start integer fill (and document
target \(N_e\)).

**Mermin free energy** (as `GpwResult::e_total`):

\[
F = E_{\mathrm{KS}} - \sigma\sum_{\sigma kb} w_k s_{\mathrm{spinor}}(\varepsilon_{\sigma kb}).
\]

Units: eigenvalues and \(\mu\) in **Ha** inside SCF; public `degauss` still **Ry**
(convert \(0.5\times\) as RKS smear / PW).

## Public API (no Config bags, no factories)

Existing members (already on `GpwDft` from RKS smear):

| Member | Type | Meaning |
| --- | --- | --- |
| `smearing` | `gpw::Smearing` | `none` / `fermi_dirac` / `methfessel_paxton` |
| `degauss` | `double` | width [Ry], `> 0` when smear on |
| `mp_order` | `int` | MP order \(N\ge 1` |
| `extra_bands` | `int` | API parity (full AO spectrum already diagonalized) |
| `n_alpha` / `n_beta` | `int` | UKS seed + \(N_e = n_\alpha + n_\beta` |

`GpwResult::{e_fermi, e_entropy, n_elec_occ}` already exist; UKS path must
populate them when smear is on (currently fixed-occ UKS leaves them 0).

## Scope / non-goals

**In:** multi-k collinear UKS GPW SCF with FD and MP occupations + free energy;
fail-fast validation updates (allow UKS+smear; keep OT+smear banned).  
**Out:** Γ-only UKS smear product path (multi-k only for this slug — metals need
k); forces/stress under smear; adding Fe/Al basis tables; non-orthorhombic cells;
OT+smear.

## Acceptance (see `acceptance.md`)

1. Unit: build **accepts** UKS+FD (was throw); spinor fill recovers \(N_e\).
2. Regression `reg_gpw-uks-smearing`: multi-k UKS + FD — converges; electron
   count; free-energy pin; identity closed-shell UKS+smear ≡ RKS+smear.
3. Dual gate: within_10x or `perf_exempt`.

## References

- Mermin, *Phys. Rev.* **137**, A1441 (1965).
- Methfessel & Paxton, *Phys. Rev. B* **40**, 3616 (1989).
- ATV: `reg_gpw-metallic-smearing`, `reg_pw-uks-fe-smearing`, `reg_gpw-kpoints-uks`.
- CP2K Quickstep `LSD` + `&SMEAR`.

# Spec — `ot-dftu-smearing`

**Loop:** CP2K method | **method:** gpw | **slug:** `ot-dftu-smearing`  
**Branch:** `wf/ot-dftu-smearing` | **Parent:** `dev`

## Capability gap

Multi-k closed-shell **RKS** already productizes:

| Combo | Status |
| ----- | ------ |
| OT + metallic FD/MP smear | ✓ `reg_ot-smearing` |
| OT + Dudarev DFT+U | ✓ `reg_ot-dftu-kpoints` |
| DIIS + smear + DFT+U | ✓ (DIIS multi-k path) |
| **OT + DFT+U + smear** | ✗ throws at multi-k OT entry |

Queue residual (wave 9): multi-k OT path throws when `smearing != none` **and**
`plus_u` non-empty.

## Product target

Public `cpu::gpw::GpwDft` members only (`use_ot`, `plus_u`, `smearing`,
`degauss`, `extra_bands`, `k_mesh` / `k_points`). No `*Config`, no factories.

- **In scope:** multi-k closed-shell RKS OT + Dudarev PLUS_U + FD (or MP)
  metallic occupations → Mermin free energy \(F = E_{\mathrm{KS}} - TS\).
- **Out of scope:** multi-k UKS OT + smear (`ot-uks-smearing`); multi-k UKS
  OT + U (`ot-uks-dftu`); forces/stress of the triple combo.

## Algorithm (implementation sketch)

Reuse the multi-k OT free-energy PR+ Armijo path (already hosts FD/MP
occupations and projected \(g = (I-SCC^\dagger)(2HCf)\)). When `plus_u` is
non-empty, `build_k_energy` already folds \(V_U\) into \(H_k\) and \(E_U\) into
\(F\). Integer-occupation multi-k OT+U keeps the existing OUTER_SCF +
`ot_minimize_fixed_k` branch; **U + smear** takes the smear PR+ path so
fractional occupations stay consistent with DIIS+U+smear.

Convergence: free-energy plateau + soft gnorm (same class as `ot-smearing`);
optional energy-only plateau after a few outer steps when \(U\) softens \(g\).

## Offline oracle

CP2K 2026.1 has **no** OT+kpoints (`OT not possible with kpoint calculations`).
Fair twin = multi-k DIIS + FERMI_DIRAC + PLUS_U (same geometry / U / smear).
Absolute CP2K \(F\) is diagnostic; decisive gate is **ATV OT ≡ ATV DIIS** free
energy + non-trivial \(E_U\) + \(N_e\).

System: H₂ toy (8 Bohr cubic, DZV-GTH-PADE, PADE, CUTOFF 100, MP 2×2×2 or
explicit κ, H \(L=0\) \(U-J=4\) eV, FD `degauss=0.02` Ry, `extra_bands≥1`).

## Acceptance

See `acceptance.md`.

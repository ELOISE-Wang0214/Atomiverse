# Spec: ot-dftu-kpoints — multi-k OT + Dudarev DFT+U

**Slug:** `ot-dftu-kpoints`  
**Domain:** gpw  
**Branch:** `wf/ot-dftu-kpoints`  
**Oracle:** CP2K offline (`tools/cp2k_oracle`)  
**Deps:** `ot-dftu` ✓ (Γ-UKS OT+U), `ot-kpoints` ✓ (multi-k RKS OT),
`gpw-dftu-kpoints` ✓ (multi-k DIIS+U)

## Goal

Productize **closed-shell multi-k Orbital Transformation (OT)** SCF with
**Dudarev PLUS_U** on the same complex k-path. Before this slug, `build()`
rejected `use_ot && !plus_u.empty()` on the k-path, and the multi-k OT
runtime threw even though `build_k_energy` already applied \(V_U\).

## Public API (no `*Config`, no factories)

```cpp
cpu::gpw::GpwDft m;
m.use_ot = true;
m.k_points = {{0.25, 0.0, 0.0}};
m.k_point_weights = {1.0};
m.n_occ = 1;
m.plus_u.push_back({molrs::Element::H, /*l=*/0, /*u_ev=*/4.0});
m.build(atoms);
const auto& r = m.run(atoms); // r.e_total, r.e_dftu [Ha]
```

Members only. Multi-k **UKS** OT remains out of scope (`ot-uks-kpoints`).

## Implementation

1. Remove k-path `build()` reject for `use_ot && plus_u`.
2. Remove multi-k OT runtime throw for `do_dftu`.
3. Soften multi-k OT gradient gate when PLUS_U is on (match Γ OT+U).
4. Shared `build_k_energy` already injects \(E_U / V_U\) (DIIS + OT).

## Acceptance

See `acceptance.md`. Dual gate: identity OT≡DIIS with U + absolute pin;
`perf_exempt` if no fair CP2K wall twin (same pattern as `ot-kpoints` /
`gpw-dftu-kpoints`).

## Out of scope

- Multi-k UKS OT (+U)
- OT + metallic smearing
- Multi-k OT+U analytic forces/stress (reuse DIIS stress paths later)

# Spec — gpw-dftu-kpoints

## Goal

Enable **Monkhorst–Pack multi-k Dudarev DFT+U** on GPW: `k_mesh` / explicit
`k_points` **with** non-empty `plus_u` (CP2K `PLUS_U` + `&KPOINTS`).

Today `GpwDft::build()` rejects that combination. Γ-UKS DFT+U energy
(`reg_gpw-dftu-nio`) and Γ DFT+U stress (`reg_gpw-dftu-stress`) are footholds —
not multi-k capability.

## Scope (this slug)

- **In:** multi-k **collinear UKS** SCF energy with Dudarev PLUS_U (primary);
  multi-k **RKS** PLUS_U if low-cost (closed-shell spin factor 2 × P/2).
- **In:** BZ-averaged Mulliken occupation per spin
  \(n^\sigma=\sum_k w_k\,\tfrac12[(P_k^\sigma S_k)+(S_k P_k^\sigma)]\)
  (Hermitian → real subspace block); energy
  \(E_U^\sigma=\tfrac{U}{2}[\mathrm{Tr}\,n-\mathrm{Tr}(nn)]\);
  potential \(V^\sigma=U(\tfrac12 I-n)\) added to every \(H_k^\sigma\).
- **In:** unit tests (helper + Γ-limit identity + build no longer rejects);
  regression identity suite ≤ 60 s; offline CP2K pin in REFERENCE.
- **Out:** multi-k DFT+U **forces / stress** (keep existing throw until a
  follow-on slug); OT+U; LOWDIN occupation (MULLIKEN only, match Γ path);
  metallic smearing; solid AFM NiO product wall if >60 s (document offline).

## Public API

No new types. Existing members:

```cpp
cpu::gpw::GpwDft dft;
dft.k_mesh[0]=dft.k_mesh[1]=dft.k_mesh[2]=2;
dft.n_alpha=…; dft.n_beta=…;
dft.plus_u.push_back({molrs::Element::Ni, /*l=*/2, /*u_ev=*/6.0});
// no PlusUConfig, no make_*
```

## Acceptance

1. `build()` accepts k-path + non-empty `plus_u` (other k-path rejects unchanged).
2. Explicit k = Γ (complex k-path) with PLUS_U matches real Γ-UKS PLUS_U within
   **1e-6 Ha** on a cheap H₂ UKS toy.
3. Multi-k (MP 2×2×2 TR) + PLUS_U converges; \(|E_U|\) non-trivial when U>0;
   \(|E(U)-E(U=0)|\) ≥ 1e-4 Ha on the same toy.
4. Offline CP2K twin cited in REFERENCE (version, input path, E_tot [Ha]).
5. Dual gate green: correctness + perf (`T_atv/T_oracle ≤ 10` or
   `perf_exempt` with reason).
6. Scoped ctest: unit tests for DFT+U + `run_reg_gpw-dftu-kpoints`.

## Science notes

- Dudarev / Anisimov DFT+U; CP2K `PLUS_U_METHOD MULLIKEN`, `U_MINUS_J` [eV].
- Occupation is **k-summed once** (not a different V per k) so \(E_U[n]\) is
  the energy functional and \(V=\partial E_U/\partial n\) is k-independent.
- MPI k-RR: Allreduce occupation contributions from local k before V/E_U.

## Non-goals

PySCF post-HF, MetalWalls MD, CUDA ports, multi-k DFT+U stress/forces.

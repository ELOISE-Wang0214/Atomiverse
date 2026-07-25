# Spec — periodic-rpa-finite-q

| Field | Value |
| ----- | ----- |
| Slug | `periodic-rpa-finite-q` |
| Method | Finite-q inter-k dRPA dielectric (GPW) |
| Branch | `wf/periodic-rpa-finite-q` |
| Oracle | CP2K offline twin (`tools/cp2k_oracle/inputs/periodic-rpa-finite-q.inp`) |
| Regression | `regressions/gpw/periodic-rpa-finite-q/` |
| Dep | `periodic-rpa-kpoints` ✓ (multi-k **k-diagonal** dRPA product) |

## Problem

Multi-k `GpwDft::direct_rpa` currently builds the Coulomb-only Casida spectrum
**k-diagonally**: same-\(k\) (vertical, \(q=0\)) pair densities only, BZ-averaged
as \(\sum_k w_k\tfrac12(\sum_s\Omega_{k,s}-\mathrm{Tr}\,A_k)\). Full
**finite-\(q\) inter-\(k\) dielectric coupling** — the product multi-k RPA
dielectric — is a declared residual in `mp2_kernels.h` /
`reg_periodic-rpa-kpoints` REFERENCE.

## Science (product)

Closed-shell **direct RPA** correlation on a multi-k GPW KS reference. For each
crystal-momentum transfer \(q\) sampled on the same closed MP mesh as \(k\):

1. Transitions \(\alpha=(k,i,a)\) with virtual at \(k+q\) (momentum
   conservation; wrap into the sampled set).
2. Pair density \(\rho_{ia}^{k,q}=\psi_{ik}^*\psi_{a,k+q}\); Coulomb potential
   via RI-GPW Poisson with transfer \(q\) (`solve_poisson_q` / `poisson_complex`
   with \(q_{\mathrm{frac}}\); \(G=0\) dropped only at \(q=0\)).
3. Coulomb-only Casida (\(\alpha_{\mathrm{hf}}=0\)):
   \(A_{\alpha\beta}=\delta_{\alpha\beta}\Delta\varepsilon_\alpha+2V_{\alpha\beta}\),
   \(B=2V\), Hermitian reduction
   \((A-B)^{1/2}(A+B)(A-B)^{1/2}\) with \(A-B=\mathrm{diag}(\Delta\varepsilon)\).
4. Furche plasmon formula, \(q\)-averaged:
   \[
   E_c=\sum_q w_q\,\tfrac12\Bigl(\sum_s\Omega_{q,s}-\mathrm{Tr}\,A_q\Bigr).
   \]

\(N_k=1\) recovers the real Γ kernel (unit gate). Total energy remains
**dRPA@KS** \(E_{\mathrm{SCF}}+E_c\) (no EXX — periodic HF blocked).

**Not this slug:** frequency-integrated minimax RI-RPA as accuracy twin;
RPAX/SOSEX; UKS; multi-k G0W0 (`periodic-gw-kpoints`); CUDA.

## Public API

No new public type. Members on concrete `GpwDft` unchanged:

```cpp
m.k_mesh[0]=2; m.k_mesh[1]=1; m.k_mesh[2]=1;
m.k_time_reversal = false;
m.mp2_ri = true;
m.mp2_n_virt = 2;
// build(); run(atoms);
const RpaResult& r = m.direct_rpa(atoms);  // q-resolved multi-k product
```

`RpaResult.omega_ha` concatenates all \(q\)-block roots then sorts (length
\(N_q\cdot N_k\cdot n_{\mathrm{occ}}\cdot n_{\mathrm{virt}}\) with \(N_q=N_k\)).

## Implementation plan

1. Replace k-diagonal body of `internal::direct_rpa_rhf_energy_kpoints` with
   full \(q\)-resolved inter-\(k\) dielectric (private `mp2_kernels`).
2. Pass \(q_{\mathrm{frac}}\) into `poisson_complex` for finite-\(q\) pair
   densities (infra already in `poisson_fft` / MP2 helpers).
3. Docs: `rpa.h`, `gpw_solver.h` multi-k comments; unit test \(N_k=1\) recovery
   + optional inter-\(k\) shape gate.
4. Regression + oracle input; update `periodic-rpa-kpoints` pin/omega count if
   the shared multi-k product path changes numbers (same science supersession).

## System (e2e)

Same deck as `periodic-rpa-kpoints`: diamond Si **bct 4-atom**, \(a=5.43\) Å,
DZVP-MOLOPT-GTH + GTH-PBE-q4, PBE, 80 Ry, MP **2×1×1**,
`k_time_reversal=false`, `n_virt=2`, ≤ 60 s.

## Acceptance

See `acceptance.md`. Dual gate: correctness pin in REFERENCE.md; perf_exempt
when no CP2K binary / no fair algorithm twin (CP2K RI-RPA ≠ plasmon on PW-RI).

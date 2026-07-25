# Spec — periodic-gw-finite-q

| Field | Value |
| ----- | ----- |
| Slug | `periodic-gw-finite-q` |
| Method | Finite-q inter-k \(W\) spectral G0W0 (GPW) |
| Branch | `wf/periodic-gw-finite-q` |
| Oracle | CP2K offline twin (`tools/cp2k_oracle/inputs/periodic-gw-finite-q.inp`) |
| Regression | `regressions/gpw/periodic-gw-finite-q/` |
| Dep | `periodic-gw-kpoints` ✓ (multi-k **k-diagonal** spectral G0W0); `periodic-rpa-finite-q` ✓ (finite-q dRPA dielectric) |

## Problem

Multi-k `GpwDft::g0w0` currently builds **k-diagonal** spectral screening:
same-\(k\) (vertical, \(q=0\)) dRPA poles and transition densities only. Full
**finite-\(q\) inter-\(k\) \(W\)** — the product multi-k spectral dielectric —
is a declared residual in `g0w0.h` / `mp2_kernels.h` /
`reg_periodic-gw-kpoints` after `periodic-rpa-finite-q` closed the energy-side
inter-\(k\) Casida path.

## Science (product)

Closed-shell **spectral G0W0@KS** on a multi-k GPW RKS reference. Screening
reuses the **finite-\(q\) inter-\(k\)** Coulomb-only Casida spectrum of
`direct_rpa_rhf_energy_kpoints`:

1. Unique transfers \(q=\mathrm{wrap}(k_a-k_i)\) over the sampled mesh.
2. For each \(q\): transitions \(\alpha=(k,i,a)\) with virtual at \(k+q\);
   pair density \(\rho_{ia}^{k,q}=\psi_{ik}^*\psi_{a,k+q}\); RI-GPW Poisson
   (same Bloch-product convention as multi-k RPA/MP2); Hermitian reduction
   with \(X+Y\) amplitudes (not energy-only).
3. Spectral transition potentials
   \(W_s^q(r)=\sum_\alpha\sqrt{2}\,(X+Y)_\alpha^s\,V[\rho_\alpha](r)\).
4. For each target \(k\) and orbital \(p\), intermediate states \(m\) at
   \(k_m=k-q\) (momentum conservation for \(W\) carrying \(q\)):
   \[
   t_{s,m,p}^{k,q}
   =\int W_s^q(r)\,\psi_m^{k-q*}(r)\,\psi_p^k(r)\,d^3r,
   \]
   and
   \[
   \mathrm{Re}\,\Sigma_{c,p}^k(\omega)
   =\sum_q\frac{1}{N_k}\sum_s\sum_m
   |t_{s,m,p}^{k,q}|^2\,
   \frac{\omega\pm\Omega_{q,s}-\varepsilon_m^{k-q}}
        {(\omega\pm\Omega_{q,s}-\varepsilon_m^{k-q})^2+\eta^2}
   \]
   (\(+\) for occupied \(m\), \(-\) for virtual \(m\); \(\eta\) spectral
   broadening as in the Γ / k-diagonal product).
5. Exact-exchange diagonals: truncated (MT) Coulomb per \(k\) (unchanged
   G=0-safe foothold). Local \(v_{\mathrm{xc}}\): same pure GGA/LDA as SCF on
   the BZ-weighted density.
6. Newton QP orbit-by-orbit, per \(k\); report min-direct QP gap.

\(N_k=1\) recovers the real Γ spectral kernel (unit gate).

**Not this slug:** frequency-integrated RI-GW (CP2K default); UKS/UGW;
periodic HF-X (EXX total energy); minimax \(W(\mathrm{i}\omega)\); CUDA.

## Public API

No new public type. Members on concrete `GpwDft` unchanged:

```cpp
m.k_mesh[0]=2; m.k_mesh[1]=1; m.k_mesh[2]=1;
m.k_time_reversal = false;
m.mp2_ri = true;
m.mp2_n_virt = 2;
// build(); run(atoms);
const G0W0Result& gw = m.g0w0(atoms);  // finite-q multi-k product
// gw.n_kpts == 2; gw.gap_ev = min_k direct QP gap [eV]
```

## Implementation plan

1. Upgrade `internal::spectral_g0w0_rhf_kpoints` to take `k_frac` and build
   finite-\(q\) inter-\(k\) \(W\) (reuse RPA q-mesh / Poisson block; keep
   \(X+Y\)).
2. Wire `GpwDft::g0w0` multi-k path to pass `result_.k_frac`.
3. Docs: `g0w0.h`, `gpw_solver.h`, `mp2_kernels.h` — product is finite-q;
   k-diagonal residual closed.
4. Unit: \(N_k=1\) multi-k ≡ real spectral QP (`test_cpu_gpw_g0w0`).
5. Regression + oracle input; update `periodic-gw-kpoints` REFERENCE to note
   product supersession (same public call; numbers may retarget).

## System (e2e)

Same deck as `periodic-gw-kpoints` / `periodic-rpa-finite-q`: diamond Si
**bct 4-atom**, \(a=5.43\) Å, DZVP-MOLOPT-GTH + GTH-PBE-q4, PBE, 80 Ry,
MP **2×1×1**, `k_time_reversal=false`, `n_virt=2`, ≤ 60 s.

## Acceptance

See `acceptance.md`. Dual gate: correctness pin in REFERENCE.md; perf_exempt
when no CP2K binary / no fair algorithm twin (CP2K RI-G0W0 ≠ spectral + MT).

# Reference — multi-k RI-GPW dRPA on diamond Si (bct 4-atom, MP 2×1×1)

## Method

Closed-shell **k-diagonal direct RPA (dRPA)** correlation on top of a multi-k
**GPW PBE** SCF. For each sampled \(k\), the Coulomb-only Casida spectrum is
built from **same-\(k\)** (vertical, \(q=0\)) occupied–virtual pair densities
via RI-GPW Poisson (complex Bloch products; linear real Poisson on Re/Im),
reduced by the Furche plasmon formula, and Brillouin-zone averaged:

\[
E_c=\sum_k w_k\,\tfrac12\Bigl(\sum_s\Omega_{k,s}-\mathrm{Tr}\,A_k\Bigr).
\]

Total energy is **dRPA@KS** \(E_{\mathrm{SCF}}+E_c\) (no EXX — periodic HF
blocked). Full inter-\(k\) / finite-\(q\) dielectric coupling is a follow-on.

Public path:

```cpp
cpu::gpw::GpwDft m;
m.k_mesh[0] = 2; m.k_mesh[1] = 1; m.k_mesh[2] = 1;
m.k_time_reversal = false;
m.mp2_ri = true;
m.mp2_n_virt = 2;
// build(); run(atoms);
const cpu::gpw::RpaResult& r = m.direct_rpa(atoms);  // multi-k, no throw
```

Kernel Γ recovery (\(N_k=1\), real-valued complex MOs ≡ real plasmon) is
unit-tested in `test_cpu_gpw_rpa` (no full SCF).

## System

| Item | Value |
|------|--------|
| Crystal | diamond Si, **bct 4-atom** cell (same as `gpw-kpoints-si`) |
| Cell | a_t = a/√2, c = a; a = **5.43 Å** |
| Positions (frac) | (0,0,0), (½,½,½), (½,0,¼), (0,½,¾) |
| Basis / PP | **DZVP-MOLOPT-GTH** / **GTH-PBE-q4** |
| XC | **PBE**, cutoff **80 Ry**, eps_default 1e-7 |
| Valence | 16 e⁻ → n_occ = 8 |
| k | **MP 2×1×1**, `k_time_reversal=false` → 2 k, w = ½ |
| SCF | DIIS 6, tol 1e-7, mix 0.3, **no level_shift** |
| RPA | RI-GPW k-diagonal dRPA, **n_virt = 2** (e2e ≤ 60 s) |

## Offline oracle status

**No CP2K binary / module on the worker host at pin time** (same status as
`periodic-mp2-kpoints` / `periodic-rpa-si`). Absolute CP2K digits are **not
invented**. Input twin:
`tools/cp2k_oracle/inputs/periodic-rpa-kpoints.inp`.

CP2K's default path is frequency-integrated RI-RPA with Gaussian `RI_AUX`
and EXX total-energy assembly — not the exact k-diagonal plasmon formula on
PW-RI \((ia|jb)\). Absolute CP2K \(E_c\) is therefore not a hard twin.

### Binding gates

1. **Hard multi-k:** SCF converged; `C_k.size()==n_kpts`; `e_corr` finite and
   **strictly negative**; all \(\Omega>0\); e_total ≡ e_scf+e_corr; n_occ/n_virt;
   `ri` and `n_kpts` flags.
2. **Kernel Γ recovery (unit):** multi-k kernel at \(N_k=1\) ≡ real
   `direct_rpa_rhf_energy`.
3. **First-green absolute band** (setup integrity) — freeze ATV digits below;
   replace with fair CP2K when available.

### Atomiverse numbers (first green, 2026-07-25)

| Quantity | Value |
|----------|--------|
| E_SCF (2k, PBE, 80 Ry) | **−15.5751782413** Ha (7 iters; stable) |
| **E_corr (2k, n_virt=2)** | **≈ −0.0163** Ha (OMP-dependent) |
| E_corr (OMP=1, deterministic) | −0.0157519268 Ha |
| E_corr (OMP=4, sample) | −0.0154 … −0.0173 Ha |
| Wall (SCF + RPA) | **~3 s** |

Regression absolute tol on E_corr: **5e-3 Ha** about −0.0163 (setup band
covering OMP orbital scatter on the multi-k SCF; not a CP2K digit twin).

## Perf

Fair wall vs CP2K when binary exists: target
\(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\). Without CP2K binary / without
fair algorithm twin (CP2K RI-RPA ≠ k-diagonal plasmon): **perf_exempt**.
Never launch CP2K from ctest.

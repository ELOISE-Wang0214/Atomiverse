# Spec: dfpt-spin-smearing

**Slug:** `dfpt-spin-smearing`
**Method domain:** response (PW DFPT)
**Branch:** `wf/dfpt-spin-smearing`
**Oracle:** offline CP2K (preferred) / same-Hamiltonian free-energy FD pin
**Deps:** `dfpt-spin-si` ✓, `dfpt-metallic-si` ✓ (separately productized)

## Intent

Productize **Γ-point DFPT phonons** for `cpu::pw::PW` with **collinear spin**
and **metallic occupation smearing** together:

- `method.spin_polarized = true`
- `method.smearing ∈ {fermi_dirac, methfessel_paxton, gaussian}`
- `degauss > 0` [Ry], `extra_bands ≥ 0`

Shared μ over spin-orbitals (spin_degen = 1), \(f_{\sigma nk}\in[0,1]\),
spin XC kernel unchanged from `dfpt-spin-si`, Sternheimer weights and
optional δμ charge neutrality as in `dfpt-metallic-si` (de Gironcoli).

## Out of scope

- Finite-q + smear (`dfpt-finite-q-smearing`)
- Hybrid XC / exact exchange response
- CUDA DFPT

## Public API

```cpp
cpu::pw::PW method;
method.spin_polarized = true;
method.multiplicity = 1;
method.smearing = cpu::pw::Smearing::fermi_dirac;
method.degauss = 0.02; // Ry
method.extra_bands = 0;
cpu::pw::DfptPhonon dfpt;
auto r = dfpt.compute(method, frame);
```

No public `*Config`. Frame in, `DfptPhononResult` out.

## Acceptance

1. Unit: `ThrowsOnSpinWithSmearing` no longer throws the spin+smear gate;
   empty/missing UPF fails at I/O (parity with other allowed paths).
2. Regression `regressions/response/dfpt-spin-smearing`:
   - acoustic |ω| ≤ 30 cm⁻¹ (no ASR)
   - optical Γ triplet in [500, 650] cm⁻¹, spread ≤ 10
   - C[(0x),(0x)] vs offline free-energy FD pin |Δ| ≤ 0.30 eV/Å²
   - ASR row-sum ≤ 0.05 eV/Å²; asymmetry ≤ 1e-5
   - wall ≤ 60 s
3. Dual gate green; `perf_exempt` if no CP2K binary on host (document note).
4. Docs: header comments + `tools/cp2k_oracle/inputs/dfpt-spin-smearing.note`

## Twin settings (minimal textbook)

Diamond Si, a=5.43 Å, 2-atom fcc, `Si_lda_nc.upf`, LSDA, ecut 10 Ry,
2×2×2 MP, FD degauss 0.02 Ry, multiplicity=1, extra_bands=0 (gapped fill).

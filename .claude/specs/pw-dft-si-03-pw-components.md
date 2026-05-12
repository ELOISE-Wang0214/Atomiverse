---
title: PW SCF 组件
slug: pw-dft-si-03-pw-components
status: approved
created: 2026-05-12
---

# PW SCF 组件

## Summary

实现六个继承 `cpu::ComponentBase` 的平面波 SCF 组件，挂载 `scf::Phase`，通过 `PWWorkspace` 共享数据：

| 组件 | Phase | Order | 功能 |
|------|-------|-------|------|
| PW_Poisson | build_hamiltonian | -100 | V_H(G) = 4π ρ(G)/G² |
| PW_LDA | build_hamiltonian | -50 | ρ(G)→FFT→evaluate_xc→FFT→V_xc(G) |
| PW_Hamiltonian | build_hamiltonian | 0 | T(G)=½|k+G|², H·ψ |
| PW_Veff | build_hamiltonian | +50 | V_eff = V_loc + V_H + V_xc |
| PW_Davidson | diagonalize | 0 | 矩阵自由迭代对角化 |
| PW_Density | build_density | 0 | ψ→FFT→ρ(r)→FFT→ρ(G) + Pulay混合 |
| PW_Convergence | convergence_check | 0 | ‖ρ_out−ρ_in‖ < tol |

复用现有 `Driver<scf::Phase, cpu::device_tag>` 编排。组件不修改 `cpu::Frame`/`cpu::Context`。

## Domain basis

- KS方程平面波形式: H(G,G') = ½|k+G|² + V_eff(G-G') (Kohn & Sham, PR 140, A1133, 1965)
- Davidson: Davidson, JCP 17, 87 (1975); Payne et al., RMP 64, 1045 (1992)
- Pulay: Pulay, CPL 73, 393 (1980)
- VWN5: Vosko, Wilk, Nusair, CJP 58, 1200 (1980)
- Slater exchange: Dirac, PCPS 26, 376 (1930)
- KB pseudopotential: Kleinman & Bylander, PRL 48, 1425 (1982)

## Design

### PWWorkspace

```cpp
struct PWWorkspace {
    xt::xarray<std::complex<Real>> rho_G, rho_G_old;
    xt::xarray<std::complex<Real>> V_H_G, V_xc_G, V_eff_G;
    xt::xarray<Real> rho_r, V_xc_r;
    // Pulay mixer buffers
    // Davidson workspace per k-point
};
```

`BuildContext` 增加 `void* pw_workspace = nullptr`。首次构建者分配，后续引用。

### 组件

**PW_Poisson**: V_H(G) = 4π ρ(G)/|G|², G=0→0.

**PW_LDA**: (1) FFT ρ(G)→ρ(r); (2) evaluate_xc(LSDA, ρ, 0)→V_xc(r); (3) FFT V_xc(r)→V_xc(G).

**PW_Hamiltonian**: build()预计算T(G)数组; launch()提供H·ψ = T·ψ + V_eff·ψ.

**PW_Veff**: V_eff(G) = V_loc(G) + V_H(G) + V_xc(G).

**PW_Davidson**: 逐k点迭代。Rayleigh-Ritz投影→残差修正。N_occ个本征对，残差<1e-8.

**PW_Density**: ψ→FFT→ρ(r)=Σ|ψ|²→FFT→ρ(G). Pulay: ρ_mix = ρ_in + Σα_i(ρ_out-ρ_in)_i.

**PW_Convergence**: diff = ‖ρ_out-ρ_in‖. tol=1e-8，预热2轮.

## Files

- `include/atomiverse/cpu/pw/pw_workspace.h` (new)
- `include/atomiverse/cpu/pw/pw_poisson.h` + `src/cpu/pw/pw_poisson.cpp` (new)
- `include/atomiverse/cpu/pw/pw_lda.h` + `src/cpu/pw/pw_lda.cpp` (new)
- `include/atomiverse/cpu/pw/pw_hamiltonian.h` + `src/cpu/pw/pw_hamiltonian.cpp` (new)
- `include/atomiverse/cpu/pw/pw_veff.h` + `src/cpu/pw/pw_veff.cpp` (new)
- `include/atomiverse/cpu/pw/pw_davidson.h` + `src/cpu/pw/pw_davidson.cpp` (new)
- `include/atomiverse/cpu/pw/pw_density.h` + `src/cpu/pw/pw_density.cpp` (new)
- `include/atomiverse/cpu/pw/pw_convergence.h` + `src/cpu/pw/pw_convergence.cpp` (new)
- `tests/cpu/pw/test_pw_components.cpp` (new)

## Tasks

- [ ] Write failing tests for all six components + PWWorkspace sharing
- [ ] Implement PWWorkspace with PulayMixer + Davidson buffers
- [ ] Implement PW_Poisson — V_H(G)=4πρ/G², G=0 guard
- [ ] Implement PW_LDA — FFT + evaluate_xc + FFT back
- [ ] Implement PW_Hamiltonian — T(G) precompute + H·ψ
- [ ] Implement PW_Veff — V_loc + V_H + V_xc assembly
- [ ] Implement PW_Davidson — matrix-free iterative eigensolver
- [ ] Implement PW_Density — FFT density + Pulay mixing
- [ ] Implement PW_Convergence — ‖ρ_out−ρ_in‖ check with warmup
- [ ] Register sources in CMakeLists.txt; run full check + test suite

## Testing strategy

- PW_Poisson: Gaussian ρ→ analytic erf solution, rel error<1e-4
- PW_LDA: uniform ρ→ VWN5 reference, abs diff<1e-10
- PW_Hamiltonian: k=0, single PW→ eigenvalue=½|G|² within 1e-12
- PW_Davidson: free-electron H, N_occ=8, residual<1e-8 in <50 iter
- PW_Density: FFT round-trip identity; Pulay ≥2x better than simple mixing
- PW_Convergence: tol=1e-8 boundary correct; warmup=2 suppresses early
- Integration: 2-atom Si, E_cut=10Ry, Gamma-only→ converge in <20 SCF steps

## Out of scope

- GPU/CUDA
- 非局域赝势 V_nl（后续组件）
- 自旋极化 LSDA
- GGA/PBE
- Anderson/DIIS/Kerker 混合器

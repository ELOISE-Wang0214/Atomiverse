---
title: PW DFT 基础设施
slug: pw-dft-si-02-pw-infra
status: approved
created: 2026-05-12
---

# PW DFT 基础设施

## Summary

实现平面波 DFT 的 CPU 端基础设施工具集：G 矢量枚举、FFTW3 3D FFT RAII 封装、UPF v1 赝势解析器、Monkhorst-Pack k 点网格、Si 金刚石晶体构造器、3D Ewald 离子能量求和。全部位于 `cpu/pw/internal/`，纯自由函数/工具类，不继承 ComponentBase。

## Domain basis

- **G 矢量**: G = h·b1 + k·b2 + l·b3, 满足 |k+G|²/2 ≤ E_cut。Martin, *Electronic Structure* (Cambridge 2004) §12.3。
- **FFT**: Frigo & Johnson, Proc. IEEE 93, 216 (2005). R2C/C2R, half-complex存储。
- **UPF v1**: Hamann, PRB 88, 085117 (2013); SG15: Schlipf & Gygi, CPC 196, 36 (2015).
- **k 点**: Monkhorst & Pack, PRB 13, 5188 (1976). 均匀网格等权重。
- **Ewald**: Ewald, Ann. Phys. 64, 253 (1921). E = E_recip + E_real + E_self。G_cut = 4*E_cut。

内部单位: Hartree 原子单位 (能量 Hartree, 长度 Bohr)，输出转 eV/Angstrom。

## Design

### GVectorList (g_vector.h/.cpp)

构造: `GVectorList(Box cell, Real ecut, Vec3 k)`. 扫描 (h,k,l) 包围盒，筛选 |k+G|²/2 ≤ E_cut，按 |G| 排序。
- `count()`, `operator[](i)`, `kinetic_energies()`, `fft_dims()`

### FFTPlan3D (fft.h/.cpp)

RAII 封装 FFTW3 plan。构造(Nx,Ny,Nz, direction)。禁止拷贝。析构释放 plan。
- `r2c(Real* in, fftw_complex* out)`, `c2r(fftw_complex* in, Real* out)`

### UPFReader (upf_reader.h/.cpp)

纯字符串解析 UPF v1 XML 风格格式。提取 Z_val, l_max, mesh, V_loc, rho_atom, KB β projectors, D coefficients。

### KPointGrid (kpoint.h/.cpp)

`generate_mp_mesh(nk1,nk2,nk3, shift) → KPointGrid {kpoints, weights}`. 分数坐标。

### Crystal (crystal.h/.cpp)

`build_si_diamond(Real a=5.43)` → `cpu::Frame`. FCC晶胞, 2原子基元 (0,0,0) + (1/4,1/4,1/4).

### EwaldSummation (ewald.h/.cpp)

构造: `EwaldSummation(positions, charges, cell, ecut)`. α自动选取。提供 `total_energy()`, `real_energy()`, `recip_energy()`, `self_energy()` (eV).

## Files

- `include/atomiverse/cpu/pw/internal/g_vector.h` + `src/cpu/pw/internal/g_vector.cpp` (new)
- `include/atomiverse/cpu/pw/internal/fft.h` + `src/cpu/pw/internal/fft.cpp` (new)
- `include/atomiverse/cpu/pw/internal/upf_reader.h` + `src/cpu/pw/internal/upf_reader.cpp` (new)
- `include/atomiverse/cpu/pw/internal/kpoint.h` + `src/cpu/pw/internal/kpoint.cpp` (new)
- `include/atomiverse/cpu/pw/internal/crystal.h` + `src/cpu/pw/internal/crystal.cpp` (new)
- `include/atomiverse/cpu/pw/internal/ewald.h` + `src/cpu/pw/internal/ewald.cpp` (new)
- `tests/cpu/pw/test_g_vector.cpp`, `test_fft.cpp`, `test_upf.cpp`, `test_kpoint.cpp`, `test_crystal.cpp`, `test_ewald.cpp` (new)

## Tasks

- [x] Write failing tests for G-vector enumeration — count in [250,300] for Si at 30Ry
- [x] Implement GVectorList
- [x] Write failing tests for FFTW3 — roundtrip error < 1e-12
- [x] Implement FFTPlan3D RAII wrapper
- [x] Write failing tests for UPF parser — Z_val==4, l_max>=1, arrays non-empty
- [x] Implement UPFReader
- [x] Write failing tests for k-point mesh — weights sum to 1, 4x4x4=64 points
- [x] Implement KPointGrid generator
- [x] Write failing tests for Si crystal — 2 atoms, correct positions
- [x] Implement build_si_diamond()
- [x] Write failing tests for Ewald — energy converges on cutoff doubling
- [x] Implement EwaldSummation
- [x] Add FFTW3 find_package + PW infra sources to CMake; add test targets
- [x] Run full check + test suite

## Testing strategy

- G矢量: Si FCC a=10.26Bohr, E_cut=30Ry → count∈[250,300]; 所有G满足|G|²/2≤E_cut
- FFT: 12³ grid随机填充→R2C→C2R→max|diff|<1e-12
- UPF: SG15 Si ONCV → Z_val=4, l_max=1, mesh non-empty
- k点: 4×4×4=64点, Σweights=1.0
- 晶体: 2原子, type_id=14, Si-Si距离=a√3/4=2.352Å
- Ewald: 2×Si(+4), E_cut vs 2×E_cut差<1e-6 Hartree

## Out of scope

- GPU/CUDA 版本
- UPF v2 / USPP / PAW
- Ewald 力
- k 点对称性约化

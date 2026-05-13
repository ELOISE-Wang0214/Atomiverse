---
title: PW DFT 驱动器和 Si 总能量基准
slug: pw-dft-si-04-pw-driver
status: approved
created: 2026-05-12
---

# PW DFT 驱动器和 Si 总能量基准

## Summary

实现 `pw::solve_si()` 集成入口：构建 Si 金刚石晶胞 → 加载 SG15 ONCV 赝势 → 生成 G 矢量 → k 点网格 → Ewald → 创建 Driver + 注册 6 个 PW 组件 → SCF 循环 → 返回总能量。CMake `ATV_BUILD_PW` 选项门控。集成测试验证 E_total ≈ -107.8 Ry ± 0.5 Ry。

## Domain basis

Si 金刚石, a=5.43 Å, SG15 ONCV 赝势 (Z_val=4), E_cut=30 Ry, 4×4×4 MP k 点, LDA (VWN5). 参考: Hamann, PRB 88, 085117 (2013); Schlipf & Gygi, CPC 196, 36 (2015). 目标值: -107.8 Ry/cell (-1467 eV / 2 atoms).

## Design

### `pw::solve_si()`

```cpp
namespace ATV::cpu::pw {
double solve_si(const std::string& upf_path,
                double ecut_ry = 30.0,
                int k_mesh = 4,
                double tol = 1e-8,
                int max_iter = 100);
}
```

流程:
1. `build_si_diamond(5.43)` → `cpu::Frame`
2. `UPFReader(upf_path)` → 赝势数据
3. `GVectorList(box, ecut_ry/2, k)` → G 矢量
4. `generate_mp_mesh(4,4,4)` → k 点
5. `EwaldSummation` → E_ewald (常数)
6. `Driver<scf::Phase, cpu::device_tag>` → add 6 组件
7. `build()` → `run(max_iter)` → `destroy()`
8. 收集 E_kin + E_H + E_xc + E_loc + E_nl + E_ewald → return

### CMake 门控

```cmake
option(ATV_BUILD_PW "Build plane-wave DFT module" OFF)
if(ATV_BUILD_PW)
    find_package(FFTW3 REQUIRED)
    target_compile_definitions(atomiverse PUBLIC ATV_HAS_PW)
endif()
```

## Files

- `include/atomiverse/cpu/pw/solve_si.h` (new)
- `src/cpu/pw/solve_si.cpp` (new)
- `CMakeLists.txt` (modify — ATV_BUILD_PW + FFTW3)
- `tests/cpu/pw/test_pw_si.cpp` (new)

## Tasks

- [ ] Write failing test for pw::solve_si() — E_total within ±0.5 Ry of -107.8
- [ ] Add ATV_BUILD_PW option + find_package(FFTW3 REQUIRED) to CMakeLists.txt
- [ ] Register pw sources under ATV_BUILD_PW guard in CMakeLists.txt
- [ ] Register test_pw_si target under ATV_BUILD_PW guard in tests/CMakeLists.txt
- [ ] Declare pw::solve_si() with Doxygen in solve_si.h
- [ ] Implement solve_si(): Si crystal + UPF + G-vector + k-points + Ewald + Driver + 6 components + SCF loop + energy sum
- [ ] Run full check + test suite

## Testing strategy

- E_total: |E + 107.8| < 0.5 Ry (Si, 30 Ry, 4×4×4)
- ATV_BUILD_PW=OFF: 无 PW 符号、无 test_pw_si target
- ATV_BUILD_PW=ON 无 FFTW3: CMake 报错含 "FFTW3"
- test_pw_si 编译链接运行无 crash

## Out of scope

- 力/应力/几何优化
- 能带/态密度输出
- 波函数文件
- k 点并行

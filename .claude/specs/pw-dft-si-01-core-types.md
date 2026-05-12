---
title: PW-DFT 核心类型层
slug: pw-dft-si-01-core-types
status: done
created: 2026-05-12
---

# PW-DFT 核心类型层

## Summary

重写 `Mat3` 为基于 `xt::xtensor_fixed` 的 OOP 类，新建 `Vec3` OOP 值类型和 `Box` 晶胞值类型，扩展 `units.h` 的 PW-DFT 单位常量，扩展 `types.h` 的类型别名。干净重写，不保留旧接口兼容。

## Domain basis

QE 约定单位常量：
- `kAngstrom2Bohr` = 1.8897259886
- `kRy2Hartree` = 0.5 (exact)
- `kRy2eV` = 13.605693

## Design

### Mat3（重写）

`class Mat3`，内部 `xt::xtensor_fixed<Real, xt::xshape<3,3>>`。

- 元素: `operator()(i,j)`
- 算术: `+`, `-`, `*(scalar)`, `*(Mat3)`, `*(Vec3)`, `/(scalar)`
- 线代: `transpose()`, `inverse()`, `det()`, `trace()`, `norm()`
- 工厂: `zeros()`, `identity()`, `diagonal(a,b,c)`

### Vec3（新建）

`class Vec3`，内部 `xt::xtensor_fixed<Real, xt::xshape<3>>`。

- 元素: `operator[](i)`, `x()`, `y()`, `z()`
- 算术: `+`, `-`, `*(scalar)`, `/(scalar)`
- 向量: `dot(v)`, `cross(v)`, `norm()`, `norm_squared()`, `normalized()`
- 工厂: `zeros()`, `unit_x()`, `unit_y()`, `unit_z()`

### Box（新建）

`struct Box`，成员 `Mat3 cell`, `Vec3 origin`, `bool pbc[3]`。

- `volume()` — |det(cell)|
- `recip_cell()` — 2π·(cell⁻¹)^T
- `lengths()` — [|a|,|b|,|c|]
- `to_frac(Vec3)` / `to_cart(Vec3)` — 坐标变换
- `wrap(Vec3)` / `wrap_batch(Real*,int)` — 折叠入胞
- `shortest_vector(Vec3,Vec3)` — 最小镜像
- `static cubic(Real a)`

### units.h 扩展

`ATV` 命名空间添加 `kAngstrom2Bohr`, `kRy2Hartree`, `kRy2eV`。

### types.h 扩展

添加 `using Mat33 = ATV::Mat3;` `using Vec3 = ATV::Vec3;`。

## Files

- `include/atomiverse/core/mat3.h` (rewrite)
- `include/atomiverse/core/vec3.h` (new)
- `include/atomiverse/core/box.h` (new)
- `include/atomiverse/core/units.h` (modify)
- `include/atomiverse/core/types.h` (modify)
- `tests/core/test_mat3.cpp` (new)
- `tests/core/test_vec3.cpp` (new)
- `tests/core/test_box.cpp` (new)
- `tests/core/test_units.cpp` (new)

## Tasks

- [x] Write failing tests for Mat3 — inverse×Mat3==I, det matches cofactor, Mat3×Vec3 matches reference
- [x] Implement Mat3 OOP class with xtensor_fixed + full operator set
- [x] Write failing tests for Vec3 — dot/cross/norm match analytic, normalized()==1
- [x] Implement Vec3 OOP class
- [x] Write failing tests for Box — to_frac/to_cart roundtrip, wrap in [0,L), shortest_vector in [-0.5,0.5], recip×cell==2πI
- [x] Implement Box value type
- [x] Write failing tests for new unit constants
- [x] Add kAngstrom2Bohr/kRy2Hartree/kRy2eV to units.h; add Vec3/Mat33 aliases to types.h
- [x] Fix all existing Mat3 call sites across codebase to use new interface
- [x] Run full check + test suite
- [x] Hygiene — /mol:simplify applied 5 fixes (3 dead includes, 1 magic number, 1 local constant dedup)

## Testing strategy

- Mat3: 5 random non-singular matrices, M×M⁻¹≈I within 1e-8 Frobenius
- Vec3: 5 random (u,v), dot/cross/norm matches analytic; cross(u,u)==0
- Box: cubic(10), recip×cell≈2πI, wrap puts coords in [0,L), shortest_vector fractions in [-0.5,0.5]
- Units: 5.43*kAngstrom2Bohr≈10.26; kRy2Hartree==0.5; kRy2eV==13.605693

## Out of scope

- 旧 Mat3 兼容接口 — 干净重写，不保留任何旧式 API
- SIMD 向量化
- CUDA 端 Mat3/Box

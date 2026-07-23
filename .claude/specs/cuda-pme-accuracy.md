---
spec: cuda-pme-accuracy
status: approved
created: 2026-07-23
chain: cuda-pme-accuracy
phase: single
grilled: false
---

# cuda-pme-accuracy — CUDA PME 倒空间能量偏差 3%

## Summary

`tests/cuda/longrange/test_pme_coul.cu` 的 `CudaPmeCoul.MadelungNaCl` 失败:CUDA PME
算出的 NaCl 静电能比解析 Madelung 参考值**低约 3.0%**,而 **CPU PME 在完全相同的参数下
只差 5e-4**。这是 **既存缺陷**,不是 `cuda-resurrect` 引入的 —— CUDA 后端在 `dev` 上根本
编不过,所以这个测试**从来没跑过**;复活后第一次跑就暴露了。

## 证据(已确认)

```
E        = -276.86183347465663
expected = -285.5524513734664      (= -0.5 * N * M / d * units::coulomb)
ratio    = 0.9696                  → 3.0% 低,容差 3e-3
```

测试参数:`d=2.82 Å`,`n=4`(64 原子,盒子 11.28 Å),`cutoff=1.9d`,
`alpha=3.5/cutoff`,PME `grid_spacing=0.5`,`order=6`,`M=1.747564594633182`。

**已排除的原因**(cuda-resurrect 期间逐一验证):

1. **不是单位改名**。`kCoulomb→units::coulomb` / `kPI→units::pi` 在 `e/e_expected`
   的比值里**完全抵消**(两边都用同一常量)。
2. **不是浮点精度**。用 `-DATV_CUDA_REAL=double` 重新构建后**仍然失败**
   (GPU job 570473)——所以不是 `real_t=float` 的输入精度问题。
3. **不是漏项**。能量项结构完整:倒空间 `src/cuda/longrange/pme.cu:137`
   (`escale * H * |S(k)|²`)、自能 `:154` + `:321`
   (`self_pref = -coulomb*alpha/√π`)。净电荷修正项对中性 NaCl 恒为 0,
   本就不该出现。
4. **不是我们的 Frame/scalars 桥**。同一条 `scalars.energy → kEnergy` 路径下,
   其它能量测试(LJ/bonded/Ewald)全绿。

## 怀疑方向(未验证,需上机数值二分)

- **B 样条归一化 / order-6 权重**:`pme_spread` / `pme_interp`
  (`pme.cu:64` / `:144`)的样条系数与 CPU 实现逐项对比。
- **网格维度**:由 `grid_spacing=0.5` 推出的 `K0/K1/K2` 是否与 CPU 一致
  (取整/取偶/FFT 友好尺寸)。
- **影响函数归一化**:`pme_influence`(`:105`)的 `escale`、`H`、以及
  k=0 与 Nyquist 分量的处理。
- **实空间/倒空间劈分一致性**:`cuda::CoulLong` 的 α 与排除项是否与 KSpacePME 严格配套。

## Tasks

- [ ] 上机二分:分别只跑 `CoulLong`(实空间)与只跑 `KSpacePME`(倒空间+自能),
      与 CPU 同参数逐项对比,定位 3% 落在哪一项。
- [ ] 按定位结果修 `src/cuda/longrange/pme.cu`(样条/网格/影响函数其一)。
- [ ] `CudaPmeCoul.MadelungNaCl` 达到 ≤ 3e-3(目标对齐 CPU 的 5e-4)。
- [ ] 确认 `test_cuda_ewald_coul` 仍绿(Ewald 与 PME 共用 α 劈分约定)。

## Testing strategy

GPU 节点 `sbatch`(配方见 memory `cuda-gpu-node-build-recipe`)。二分用
`ctest -R test_cuda_pme_coul --output-on-failure`;对照值取 CPU 的
`tests/cpu/longrange/test_pme_coul.cpp` 与 `regressions/md/npt_pme_nacl`
(REFERENCE.md 记录同网格/阶数下 CPU 达 5e-4)。

## Out of scope

- 改 `ATV_CUDA_REAL` 默认值(保持 `float`)。
- CPU PME/Ewald 任何改动。
- 放宽测试容差来"修"失败 —— 容差是物理门禁,CPU 同参数能达标,CUDA 也应该。

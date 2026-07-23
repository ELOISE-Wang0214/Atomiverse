---
slug: analytic-tda-forces
status: approved
created: 2026-07-23
track: A-response
grilled: false
---

# analytic-tda-forces

## Summary

为闭壳层 **RHF TDA（CIS）** 提供 **解析核梯度** \(-\nabla_R E_n^{\mathrm{TDA}}\) [eV/Å]，作为 `RHF::tda_forces` 的默认/可选解析路径；现有中心差分实现保留为对照与回归 oracle。使激发态势能面可在不付出 6N 次 SCF+TDA 的代价下驱动 Minimize / 为 NAC 提供解析力差。

## Domain basis

CIS/TDA 激发能对核坐标的梯度（Furche–Ahlrichs / Maurice–Head-Gordon 路线；闭壳层）：

\[
\frac{\partial E_n}{\partial R} =
\frac{\partial E_{\mathrm{HF}}}{\partial R}
+ \mathbf{X}_n^\top \frac{\partial \mathbf{A}}{\partial R}\mathbf{X}_n
+ \text{(orbital-response / Z-vector if needed)}
\]

对纯 CIS（无 DFT 核），常用高效形式用未弛豫密度 + 拉格朗日乘子消除 CPHF 显式求解，或走一次 Z-vector。参考：Furche & Ahlrichs, *J. Chem. Phys.* **117**, 7433 (2002), DOI 10.1063/1.1508368；Maurice & Head-Gordon, *Mol. Phys.* **96**, 1533 (1999)。单位：内部 a.u.；公开力 eV/Å（`core/units.h`）。

## Design

- **API**：扩展 `RHF::tda_forces(frame, state, …)` — 默认解析；可保留 `h` 参数仅用于 FD 对照重载，或加枚举/bool `analytic=true`。禁止新 factory；成员方法 only。
- **私有工作区**：TDA 向量 \(\mathbf{X}\)、过渡密度、（若需）Z-vector 求解器放在 `cpu/scf/internal/`，不出现 `IntegralPack` 公共类型。
- **正确性门**：与现有 FD `tda_forces`（`reg_h2_tda_force`）在 H2/甲醛小基上分量一致（相对或绝对容差文档化）。
- **错误**：未 `run()` / 无效 state / 无虚轨道 → throw，无静默零力。

### Reuse decision

- `reuse` `RHF::tda` / `tda_forces` 签名落点（`rhf.h`）— 解析实现替换或并列 FD。
- `reuse` 基态 `RHF::forces` / 梯度基础设施（`gradient` / 积分导数路径）。
- `reuse` `reg_h2_tda_force` FD 金标作交叉验证。
- `pattern` `mp2.h` 结果型 + 方法成员（非 Component）。
- `generalize` 若 FD 与解析共用“选根 + 恢复几何”外壳 → 抽 private helper 一次服务两者。

## Files to create or modify

- `include/atomiverse/cpu/scf/rhf.h`
- `src/cpu/scf/rhf.cpp`（及 TDA/梯度相关 TU）
- `include/atomiverse/cpu/scf/internal/` — analytic TDA gradient workspace (new as needed)
- `tests/cpu/scf/test_tda_forces.cpp` — 扩展解析 vs FD
- `tests/cpu/scf/test_analytic_tda_forces.cpp` (new)
- `regressions/scf/h2_analytic_tda_force/reg_h2_analytic_tda_force.cpp` (new) + `REFERENCE.md` (new)

## Tasks

- [ ] Write failing unit tests for analytic TDA force vs FD (tests/cpu/scf/test_analytic_tda_forces.cpp)
- [ ] Implement analytic CIS/TDA nuclear gradient in cpu/scf (private workspace + RHF::tda_forces path)
- [ ] Generalize shared root-selection/geometry restore helper if FD and analytic share it
- [ ] Document API units [eV/Å] and analytic default in rhf.h (Doxygen)
- [ ] Add regression regressions/scf/h2_analytic_tda_force (hard-coded force components or max|ΔF_analytic−FD|)
- [ ] Verify H2 S1 analytic forces match FD within documented tol
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：解析力 vs 同几何 FD（h=1e-3 Å）max 分量误差门槛（如 1e-3 eV/Å 量级，以实现时标定并硬编码）。
- **Unit**：平移不变性 \(\sum_a \mathbf{F}_a \approx 0\)。
- **Regression**：`h2_analytic_tda_force` 公开 `RHF` API；无第三方运行时。
- 不在 `tests/` 跑完整表面扫描。

## Out of scope

- TDDFT / ALDA 力（仅 RHF TDA/CIS）
- 真 NAC 向量（`true-nac-vectors`）
- 解析 Hessian / 红外 ES
- CUDA
- 开壳层 UHF-CIS 力

---

**Librarian note：** placement `cpu/scf/`；reuse FD path + HF gradient；confidence high。

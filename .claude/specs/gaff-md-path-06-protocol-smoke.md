---
title: GAFF 分段协议冒烟 — min + 短 NVT/NPT
slug: gaff-md-path-06-protocol-smoke
status: approved
chain: gaff-md-path
link: 06/07
created: 2026-07-21
depends_on: gaff-md-path-05-xcheck
---

# GAFF 分段协议冒烟 — min + 短 NVT/NPT

## Summary

用 05 的力场接线在 **短步数** 上跑通 Gudla–Zhang 风格四段协议的 **形状**：`Minimize`（LBFGS 或 CG）→ `md::MD` NVT → NPT 升温斜坡 → NPT 恒温，证明带电荷 GAFF 拓扑在 ATV 驱动下不 NaN、温度/密度有限。这不是 paper 的 10 ns×3 段，而是端到端冒烟（建议总步数 O(10²–10³)），为后续换真实时长与 paper 对照留接口。

## Domain basis

- 协议形状（Gudla & Zhang, *J. Phys. Chem. B* 2024, DOI: 10.1021/acs.jpcb.4c06018）：  
  (0) minimize (1) NVT \(T=400\,\mathrm{K}\) (2) NPT anneal \(400\to1000\to400\,\mathrm{K}\) (3) NPT equil \(400\,\mathrm{K}\)，\(P=1\,\mathrm{atm}\)。  
  本环保留温度设定与 ensemble 切换，**大幅缩短**步数；`dt` 可用 1 fs（与 paper）或 0.5 fs 若稳定性需要。
- 积分：NVT 用 `NoseHoover` 或 MTK-NPT 的 T 通道；NPT 用 `NoseHooverNPT` 各向同性；T-ramp 用已有 `T_target`/`T_stop`。
- 最小化：paper 用 steepest descent；ATV 用 LBFGS/CG（roadmap 已说明）— REFERENCE 注明差异，不要求同路径。
- 成功判据：全程无 NaN/Inf；结束时 \(T>0\)、盒体积有限、`kEnergy` 有限；可选 \(|\Delta E|\) 在 NVE 短探针可接受（非硬门）。

## Design

- **回归**：`regressions/md/gaff_protocol_smoke/`  
  - 复用 05 fixture 或略大一点的周期盒  
  - 四段 `run` 封装在单一 cpp（或小型 helper 仅限该目录）  
  - 退出码 0 + 打印 thermo 摘要
- **不**实现通用 Protocol 类进 public API（避免 LARGE）；协议逻辑可在 reg 文件内联。
- 若 05 fixture 太小导致 NPT 不稳：固定体积 NVT-only 三段仍可；但至少 **一段 NPT** 必须跑通（acceptance 写清）。

## Files to create or modify

- `regressions/md/gaff_protocol_smoke/reg_gaff_protocol_smoke.cpp` (new)
- `regressions/md/gaff_protocol_smoke/REFERENCE.md` (new)
- `regressions/CMakeLists.txt`
- 可复用 05 的 data/coeffs（相对路径 include 或 CMake 复制）

## Tasks

- [ ] Write failing smoke regression that assembles GAFF forces and attempts staged runs
- [ ] Implement minimize + short NVT + short NPT (with T ramp if feasible) without NaN
- [ ] Assert finite T, energy, volume at end; document step counts in REFERENCE.md
- [ ] Register in `regressions/CMakeLists.txt` with reasonable wall-time budget (≤120 s serial)
- [ ] Run full check + test suite

## Testing strategy

- **快乐路径**：reg 退出 0；日志含四段或文档化的缩减段。
- **边界**：初始重叠构型可先 minimize；若仍炸，REFERENCE 记录所需 soft start（不实现 soft-core）。
- **领域**：不比对 paper 密度平台；只证明驱动栈可承载 GAFF。

## Out of scope

- 10 ns 级生产与 multi-replica。
- \(R_g\) / 密度与 paper 图数值对齐。
- Packmol/molpack 打包集成。
- 公共 `Protocol` 组件库。

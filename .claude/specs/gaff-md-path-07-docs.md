---
title: 文档与路线图 — GAFF/Amber 经典路径状态
slug: gaff-md-path-07-docs
status: approved
chain: gaff-md-path
link: 07/07
created: 2026-07-21
depends_on: gaff-md-path-06-protocol-smoke
---

# 文档与路线图 — GAFF/Amber 经典路径状态

## Summary

更新力场/MD 路线图与用户可见文档，明确 Atomiverse 现已支持 **Amber special_bonds 缩放、full data 电荷读入、PME 实空间 specials、单 term Fourier→Charmm、GAFF 片段回归**，并列出仍缺项（多 term Fourier、完整 `.ff` 解析、CUDA Amber、paper 级 30 ns）。避免路线图继续写「LJ 无多 type / 无 impropers」等过时句导致能力误判。

## Domain basis

无新物理。文档准确性本身是交付物：与 01–06 已实现行为一致。

## Design

- `docs/roadmap/force-fields.md`：修正 special_bonds / multi-type LJ / CoulLong 叙述；增加 GAFF/Amber 状态小节。
- `docs/roadmap/molecular-dynamics.md`：如有「无电荷聚合物」暗示则更正。
- `docs/roadmap/index.md` Classical force fields 子弹要点。
- 教程：在 `docs/tutorial/ch10-classical-md.md` 或 ch11 末增加 **短节**「Charged molecular FF (GAFF path)」链到回归目录与 Amber specials API；不新开整章除非篇幅需要。
- `docs/tutorial/index.md` 一行索引（若新节）。

## Files to create or modify

- `docs/roadmap/force-fields.md`
- `docs/roadmap/molecular-dynamics.md`
- `docs/roadmap/index.md`
- `docs/tutorial/ch10-classical-md.md`（或 ch11，择一）
- `docs/tutorial/index.md`（若需要）

## Tasks

- [ ] Update `docs/roadmap/force-fields.md` with Amber specials, multi-type LJ, CoulLong specials, GAFF fragment regression status
- [ ] Update `docs/roadmap/index.md` classical FF bullets for accuracy
- [ ] Add short tutorial subsection linking API + `regressions/forcefield/gaff_peo_fragment`
- [ ] Cross-check docs against implemented APIs (no claiming multi-term fourier or 30 ns production)
- [ ] Run full check + test suite (docs-only: build still green)

## Testing strategy

- **docs 型**：人工/验收读 roadmap 与代码 API 一致；不引入新 runtime 测试。
- 构建不破（无死链强制工具则抽查相对路径）。

## Out of scope

- 新功能实现。
- 完整用户手册式 GAFF 教程长文。
- 更新 molpack-paper。

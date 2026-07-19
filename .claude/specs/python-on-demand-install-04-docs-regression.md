---
slug: python-on-demand-install-04-docs-regression
title: docs and primary install regression
status: approved
grilled: true
created: 2026-07-19
chain: python-on-demand-install
phase: 04-docs-regression
---

# docs and primary install regression (04)

## Summary

在 01–03 落地后，重写 `python/docs/getting_started.md` 安装章节（根 `README.md` 加指针），说明两阶段安装、`--config` XOR 内联 `--md`/`--global`、版本锁、`--src`/`--url`、模块列表、`--hf` 拒绝、单环境单二进制覆盖；并新增 `regressions/python/install_cli` 库外用例（help/互斥/版本 mismatch/`--hf`，无完整原生编译）。文档字段与 03 schema 一致。

## Domain basis

无物理方程；literature N/A。

## Design

| 实体 | 角色 |
|------|------|
| `python/docs/getting_started.md` | 重写 Installation |
| `README.md` | 指针到 getting_started |
| `regressions/python/install_cli/` | 链上主 e2e（shell + REFERENCE + fixtures） |
| `regressions/CMakeLists.txt` | `run_reg_python_install_cli` |

文档契约：纯 pip → install；互斥；版本锁；源定位；模块枚举；`--hf` 拒绝；覆盖语义。

### Reuse decision

- **generalize** `python/docs/getting_started.md` Installation
- **new** `regressions/python/install_cli/`（launcher shell pattern only）

## Files to create or modify

- `python/docs/getting_started.md`
- `README.md`
- `regressions/python/install_cli/reg_install_cli.sh` (new)
- `regressions/python/install_cli/REFERENCE.md` (new)
- `regressions/python/install_cli/fixtures/minimal_ok.toml` (new)
- `regressions/CMakeLists.txt`
- `regressions/README.md`

## Tasks

- [ ] Rewrite Installation section in `python/docs/getting_started.md` for two-phase install
- [ ] Add brief Python on-demand install pointer in `README.md`
- [ ] Write `REFERENCE.md` in `regressions/python/install_cli/`
- [ ] Add regression example `regressions/python/install_cli/reg_install_cli.sh` (public CLI only; no full native compile)
- [ ] Wire `run_reg_python_install_cli` in `regressions/CMakeLists.txt` and row in `regressions/README.md`
- [ ] Verify docs field names match 03 unified config schema
- [ ] Run full check + test suite

## Testing strategy

文档可复制命令；回归 ≤60s：help、互斥、`--hf`、version mismatch；无完整编译。

## Out of scope

实现 CMake/CLI；manylinux；启用 HF；改 MD/AM1 教程主体

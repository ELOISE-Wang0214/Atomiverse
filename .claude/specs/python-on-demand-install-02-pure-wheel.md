---
slug: python-on-demand-install-02-pure-wheel
title: pure-Python wheel and soft extension import
status: approved
grilled: true
created: 2026-07-19
chain: python-on-demand-install
phase: 02-pure-wheel
---

# pure-Python wheel and soft extension import (02)

## Summary

`pip install atomiverse` 交付 **纯 Python wheel**（安装期不编译 C++）：保留 scikit-build-core 为 build-backend，默认 `wheel.cmake = false`，并将 nanobind 移出默认构建依赖。缺少 `_atomiverse` 时包仍可导入、暴露 `__version__` 与 `python -m atomiverse`；计算路径抛出含 `python -m atomiverse install` 的 `ImportError`。扩展存在时行为与现网一致（含 PW `hasattr` 可选导出）。完整 install CLI 由 03 实现；本 spec 固定 pure 交付与 soft-import 契约，`__main__` 可为 stub。

## Domain basis

无物理方程。验证对象是导入与打包契约。

## Design

| 实体 | 变更 |
|---|---|
| `python/pyproject.toml` | 默认 pure；保留 `cmake.define` 映射 |
| `atomiverse._binding` | soft-import + `require_cpp` / 哨兵 `_cpp` |
| 包根 | soft-load；始终 `__version__`；PEP 562 `__getattr__` |
| `__main__.py` | stub（03 再填子命令） |

**错误文案**：说明扩展未装；主提示 `python -m atomiverse install`；cmake 手编路径可作次要补充。

### Reuse decision

- **reuse** `[tool.scikit-build]` / build-backend
- **reuse** `cmake.define` 映射（pure 跳过执行）
- **generalize** `_binding.py`
- **pattern** PW `hasattr`
- **new** `__main__.py` stub — 非 app/ MPI launcher

## Files to create or modify

- `python/pyproject.toml`
- `python/src/atomiverse/_binding.py`
- `python/src/atomiverse/__init__.py`
- `python/src/atomiverse/__main__.py` (new)
- `python/tests/test_pure_import.py` (new)
- `regressions/python/pure_import/reg_pure_import.sh` (new)
- `regressions/python/pure_import/REFERENCE.md` (new)
- `regressions/CMakeLists.txt`

## Tasks

- [ ] Write failing tests for soft-binding and pure package import (`python/tests/test_pure_import.py`)
- [ ] Generalize `_binding.py` soft-import, `extension_available` / `require_cpp`, and install-message `ImportError`
- [ ] Soft-load package root with `__version__` and PEP 562 compute `__getattr__`; preserve PW `hasattr` when extension present
- [ ] Add `__main__.py` stub so `python -m atomiverse` prints install-oriented guidance
- [ ] Configure pure-Python wheel defaults in `python/pyproject.toml` (`wheel.cmake = false`; drop nanobind from `build-system.requires`; keep `cmake.define`)
- [ ] Add regression example `regressions/python/pure_import/reg_pure_import.sh` plus `REFERENCE.md`
- [ ] Wire Python-only `run_reg_python_pure_import` into `regressions/CMakeLists.txt`
- [ ] Run full check + test suite

## Testing strategy

扩展缺失：import 成功、`__version__` 可读、计算路径 ImportError 含 install 文案、`python -m atomiverse` 可启动。扩展存在：现有 surface 测试不回退。纯 wheel 构建无 CMake/无 so。Regression：隔离纯 Python 路径。

## Out of scope

- CMake 模块门（01）、完整 install CLI（03）、manylinux、改 C++ 绑定表面

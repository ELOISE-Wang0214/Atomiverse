---
slug: python-on-demand-install-03-install-cli
title: on-demand install CLI
status: approved
grilled: true
created: 2026-07-19
chain: python-on-demand-install
phase: 03-install-cli
---

# on-demand install CLI (03)

## Summary

用户在纯 Python 轮子（02）上运行 `python -m atomiverse install`，按统一配置（`--config path.toml` **XOR** 内联 `--md k:v,…` / `--global k:v,…`）解析友好键、闭包自动打开依赖、映射到 01 的 `ATV_*`，再经 scikit-build/`pip install --force-reinstall --no-deps` 从与已装包版本一致的源码把单一 nanobind 扩展编入 site-packages。空安装、`--hf`、绑定键走 `-D`、模式混用一律 fail-loud。`status` 汇报模块与 define 指纹。默认测试覆盖解析/互斥/帮助，不强制完整 C++ 编译。

## Domain basis

无物理方程。

## Design

**入口**：`__main__.py` → `cli.main`。CLI **不得**在 import 时加载 `_binding`/`_cpp`。

### 互斥模式

| 模式 | 入口 | 禁止 |
|------|------|------|
| A | `install --config PATH.toml` | 任何 `--md`/`--pw`/`--global` |
| B | `install --md [k:v,…] --pw … --global k:v,…` | `--config` |

空 install → fail-loud。

### 统一配置 schema

```toml
[modules]
md = true
pw = false

[md]
mpi = true

[global]
# real_storage = "double"

[cmake.unbound]
# ONLY unbound keys
```

内联：`--md mpi:yes,key:value`；逗号分项；**首个 `:`** 分键值；bool：`yes/no/true/false/on/off/1/0`。模块 flag 出现即启用。

**绑定表** `BINDING_TABLE`：友好键 → `ATV_*`（01 名）。跨模块只走 `global`/`--global`。绑定键禁止 `-D`/unbound；表外键才可 unbound。`--hf` 拒绝直至 SCF Python 绑定。

**闭包**：与 01 CMake DAG 同语义（md→potentials+neighbor 等），表驱动，注释指向 01。

**源码**：`importlib.metadata.version("atomiverse")`；无 `--version`；默认  
`https://github.com/Roy-Kid/Atomiverse/archive/refs/tags/v{version}.tar.gz`；`--src`/`--url` 覆盖；校验源树 `python/pyproject.toml` version 一致。

**构建**：仅  
`pip install --force-reinstall --no-deps` + `-C cmake.define.*`，强制 `ATV_BUILD_PYTHON=ON`，覆盖单一 `_atomiverse`。

**status**：启用模块 + 稳定 define 指纹。

### Reuse decision

- **reuse** scikit-build cmake.define 路径
- **reuse** 单一 `_atomiverse` so
- **new** `python/src/atomiverse/cli/`
- **reuse** 02 soft `_binding`（CLI 热路径不依赖）
- **reuse** 拒绝 HF 政策

## Files to create or modify

- `python/src/atomiverse/__main__.py` (generalize 02 stub)
- `python/src/atomiverse/cli/__init__.py` (new)
- `python/src/atomiverse/cli/main.py` (new)
- `python/src/atomiverse/cli/config.py` (new)
- `python/src/atomiverse/cli/source.py` (new)
- `python/src/atomiverse/cli/install.py` (new)
- `python/src/atomiverse/cli/status.py` (new)
- `python/tests/test_install_cli.py` (new)
- `regressions/python/install_cli_dryrun/reg_install_cli_dryrun.py` (new)
- `regressions/python/install_cli_dryrun/REFERENCE.md` (new)

## Tasks

- [ ] Write failing tests for InstallConfig resolve, mutual exclusion, HF/empty/bound-`-D` guards (`python/tests/test_install_cli.py`)
- [ ] Implement binding table, TOML/inline parse, auto-enable closure, and `resolve_defines` in `python/src/atomiverse/cli/config.py`
- [ ] Implement SourceSpec (importlib version, default tarball, `--src`/`--url`, pyproject version gate) in `python/src/atomiverse/cli/source.py`
- [ ] Implement install orchestration (`pip` + cmake.define, always `ATV_BUILD_PYTHON=ON`) in `python/src/atomiverse/cli/install.py`
- [ ] Implement status subcommand in `python/src/atomiverse/cli/status.py`
- [ ] Wire `python -m atomiverse` argparse entry (`install` / `status`) in `__main__.py` and `cli/main.py`
- [ ] Add regression example `regressions/python/install_cli_dryrun/reg_install_cli_dryrun.py` (no C++ compile)
- [ ] Run full check + test suite

## Testing strategy

Unit：解析、闭包、互斥、HF、版本 mismatch、绑定 vs -D。Regression dry-run：help、互斥、空/HF、define 集合含闭包与 PYTHON=ON；不调用真实 pip/C++。

## Out of scope

manylinux；多 profile 缓存；手写 cmake；改 MPI launcher；实现 SCF 绑定；改 01 选项定义本身；默认 ctest 全量编译

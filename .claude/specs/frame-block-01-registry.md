---
title: Frame block 化第一步 — 字段注册表 + 编译期 ColKey + mapping API（成员暂不删）
slug: frame-block-01-registry
status: draft
science_required: false
created: 2026-07-14
acceptance: frame-block-01-registry.acceptance.md
depends_on: frame-organization (closed 2026-07-14, verdict in docs/developer/frame-block-redesign.md)
---

# frame-block-01 — registry + ColKey + mapping API

## Summary

执行已签字的 Frame 裁决（变体 C 路线）的第一步：**建立字段注册表与
mapping API，物理存储一个字节不动**。`cpu::Atom` 的命名成员原地保留；
注册表把它们注册成非拥有的 Column 视图，Frame 获得 molrs 同构的
mapping 语义（存在性检查、类型化取列、扩展列）。所有"每处各抄一遍
字段清单"的站点改为遍历注册表——probe 4 的 grow 数据毁灭 bug 这一类
错误从机制上消灭。

零热路径风险：447 个字段访问点**一个都不改**（那是 02+ 的逐目录迁移）。

## Design（细节见 docs/developer/frame-block-redesign.md 变体 A/B 节）

- **字段表**（`core/` 下，backend 无关）：每字段一行——molrs schema 键名
  （`"atoms"/"x"` …，与成员 1:1）、dtype（ATV 精度档：StorageReal/AccReal/
  TagInt/LocalInt/int8）、stride（1D/2D 宽列）、语义标签（随行迁移？入
  exchange payload？入 write_back 契约？拓扑 lazy 组？）。
- **`ColKey<T>`**：编译期常量键对象（带 dtype），`frame.span(kX)` 在 launch
  边界解析一次返回 span；裸字符串查询仅用于泛化路径（FFI、拓展列）。
- **CPU 视图注册**：`cpu::Atom` 构造时把成员注册进映射（非拥有）；扩展列
  （未来的 per-atom 载荷）由注册表拥有——`has_bonded_topology_`/
  `ensure_bonded_storage()` 语义并入"拓扑 lazy 组"的存在性。
- **簿记 = typed meta，一个区域**（2026-07-14 裁决修订）：Frame 形状与 molrs
  完全同构 `{ blocks, meta, simbox }`；`meta` 是 molrs `string→string` 的类型
  扩展（`string→MetaValue` variant），装 nlocal/nghost/nmax、SystemScalars、
  tag→row 索引（引擎本地标签，永不序列化）；Box 走 simbox 槽。没有第二个
  簿记区域。本步（01）先定义 `MetaValue` + meta 骨架并把 SystemScalars/counts
  迁入；热读取（`nlocal` 等）经 launch 边界一次解析的引用，meta 查键计入
  debug 解析计数断言。上游（molrs typed meta）非阻塞，FFI 导出降级 string 兼容。
- **单清单派生**：`grow`/`set_capacity`/`copy_row`/`pack_exchange`/
  `unpack_exchange`/`materialize`/`write_back` 改为遍历注册表；每站点的
  手抄字段列表删除。CUDA/torch 的镜像迁移在 02+。
- **纪律机器化**：debug 构建下 Frame 维护每 launch 解析计数，
  `step_once` 断言 ≤ K——循环内查键在任何 debug 测试跑中立刻炸。

## Files

- `core/frame_schema.h` — 字段表（molrs 键名 + dtype + 语义标签）
- `cpu/frame/column.h` — `ColKey<T>`、Column 视图、dtype variant
- `cpu/frame/atom.h/.cpp` — 视图注册；grow/set_capacity/copy_row 遍历注册表
- `cpu/comm/atom_pack.cpp` — pack/unpack 遍历注册表（payload 标签筛选）
- `cpu/frame/materialize.h/.cpp` — 泛化键拷贝 + 拓扑特例（特例保留，裁决已知）
- `cpu/driver/`（debug 解析计数断言挂点）
- tests：注册表完整性（每成员有表项、无孤儿）、ColKey 取列、扩展列
  插入/存在性、grow/copy/pack 与手写旧实现的等价金标

## Tasks

- [ ] `core/frame_schema.h` 字段表 + `ColKey<T>`（molrs 键名、ATV dtype 档）
- [ ] CPU Column 视图 + `Frame::has/span/insert` mapping API
- [ ] `grow`/`set_capacity`/`copy_row` 改遍历注册表（删手抄清单）
- [ ] `pack_exchange`/`unpack_exchange` 改遍历注册表（payload 标签）
- [ ] `materialize`/`write_back` 泛化键拷贝（拓扑转换特例保留）
- [ ] debug 每 launch 解析计数断言
- [ ] 单测（注册表完整性 / API / 等价金标）
- [ ] 三道门验收（LJ ⚡ / memory_bytes / test_materialize）

## Out of scope

- 删除任何命名成员、改任何 kernel 字段访问（02+ 逐目录：torch → cuda →
  cpu/md、minimizer、path → cpu/potentials）
- per-atom ⇄ 扁平拓扑转换的泛化（裁决已知特例）
- molrs 侧改动

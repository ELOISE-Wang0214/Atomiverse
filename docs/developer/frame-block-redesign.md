# Frame = mapping[string, Block] 重设计调研报告

2026-07-14 · frame-organization 第三轮方向的专项调研
（用户指令：风格与 molrs 完全对齐，Frame 完全是 mapping[string, Block]，
包括 x/v/f 在内全部 block 化；多后端实现与性能需重新分析）

---

## 1. 对齐目标实录：molrs 的 Frame/Block 到底长什么样

读的是 `molrs/src/core/store/{frame.rs, block/mod.rs, block/column.rs, block/dtype.rs}`（dev 主树，非 worktree）。

```
Frame  { map: HashMap<String, Block>,  meta: HashMap<String,String>,  simbox: Option<SimBox> }
Block  { map: HashMap<String, Column>, nrows: Option<usize>, shape: Option<Vec<usize>> }
Column = enum { Float(Arc<ColumnHolder<F>>), Int(…I), Bool, UInt, U8, String }
ColumnHolder<T> = ManuallyDrop<ArrayD<T>> + Option<foreign_keeper>   // COW、O(1) clone
```

要点：

- **两级映射**。Frame 是 `string → Block`；Block 是 `string → Column` 且强制
  所有列共享 axis-0 长度（`nrows`）。Frame 层**不**强制跨 block 一致性。
- **Column 是 dtype 枚举 + `Arc` 共享 + ndarray 动态维**（`ArrayD<T>`）。
  clone 是引用计数；写走 `realize_owned_mut`（对 foreign 缓冲 copy-on-write，
  这是 numpy 零拷贝借入的机制）。
- **dtype 集合**：Float(F)/Int(I)/Bool/UInt/U8/String，F/I 是编译期标量别名。
- **shape 元数据**：Block 可声明自己是 N 维（体积网格），列仍按行主序平铺。
- 交换 schema 的既定键名（`materialize.cpp` 实测）：block `"atoms"` 的列
  `x/y/z/type/vx/…/tag/image_x…`；block `"bonds"`/`"angles"`/`"dihedrals"`。
  **与 `cpu::Atom` 字段名 1:1** —— 对齐没有命名迁移成本。
- 一个偏差点必须记录：molrs 的 `Frame.meta` 是 `string→string`；ATV 的
  `nlocal`/`nghost`/`nmax`、`tag→row map`、`SystemScalars` 不是列也不是字符串
  元数据。**全对齐形态下这些要有一个 typed-attrs 位置**（见 §5 开放问题）。

## 2. 消费面盘点（迁移成本的硬数字）

对 `src/`（不含测试）的实测统计：

| 指标 | 数字 | 含义 |
|---|---|---|
| 直接 `atoms.<field>[i]` / `(i,d)` 访问点 | **447 处**（.cpp/.cu） | 全 block 化的改写面 |
| launch 边界 `.data()` 裸指针提升 | **83 处** | 已遵守"解析一次"纪律的部分 |
| 含字段访问的头文件 | 8 个 | 模板/内联路径 |
| 手抄字段清单的站点 | **7+**：`grow`/`set_capacity`/`copy_row`/`pack_exchange`/`unpack_exchange`/`materialize`/`write_back`/`to_cuda`/`to_cpu` | probe 4 的 grow 数据毁灭 bug 的温床 |
| 访问最重的文件 | `path/neb.cpp` 42、`md/init.cpp` 33、`minimizer/cg.cpp` 24、`nose_hoover_npt.cpp` 24、`cuda/comm/multi_gpu.cu` 23 … | 迁移排序参考 |

热路径访问模式二分：

- **纪律良好**（`pair_lj_cut.cpp:264-270` 典型）：launch 开头 `.data()` 解析成
  裸指针，内层循环只见指针。xsimd ket loop、邻居表、comm pack 都是这一类。
- **散布访问**（积分器、minimizer、NEB、init）：循环体内直接 `atoms.x[i]`。
  xtensor `operator[]` 内联后同样是指针算术，性能等价——但**换成 map 查找就不
  等价了**，这些是全 block 化时必须逐处提升的地方。

三个后端的现状表示：

| 后端 | 存储 | 备注 |
|---|---|---|
| `cpu::Atom` | xtensor 1D/2D 命名成员 | 拓扑宽数组 per-atom 布局 + 运行时 stride |
| `cuda::Frame` | 裸设备指针命名成员（`real_t* x` …） | kernel 以指针实参接收——**天然就是 launch 边界解析** |
| `torch::Frame` | `at::Tensor` 命名成员 | `at::Tensor` 本身就是一个 refcounted、dtype/device 标注的 Block 列 |

## 3. 先例调研

- **LAMMPS**（[fix property/atom](https://docs.lammps.org/fix_property_atom.html)、
  [Atom/AtomVec 开发者文档](https://docs.lammps.org/Classes_atom.html)、
  [per-atom 数据访问](https://docs.lammps.org/Developer_atom.html)）：
  核心字段（x/v/f/type/tag/image）是**固定成员**，kernel 直接访问；
  运行时扩展列走 `Atom::add_custom`（名字注册、索引稳定、可选 ghost 通信/重启）。
  **关键先例**：现代 `AtomVec` 用**字符串字段清单**（`fields_grow`,
  `fields_copy`, `fields_exchange`, `fields_border`, …）驱动 grow/copy/pack 的
  生成——LAMMPS 用注册表解决了与我们 7 处手抄完全同款的问题，
  **但没有把热路径访问 block 化**。
- **HOOMD-blue**（[ParticleData.h](https://github.com/joaander/hoomd-blue/blob/master/hoomd/ParticleData.h)、
  [hoomd.data 文档](https://hoomd-blue.readthedocs.io/en/v4.2.0/module-hoomd-data.html)、
  [统一内存论文](https://par.nsf.gov/servlets/purl/10123335)）：
  内部按属性 SoA（GlobalArray，CPU/GPU 双驻留），**命名/映射语义只出现在快照
  API 边界**（`LocalSnapshot` 暴露 named arrays）；热路径仍是句柄+裸指针。
- **TensorDict**（[overview](https://docs.pytorch.org/tensordict/main/overview.html)、
  [GitHub](https://github.com/pytorch/tensordict)）：纯 mapping 模型的代表。
  官方文档明说 `get` 已成"expensive operation"、热循环中反复取键有可测开销，
  依赖 `torch.compile` 消除——**C++ 没有 JIT，等价物只能是人工的 launch 边界
  提升纪律**。
- 共同结论：**没有一个生产级 MD 引擎在 kernel 内层用字符串键取数据**；
  mapping 语义都放在 (a) 簿记（grow/copy/pack 生成）和 (b) API 边界（快照/
  交换 schema）。这不否定全 block 化——它指明全 block 化的**安全实现形态**：
  物理存储可以是 mapping，热路径访问必须经过解析边界。

## 4. 微基准（本机实测，GCC 14.3 `-O3 -march=native`，1M 原子，50 次取均值）

```
launch-resolved raw pointers        0.483 ms/pass    ← 与今天的直接成员完全同速
per-access map lookup (forbidden)   6.600 ms/pass    ← 13.7× 慢
single string-key resolution        7.8 ns/lookup
```

- **每 launch 解析 ~15 列 ≈ 0.12 μs**，对 ms 级的步长是 1e-4 量级——归零。
- **失守一处**（循环体内查 map）= 该 kernel 13.7×。这是全 block 化唯一的
  性能悬崖，且是纪律问题而非架构问题。

## 5. 设计变体

### 变体 A —— registry-backed struct（LAMMPS AtomVec 形态，映射在 API 层）

- **存储**：`cpu::Atom` 的命名 xtensor 成员原地不动。
- **新增**：一张 constexpr 字段表（名字、dtype、stride、语义标签：随行迁移？
  入 exchange payload？入 write_back 契约？），每个后端把成员注册成非拥有的
  `Column` 视图 → `Frame` 获得 `has("atoms","rho")` / `get<f64>("atoms","x")` /
  `insert(...)`（扩展列由注册表拥有）。
- **grow/set_capacity/copy_row/pack/unpack/materialize 全部改为遍历注册表**
  ——7 处手抄清单变 1 处。
- 热路径：**447 个访问点一个都不改**。
- 定位：mapping 语义完整可用（存在性检查、扩展列、FFI 泛化），但物理上核心
  字段仍是成员——**不满足"完全是 mapping"的字面要求**，是它的 API 前缀。

### 变体 B —— 纯 mapping（物理形态即 mapping[string, Block]，用户方向的字面形态）

- **存储**：`Frame = flat_map<string, Block>`；`Block = flat_map<string, Column> + nrows`；
  `Column` 是后端各自的 dtype×存储 variant：
  - cpu: `std::variant<xt::xtensor<StorageReal,1>, xt::xtensor<AccReal,1>, xt::xtensor<TagInt,1>, xt::xtensor<int,1>, xt::xtensor<int8_t,2>, …>`
  - cuda: `{DType, void* dev_ptr, extent}` 的 tagged struct（kernel 收裸指针实参，现状即如此）
  - torch: `at::Tensor` **本身就是 Column**——torch 后端几乎免费
- **访问 API**：强类型键常量 + launch 边界解析：
  ```cpp
  // 键是编译期常量对象（带 dtype），不是裸字符串——查错在类型系统里死
  inline constexpr ColKey<StorageReal> kX{"atoms", "x"};
  auto x = frame.span(kX);            // launch 开头，一次哈希，返回 span
  ```
  内层循环只见 `span`/裸指针（微基准 (a) 路径，零开销）。
- **纪律的机器化**（对付 13.7× 悬崖）：debug 构建下 Frame 维护"本 launch 的
  解析计数"，`step_once` 断言每 launch 解析次数 ≤ K（K≈列数×2）；任何在循环
  里查键的代码在单测的 debug 跑里立刻炸。加上既有 LJ ⚡ 性能门兜底。
- **grow/copy_row/pack/unpack**：遍历 mapping 的单一实现（语义标签同 A），
  grow-bug 类错误机制性消灭。
- **拓扑**：宽 2D 列（`bond_atom[nmax,k]`）用 Block 的 shape/stride 元数据表达，
  re-stride = 列替换；`has("bonds_per_atom")` 取代 `has_bonded_topology_`。
- **簿记的家**（molrs 对齐的必要偏差）：`nlocal/nghost/nmax`、`tag→row map`、
  `SystemScalars`、`Box` 不是列。放 Frame 的 typed attrs 区（molrs `meta` 是
  string→string，装不下）。这是"完全对齐"唯一需要明说的偏差。
- **迁移面**：447 访问点 + 83 hoist 全部改写；8 个头文件的内联路径重排。
  CUDA 最容易（kernel 已经只吃指针实参，改的只是 host 侧组装）；torch 次之；
  CPU 是大头。

### 变体 C —— A ⊂ B 的分步路线（推荐的执行策略，不是第三个终态）

A 的注册表、字段语义标签、`ColKey` API、三道验收门与 B **完全同构**——A 就是
B 去掉"删除命名成员"这最后一步。路线：

1. **A 落地**（一个 spec）：字段表 + registry 视图 + mapping API + 簿记单点化。
   零热路径风险，立即消灭 grow-bug 类，FFI/materialize 泛化。
2. **B 逐目录推进**（每目录一个小 spec，按上面"访问最重文件"倒序——先易后难，
   `torch/` → `cuda/` → `cpu/md、minimizer、path` → `cpu/potentials`）：
   成员访问改 `span(kX)`，目录清空后该成员从 struct 删除、所有权移交 mapping。
3. **终态 B**：`cpu::Atom` 只剩 mapping + attrs，"Frame 完全是 mapping[string,
   Block]"字面成立；`test_materialize` 金标从第一步起每步不动。

## 6. 性能与内存对照表

| 维度 | 今天 | A | B（守纪律） | B（失守一处） |
|---|---|---|---|---|
| 内层循环每原子成本 | 指针算术 | 同 | **同**（微基准 0.483 vs 0.483 ms） | 该 kernel **13.7×** |
| 每 launch 额外成本 | `.data()`×K | 同 + 0 | +7.8 ns×K ≈ 0.1 μs | 同左 |
| xsimd ket loop / min-image 内联 | — | 不动 | 不动（span→裸指针后同码） | 不动 |
| LJ ⚡ 门（≤2%） | — | 稳过 | 过（条件：解析计数断言在） | **必挂**（门兜住） |
| grow/copy/pack 维护 | 7+ 处手抄 | 1 处 | 1 处 | 1 处 |
| 内存 | — | +注册表 ~1 KB | 哈希节点+描述符 ~数 KB（vs 数据 GB 级） | 同左 |
| 迁移面 | — | **0 个访问点** | 447+83 全量 | — |

## 7. 推荐与最大风险

**推荐：按变体 C 的路线走向 B 终态。** 用户方向（完全 mapping）在性能上是
可行的——前提被微基准与先例双重确认：**解析必须停在 launch 边界**。A 一步
先行不是妥协，是让 447 点迁移在注册表、类型化键与三道门都就位之后进行，
每目录的 diff 纯机械且被 debug 断言 + ⚡ 门双重看护。

**最大风险**（按序）：

1. **纪律悬崖**：一处循环内查键 = 13.7×。缓解：debug 解析计数断言（机器化，
   不靠 review）+ LJ ⚡ 门。
2. **簿记偏差**：nlocal/map/scalars/Box 需要 typed attrs，是对 molrs 字面
   对齐的唯一结构性偏差——需要用户认可这个偏差（或推动 molrs 侧加 typed meta）。
3. **拓扑转换特例**：per-atom LAMMPS 布局 ⇄ molrs 扁平 bonds 的转换逻辑
   不能被"泛化键拷贝"吸收，仍是 materialize 里的手写段（今天也是）。
4. **迁移长尾**：447 点分目录推进期间，成员与 mapping 双所有权并存的窗口要
   靠"目录清空即删成员"的节奏压短，避免 A/B 混态长期化。

## 8. 开放问题（裁决前要用户表态）

1. typed attrs 偏差可接受吗？还是推动 molrs 加 typed meta 后再全对齐？
2. dtype 策略：molrs `F` 固定 f64；ATV `StorageReal` 可配 f32。Block 列 dtype
   跟 ATV 精度档（转换在 FFI 边界，现状），还是跟 molrs（放弃 f32 档）？
3. B 终态里 `AccReal`（force 累加器恒 f64）与 `StorageReal` 列共存于同一
   Block——dtype variant 里两种 float 列并存，接受？
4. 键常量的命名域：`"atoms"/"x"` 沿用 molrs 交换 schema（建议，零迁移），
   还是另立 ATV 命名空间？

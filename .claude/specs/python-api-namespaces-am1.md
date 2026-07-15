---
title: Python API 命名空间化重构 + AM1 半经验单点绑定
status: approved
created: 2026-07-15
---

# Python API 命名空间化重构 + AM1 半经验单点绑定

## Summary

本 spec 将 Python 包 `atomiverse` 从扁平布局重构为镜像 C++ 后端目录分层的深层命名空间（`atv.potentials.pair.LJCut`、`atv.md.integrator.VelocityVerlet`、`atv.neighbor.NeighborList` 等），深层路径是唯一实现真源（single source of truth），顶层 `atv.PairLJCut` 等 LAMMPS 规范复合名仅作为同一类对象的便捷再导出保留——包未发布、无外部用户，旧扁平实现布局直接删除，不设 DeprecationWarning 层、不留双实现。同时新增 AM1 半经验方法的 Python 绑定（仅单点：`AM1Config` + `AM1Evaluator` + `single_point` 入口；`AM1Evaluator` 能力为 Energy|Charge|DensityMatrix、无 Force，AM1 梯度/AM1-MD 显式出界），能量落 `frame.scalars.energy` [eV]、Mulliken 电荷经 write_back 契约写回 molrs.Frame 的 charge 列、精细量经 `last_result()` 读取。后端选择以 `device="cpu"` 构造参数表达（PyTorch 约定），不引入 `atv.cuda.*` 命名空间轴；v1 仅接受 `"cpu"`，其余值即刻报错。本 spec 独立可交付（不依赖 launcher 链），但必须先于 `mpi-launcher-embedded-python` 第 04 环合入——04 环的 TaskSpec 镜像记录点在绑定 lambda 内，本 spec 保持 `bind_md.cpp` 每个 `add_*` 为离散 lambda 位点且**零改动**。

## Domain basis

命名空间重构部分为纯软件结构，不含物理。AM1 绑定部分暴露的是**已在 C++ 层验收过的** AM1 半经验方法（Dewar, Zoebisch, Healy, Stewart, *J. Am. Chem. Soc.* **1985**, 107, 3902–3909, DOI [10.1021/ja00299a024](https://doi.org/10.1021/ja00299a024)）；Python 层不新增任何物理数学，只暴露既有 `cpu::semiempirical::AM1Evaluator`。域校验沿用 `tests/cpu/semiempirical/test_am1_scaffold.cpp` 已钉死的参考 oracle（AM1 参考实现 23.2.5 `AM1 1SCF XYZ` 与 AmberTools sqm 单点电荷）：正四面体 CH₄（C–H 沿体对角、h=0.629118 Å，即 C–H≈1.0897 Å）的 Mulliken 电荷 C = −0.262125 e、H = +0.065531 e。Python 侧验收**不引入新参考值**（与 `test_surface.py` 的既定惯例一致），容差取 2.5e-3 e——覆盖 `ATV_REAL_STORAGE=float` 轮子下 frame charge 列的单精度存储损失（`last_result().mulliken_charges` 恒为 double，不受影响）。单位：`frame.scalars.energy` [eV]；`last_result()` 上 `electronic_energy_hartree` [Hartree]、`heat_of_formation_kcal_mol` [kcal/mol]；电荷 [e]。

## Design

### 命名空间：深层为真源，顶层仅再导出

深层路径逐一镜像 CLAUDE.md「Backend directory layering」的目录分层（不是 C++ namespace——namespace 与目录尚未对齐是 C++ 侧自己声明的待办；Python 按用户裁决对齐**目录**）。深层用去前缀名，顶层用 LAMMPS 规范复合名，**两者是同一个类对象**（`atv.PairLJCut is atv.potentials.pair.LJCut == True`）。映射表（本表同时写入 `__init__.py` docstring）：

| 深层规范路径 | 顶层再导出 | C++ 符号（目录） |
|---|---|---|
| `atv.md.MD` | `atv.MD` | `cpu::md::MD`（md/） |
| `atv.md.integrator.VelocityVerlet` | `atv.VelocityVerlet` | `cpu::md::VelocityVerlet`（md/） |
| `atv.md.thermostat.NoseHoover` | `atv.NoseHoover` | `cpu::md::NoseHoover`（md/） |
| `atv.md.InitVelocity` | `atv.InitVelocity` | `cpu::md::InitVelocity`（md/init.h，无更细类目） |
| `atv.neighbor.NeighborList` | `atv.NeighborList` | `cpu::md::NeighborList<>`（neighbor/ 是后端根目录，不挂 md 下——挂 md 会复刻 C++ 自认的 namespace 错位） |
| `atv.potentials.pair.LJCut` | `atv.PairLJCut` | `cpu::PairLJCut`（potentials/） |
| `atv.potentials.bond.Harmonic` | `atv.BondHarmonic` | `cpu::BondHarmonic`（potentials/） |
| `atv.potentials.angle.Harmonic` | `atv.AngleHarmonic` | `cpu::AngleHarmonic`（potentials/） |
| `atv.potentials.dihedral.MultiHarmonic` | `atv.DihedralMultiHarmonic` | `cpu::DihedralMultiHarmonic`（potentials/） |
| `atv.recorder.Thermo` | `atv.Thermo` | `cpu::Thermo`（recorder/） |
| `atv.pw.pw_scf` / `atv.pw.PWResult` | `atv.pw_scf` / `atv.PWResult` | `cpu::pw`（pw/，既有模块原样保留） |
| `atv.semiempirical.AM1Config` | `atv.AM1Config` | `cpu::semiempirical::AM1Config`（semiempirical/） |
| `atv.semiempirical.AM1Evaluator` | `atv.AM1Evaluator` | `cpu::semiempirical::AM1Evaluator`（semiempirical/） |
| `atv.semiempirical.AM1Result` | `atv.AM1Result` | `cpu::semiempirical::AM1Result`（semiempirical/） |
| `atv.driver.single_point` | `atv.single_point` | `cpu::single_point`（driver/single_point.h） |

实施形态：`atomiverse/md.py` 移入 `atomiverse/md/__init__.py`（MD fluent 包装类 + `InitVelocity`，import 路径 `atomiverse.md` 保持可用），`integrator.py`/`thermostat.py` 为其子模块；`potentials/` 四子模块各持一个 `_cpp` 类的去前缀别名；`neighbor.py`、`recorder.py`、`driver.py`、`semiempirical/` 同理。顶层 `__init__.py` 只做 `from atomiverse.potentials.pair import LJCut as PairLJCut` 式再导出，`__all__` 恰为 16 个顶层名（上表右列 + `single_point`）。**nanobind 扩展 `_atomiverse` 保持扁平**：`_cpp` 是私有实现细节，类名保持 LAMMPS 规范名，这同时保证 04-embed 的 `MdHandle` 计划与「每个 add_* 是离散 lambda 位点」不受本 spec 扰动——`bind_md.cpp` 本 spec 零改动。

### device= 是构造参数，不是命名空间轴（用户裁决，今日给定）

`atv.md.MD(dt=None, device="cpu")`、`atv.driver.single_point(frame, evaluator, device="cpu")`。理由记录：C++ 分 `cpu::`/`cuda::` 是因为它们是不同的编译期类型；Python 若以命名空间分后端，每个 import 路径都要按后端翻倍，而 PyTorch 的 device 字符串约定已是用户的既成预期；既有 `python/src/atomiverse/_device.py` 的 molpy device-handler 桩正预示了这一方向；且与 launcher 链的 TaskSpec 模式同构（后端无关参数 + 物化器在某后端实例化，未来 TaskSpec 加 device 字段即可扩展 CUDA 分发）。v1 仅 `"cpu"` 合法；传 `"cuda"`（或任何其他值）抛 `ValueError`，消息明确说明 CUDA 绑定尚未提供（no-fallback 规则；`cuda::md::MD` 在 C++ 存在但未绑定，CUDA 绑定是未来独立 spec）。校验在 Python 包装层完成（v1 无 C++ 改动必要）。Frame 的设备驻留遵循 molrs/molpy 的 device 语义，atv 不自造一套。

### AM1 绑定：绑 `cpu::single_point`，不绑 `Evaluate`（取最小者）

新增 `python/src/bindings/bind_semiempirical.cpp`（无条件编译——semiempirical 源无 `ATV_BUILD_*` 门）：

- `AM1Config`（即 `NDDOConfig`，`nb::class_` 名 "AM1Config"）：kwargs 构造 `AM1Config(charge=0, multiplicity=1, max_iter=200, density_tol=1e-8, energy_tol=1e-10)` + 各字段 `def_rw`；`method` 不暴露（`AM1Solver` 构造时固定为 AM1）。
- `AM1Evaluator`：`nb::init<AM1Config>`（含默认构造），`last_result()`、`config()` 只读访问器。
- `AM1Result`（即 `NDDOResult`）：只读属性 `converged`、`iterations`、`electronic_energy_hartree`、`heat_of_formation_kcal_mol`、`mulliken_charges`（→ list[float]）、`reference`。
- `m.def("single_point", ...)`：lambda 接收 `(nb::object frame_obj, AM1Evaluator& ev)`——`borrow_molrs_frame` → `cpu::materialize` → GIL 释放窗口内 `cpu::single_point(work, ev)` → `cpu::write_back(work, *fref)`。与 `bind_pw.cpp` 的差异是**必须 write_back**：AM1 声明 `Capability::Charge`，电荷按已钉死的 write_back round-trip 契约（charge 在回写集内，golden-tested）回到同一 molrs.Frame。evaluator 以左值借入，调用后用户可读 `ev.last_result()`。

选 `single_point` 而非绑 `Evaluate` 驱动的理由：(1) `cpu/driver/single_point.h` 自己就钦定它是 `Evaluate + Evaluator` 模式的合法糖，且明令禁止 `SinglePoint` 类——单点就是一次 Evaluate 任务；(2) 绑 `Evaluate` 需要 `evaluate::Phase` 枚举 + 新 PhaseBuilder + `add_am1_evaluator` 三件套，而扁平 `_cpp` 里 `Phase` 名已被 `md::Phase` 占用，被迫改名或重构扩展模块结构，收益为零；(3) 一个 `m.def` lambda 对三个绑定类，最小面即诚实面；未来 SCF/PW evaluator 可按类型加 `single_point` 重载，无需返工。Python 侧 `atv.driver.single_point(frame, evaluator, device="cpu")` 包装 `_cpp.single_point`。

### 测试基础设施（已核实）

Python 侧既有 pytest 于 `python/tests/`（conftest + test_md + test_surface，经 `pyproject.toml testpaths=["tests"]` 驱动，**不经 ctest**）——nanobind 级检查与 molpy 互操作检查都落这里；molpy 互操作用 `pytest.importorskip("molpy")` 守门，验收在 `~/work/.x86_64` venv（x86_64 登录节点；GPU 节点用兄弟 `.aarch64`，绝不混用）实跑。C++ `tests/` 三层铁律不受触碰（本 spec 不向 ctest 注册任何东西）。互操作的交换通货是双方共同包装的 molrs Frame：molpy 建体系 → 取其 molrs.Frame → atv 跑短 MD。

## Files to create or modify

- `python/src/atomiverse/__init__.py` — 重写：16 名顶层再导出 + 映射表 docstring
- `python/src/atomiverse/md.py` — **删除**（实现移入 `md/` 子包）
- `python/src/atomiverse/md/__init__.py` (new) — MD 包装类（+ `device=` 参数）+ `InitVelocity`
- `python/src/atomiverse/md/integrator.py` (new) — `VelocityVerlet`
- `python/src/atomiverse/md/thermostat.py` (new) — `NoseHoover`
- `python/src/atomiverse/neighbor.py` (new) — `NeighborList`
- `python/src/atomiverse/potentials/__init__.py` (new)
- `python/src/atomiverse/potentials/pair.py` (new) — `LJCut`
- `python/src/atomiverse/potentials/bond.py` (new) — `Harmonic`
- `python/src/atomiverse/potentials/angle.py` (new) — `Harmonic`
- `python/src/atomiverse/potentials/dihedral.py` (new) — `MultiHarmonic`
- `python/src/atomiverse/recorder.py` (new) — `Thermo`
- `python/src/atomiverse/driver.py` (new) — `single_point(frame, evaluator, device="cpu")`
- `python/src/atomiverse/semiempirical/__init__.py` (new) — `AM1Config`/`AM1Evaluator`/`AM1Result`
- `python/src/bindings/bind_semiempirical.cpp` (new)
- `python/src/bindings/bind.cpp` — 注册 `bind_semiempirical(m)`
- `CMakeLists.txt` — `nanobind_add_module(_atomiverse ...)` 源列表（约 L749）加 `bind_semiempirical.cpp`
- `python/tests/test_namespaces.py` (new)
- `python/tests/test_semiempirical.py` (new)
- `python/tests/test_molpy_interop.py` (new)
- `python/tests/test_surface.py` — `__all__`/身份契约更新到新面
- `python/tests/test_md.py` — 深层规范路径改写 + device= 用例
- `python/tests/conftest.py` — fixture 改用深层规范路径
- `python/docs/getting_started.md` — 示例更新为规范深层路径 + 顶层便捷两写法

## Tasks

- [ ] Write failing tests for the namespaced package（python/tests/test_namespaces.py：全映射表逐条 `is` 身份断言、深层模块可 import、`atomiverse/md.py` 扁平实现不再是类定义处、`MD(device="cpu")` 默认可跑 / `device="cuda"` 抛 ValueError 且消息点名 CUDA 绑定不可用；同步更新 python/tests/test_surface.py 的 `__all__`==16 名契约）
- [ ] Implement the namespaced package layout（新建 md/{__init__,integrator,thermostat}.py、neighbor.py、potentials/{__init__,pair,bond,angle,dihedral}.py、recorder.py、driver.py；删除扁平 md.py；重写顶层 __init__.py 再导出 + 映射表 docstring）
- [ ] Implement the device= parameter in atv.md.MD and atv.driver.single_point（默认 "cpu"，其余值抛 ValueError 并说明 CUDA 绑定为未来 spec；纯 Python 层校验）
- [ ] Write failing tests for the AM1 Python surface（python/tests/test_semiempirical.py：AM1Config kwargs 构造与字段、single_point 后 CH₄ 电荷 vs 钉死 oracle（容差 2.5e-3 e，frame charge 列与 last_result().mulliken_charges 两处）、frame.scalars.energy 为有限负 eV、不支持元素 Z=99 抛 RuntimeError、device="cuda" 抛 ValueError）
- [ ] Implement bind_semiempirical.cpp（AM1Config/AM1Evaluator/AM1Result 绑定 + single_point lambda：borrow → materialize → GIL 释放计算 → write_back；bind.cpp 注册 bind_semiempirical(m)；CMakeLists.txt 源列表接线）
- [ ] Implement atomiverse/semiempirical/__init__.py and wire the top-level re-exports（AM1Config/AM1Evaluator/AM1Result/single_point 进顶层 __all__）
- [ ] Add docstring per doxygen with units（bind_semiempirical.cpp 全部新公开绑定符号：@brief/@param 带单位 [e, Hartree, kcal/mol, eV]/@return；Python 模块 docstring 同步）
- [ ] Update existing Python tests and docs to the canonical deep paths（python/tests/conftest.py、python/tests/test_md.py、python/docs/getting_started.md；顶层便捷名用例保留一组证明再导出等价）
- [ ] Write the molpy interop check and verify in the ~/work/.x86_64 venv（python/tests/test_molpy_interop.py：importorskip("molpy")，molpy 建体系 → 取共同包装的 molrs Frame → atv 跑 10 步 LJ MD，能量有限；在 venv 内实跑 `pytest python/tests` 全绿）
- [ ] Run full check + test suite（C++ 侧 cmake build + ctest 不受影响地全绿；确认 git diff 中 bind_md.cpp 零改动）

## Testing strategy

- **正常路径（pytest, python/tests/）**：映射表 16 条逐一 `atv.<顶层名> is <深层路径>` 身份成立；`import atomiverse.potentials.pair` 等深层 import 全通；既有 LJ MD smoke（conftest 参数不变）经深层路径与顶层路径各跑一次结果一致地可运行；AM1 `single_point` 于正四面体 CH₄ 收敛且能量/电荷双通道正确。
- **边界/失败路径**：`device="cuda"`/任意非 "cpu" 值 → `ValueError` 点名 CUDA 绑定不可用（MD 与 single_point 两处）；Z=99 → C++ `std::runtime_error` 干净越过 nanobind 到 Python（不崩溃）；`test_surface.py` 的死符号 grep 门与 fail-loud 契约在新布局下继续通过。
- **域校验（science.required）**：AM1 CH₄ Mulliken 电荷对钉死参考 oracle（AM1 1985 论文方法、AM1 参考实现 23.2.5 输出，C −0.262125 / H +0.065531 e）容差 2.5e-3 e——沿用 C++ 单测已引证的值，Python 侧不新增参考值；`last_result().electronic_energy_hartree`/`heat_of_formation_kcal_mol` 可读且有限。
- **互操作（验收级，venv 实跑）**：`~/work/.x86_64` venv 内 molpy（源 `~/work/molcrafts/molpy`）与 atomiverse 混合 import，molpy 建体系经 molrs Frame 通货交给 atv 跑短 MD——`pytest python/tests` 全绿即验收；无 molpy 的环境该文件 skip，不阻塞 C++ 流水线。
- **不回归**：C++ `tests/` 不新增不修改（三层铁律：ctest 侧零第三方读取）；ctest 全绿；`bind_md.cpp` 零 diff（保护 04-embed 的 lambda 位点计划）。

## Out of scope

- **AM1 梯度 / AM1-MD / AM1 几何优化**——`AM1Evaluator` 无 Force 能力，绑进 MD/Minimize 无意义；未来独立 spec。
- **CUDA Python 绑定**（`cuda::md::MD` 等）与 **TaskSpec 的 device 字段**——`device=` 参数今日只收 "cpu"，CUDA 分发沿 TaskSpec 扩展，均为后续 spec。
- **SCF / Minimize / MC / Path 的 Python 绑定**——维持 `bind.cpp` 既有注释的立场，各自成 spec。
- **`_atomiverse` 扩展模块的子模块化**（nb 子模块镜像命名空间）——扁平 `_cpp` 是私有细节，重构它对用户面零收益且干扰 04-embed。
- **`_device.py` 的 molpy device-registry 实装**——本 spec 只引用其方向佐证 device= 决策，不实现 handler。
- **mpi-launcher 链 01–05 的任何内容**——本 spec 与其唯一的交互约束是 `bind_md.cpp` 零改动、绑定 lambda 位点结构不变。
- **DeprecationWarning / 兼容 shim**——用户裁决明确不做（未发布包）。

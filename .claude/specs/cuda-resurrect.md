---
spec: cuda-resurrect
status: approved
created: 2026-07-22
chain: cuda-resurrect
phase: single
grilled: false
---

# cuda-resurrect — 复活 CUDA 后端:镜像 CPU 现有 API

## Summary

CUDA 后端在 `dev`(以及 `backend-mirror` 分支)上**根本编不过** —— GH200 节点
`-DATV_BUILD_CUDA=ON` 构建报 ~180 个错误,全是**既存腐烂**:CUDA 侧从未跟上之前
CPU 侧的三次重构。本 spec **只改 CUDA 侧**,把它对齐到 CPU 当前的公共 API,使
`cuda/` 重新在 GPU 节点上**配置 + 编译 + 通过 ctest**。**不改 CPU 代码、不改物理、
不加 CUDA 新功能** —— 纯粹是让死掉的后端复活,镜像 CPU 现在的样子。

背景/来源:`.claude/notes/notes.md` [2026-07-22] 条目 + memory `cuda-backend-rotted-on-dev`。
错误在 backend-mirror-04 的 GPU 验证中暴露,与 backend-mirror 重构正交(改名 grep 一致,
无一错误来自改名)。

## Domain basis

无新物理。CUDA kernel 的计算数学不变;只有**主机边界 API**(Frame 访问、精度类型、
单位常量)需要重新对齐到 CPU 现有表面。所有改动是机械的 API 迁移,不涉及数值算法。

## Design

三类错误,修法就是**照 CPU 现在怎么写**:

| 错误类别(~180) | 受影响 cuda 文件数 | 镜像 CPU 的修法 |
|------------------|--------------------|-----------------|
| **Frame 旧成员** `.atoms` / `.box` / `.scalars`(~150 次) | ~50 | Frame 现在是 `blocks + meta + simbox`。改用 `frame.span(kAtom_x)` 等 typed span、`frame.meta.get(kEnergy)`、`frame.simbox`。镜像 `cpu::Frame` 现有访问(见 `cpu/frame/frame.h` 与任一 cpu 消费者)。 |
| **已删精度类型** `Real` / `core/real.h`(~9 文件) | 9 | `core/real.h` 与 catch-all `Real` 已删除。主机侧用 `StorageReal`/`AccReal`(`core/precision.h`),设备侧用 cuda 自有 `real_t`/`work_t`/`prec_t`(`cuda/frame/types.cuh`)。`transfer.cu` 的 `const Real*` 主机指针 → `const StorageReal*`(或按语境 `AccReal`)。 |
| **已改名单位常量** `kFtm2v`/`kMvv2e`/`kBoltzmann`/`kPI`/`kCoulomb`(~19 文件) | 19 | CPU 已把 LAMMPS 常量改成描述性名:`kFtm2v`→`units::force_to_accel`、`kMvv2e`→`units::mass_v2_to_energy`、`kBoltzmann`→`units::` 玻尔兹曼常量、`kPI`→`std::numbers`/`M_PI` 或 units 内常量、`kCoulomb`→ 静电常量。对齐 `core/units.h` 现有名(`kAngstromInBohr`/`kEvInHartree`/`kB_ev_per_k`/…);MD 转换见 `cpu/dynamics/velocity_verlet.h`、`fix_langevin.h`(`units::force_to_accel` 即 ftm2v 等价)。 |

**依赖顺序(编译前后依赖):** 精度/单位(叶子,被广泛 include)→ frame 边界(`transfer`/`context`,cuda 读 cpu::Frame 的入口)→ 各 domain `.cu`(dynamics/potentials/observe/record/scf/planewave/longrange/minimize/neighbor/evaluate)。CUDA 是**全有或全无编译**:改完全部才能编过,所以本 spec 是**单个 spec 分阶段任务**,而非可中间绿的 chain。

**验证只能在 GPU 节点。** 登录节点无 CUDA。用 `sbatch -p gpu -A naiss2026-4-715-gpu --gpus 1`
+ memory `cuda-gpu-node-build-recipe` 的配方(aarch64 Rust toolchain override + cache
wipe + arch-split cargo PATH + BLAS hints)。每轮 ~6-15 min。

**镜像纪律:** CUDA 的每个改动都应能在 CPU 侧找到对应写法。若某处 CPU 无对应(cuda 独有
kernel 逻辑),保持 kernel 数学不变,只迁移它读写 Frame/精度/单位的边界。

### Reuse decision

| Candidate | Decision |
|-----------|----------|
| CPU `cpu/frame/frame.h` span/meta/simbox API | **pattern** — cuda 镜像其访问方式 |
| `core/precision.h` 精度分层 + `cuda/frame/types.cuh` `real_t` | **reuse** — 主机 StorageReal/AccReal,设备 real_t/work_t/prec_t |
| `core/units.h` 现有常量 + `units::force_to_accel` 等 | **reuse** — cuda 对齐现名 |
| CUDA kernel 计算体 | **reuse** — 数学不变,只改边界 |

## Files to create or modify

**精度/单位(叶子,先修):**
- `cuda/potentials/pair_kernel.h`, `cuda/potentials/coul_long.cuh`,
  `cuda/planewave/internal/{device_solver,fft,xc_lda}.cuh`,
  `src/cuda/planewave/{internal/device_solver,pw_scf}.cu`,
  `src/cuda/scf/{cublas_dispatch.cuh,device_integrals.cu}` —— `Real`/`core/real.h` → 精度分层
- `cuda/dynamics/{init,kinetic_energy,nve_limit}.cuh`, `cuda/minimize/sd.cuh`,
  `cuda/observe/stress.cuh` + 对应 `.cu`,`src/cuda/{communication/multi_gpu,
  dynamics/*,longrange/{ewald,pme},minimize/sd,observe/stress,potentials/coul_long}.cu`
  —— 单位常量对齐

**Frame 边界:**
- `cuda/frame/transfer.cuh` + `src/cuda/frame/transfer.cu`, `cuda/component/context.cuh`
  —— `.atoms/.box/.scalars` → blocks/meta/simbox

**Domains(bulk,~45 文件):**
- `cuda/{evaluate,dynamics,minimize,observe,record,scf,potentials}/*.cuh` +
  `src/cuda/{driver,dynamics,minimize,observe,potentials,scf,neighbor,communication,
  longrange}/*.cu` —— Frame API + 单位对齐

(精确清单见 `git grep` 于本 spec 的 grep 模式;impl 时逐子系统扫。)

## Tasks

- [ ] 阶段 A(精度/单位):把 cuda 树里 `Real`/`core/real.h` 改为主机 `StorageReal`/`AccReal` + 设备 `real_t/work_t/prec_t`;把 `kFtm2v/kMvv2e/kBoltzmann/kPI/kCoulomb` 改为 `core/units.h` 现名(`units::force_to_accel` 等)。逐子系统 grep 清零。
- [ ] 阶段 B(frame 边界):`transfer.cu/.cuh` + `component/context.cuh` 改用 `span/meta/simbox`,镜像 `cpu::materialize`/`to_cuda` 的现有 Frame 访问。
- [ ] 阶段 C(domains):逐 cuda 子系统(dynamics → potentials → neighbor → observe → record → minimize → evaluate → scf → planewave → longrange)把 `.atoms/.box/.scalars` 迁到新 Frame API。每个子系统对照同名 CPU 实现镜像。
- [ ] 阶段 D(GPU 绿):在 GH200 节点 `sbatch` 完整 `-DATV_BUILD_CUDA=ON -DATV_BUILD_TESTS=ON`(+HF/KSPACE)configure+build+ctest,迭代到全绿。harmless `[CUDA] Kernel ... sm20 ...` warning 忽略。
- [ ] 加 `tests/architecture/cuda_no_stale_api_scrub.sh`:grep cuda 树无 `\.atoms/\.box/\.scalars` Frame 访问、无 `\bReal\b`/`core/real.h`、无旧单位常量名。注册 ctest(纯 grep,登录节点也能跑)。
- [ ] 运行:登录节点 CPU 套件仍 236/236 绿(证明没碰坏 CPU);GPU 节点 cuda ctest 绿。

## Testing strategy

- **Grep 门禁(登录节点可跑):** `cuda_no_stale_api_scrub.sh` 锁死"cuda 不再引用已删的 Frame 成员/精度类型/单位常量"。进 ctest。
- **GPU 节点构建+ctest:** 唯一能验证 CUDA 真编过 + kernel 跑对的地方。这是本 spec 的核心验收,不进登录节点 ctest(无 CUDA)。
- **CPU 回归:** 登录节点 236/236 保持绿,证明零 CPU 改动。
- **无新 e2e / 无第三方 oracle。** kernel 正确性由既有 cuda 单测(GPU 节点)覆盖。
- **Regression:** `regressions/layout/cuda_mirrors_cpu/` —— 可选,记录 cuda 与 cpu 目录镜像(布局层,登录节点可编的 host 部分)。

## Out of scope

- 任何 CPU 代码改动(纯 CUDA 侧)。
- 命名空间对齐目录(backend-mirror 已推迟的独立后续)。
- CUDA 新功能 / 新 kernel / SCF 梯度(CUDA SCF 仍能量/密度 only)。
- Torch 后端。
- 精度默认值调整(保持 `ATV_CUDA_REAL=float` 等现有默认)。

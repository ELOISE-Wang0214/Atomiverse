# Specs index

- [scf-driver-refactor-01-plumb](scf-driver-refactor-01-plumb.md) — 门面管道铺设：scf::Phase(fock/mix/diagonalize/density/check) + PhaseKind::SCF + PhaseTraits + LoopMetrics 收敛字段 + ctx.scf 切片；删除 scf_phase_scrub.sh 门 [approved]
- [scf-driver-refactor-02-stages](scf-driver-refactor-02-stages.md) — 阶段组件化：五 stage 组件（构造器直传参数、无 config）+ ScfWorkspace + ScfConverged handler；旧循环仍在位 [approved]
- [scf-driver-refactor-03a-convert](scf-driver-refactor-03a-convert.md) — 方法驱动化：MeanFieldMethod<D> 派生 CPUDriver<scf::Phase>、drive_scf_ 私有漏斗、按方法拆文件；旧 SCFSolver 留树作逐位对拍活参考 [approved]
- [scf-driver-refactor-03b-retire](scf-driver-refactor-03b-retire.md) — 退役旧机器：删 SCFSolver/SCFConfig/make_config、梯度解耦、首个文献值回归 regressions/scf/h2o_rhf_sto3g [approved]
- [scf-driver-refactor-04a-nest](scf-driver-refactor-04a-nest.md) — 原生 Driver 节点嵌套：DriverNode 擦除基 + PendingBinding 变体 + on(phase).add(driver) 重载；删 sub_driver.h [approved]
- [scf-driver-refactor-04b-wire](scf-driver-refactor-04b-wire.md) — 外层接线：HF 纯多态 tag、删 build/launch（run/launch 债务终结）、MeanFieldOutputs 收窄、scf_evaluator 消亡 [approved]
- [scf-driver-refactor-05-rename](scf-driver-refactor-05-rename.md) — 编译开关改名：ATV_BUILD_RHF→ATV_BUILD_HF / ATV_HAS_RHF→ATV_HAS_HF，两后端 ~67 文件机械替换 [approved]
- [scf-driver-refactor-06-criterion-unify](scf-driver-refactor-06-criterion-unify.md) — 终止即组件：三内置判据迁 Criterion、删 LoopHandler/add_handler、统一 verdict 轮询 [approved]

Chain scf-driver-refactor complete (8 links, all approved + grilled): 01-plumb → 02-stages → 03a-convert → 03b-retire → 04a-nest → 04b-wire → 05-rename → 06-criterion-unify. Implement in order via /mol:impl-all scf-driver-refactor.

Completed work lives in the git log (campaigns: frame-block, cpu-breadth,
pw-breadth, se-breadth / RM1, python-on-demand-install).

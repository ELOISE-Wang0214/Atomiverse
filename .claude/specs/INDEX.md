# Specs index

- [scf-driver-refactor-04a-nest](scf-driver-refactor-04a-nest.md) — 原生 Driver 节点嵌套：DriverNode 擦除基 + PendingBinding 变体 + on(phase).add(driver) 重载；删 sub_driver.h [approved]
- [scf-driver-refactor-04b-wire](scf-driver-refactor-04b-wire.md) — 外层接线：HF 纯多态 tag、删 build/launch（run/launch 债务终结）、MeanFieldOutputs 收窄、scf_evaluator 消亡 [approved]
- [scf-driver-refactor-05-rename](scf-driver-refactor-05-rename.md) — 编译开关改名：ATV_BUILD_RHF→ATV_BUILD_HF / ATV_HAS_RHF→ATV_HAS_HF，两后端 ~67 文件机械替换 [approved]
- [scf-driver-refactor-06-criterion-unify](scf-driver-refactor-06-criterion-unify.md) — 终止即组件：三内置判据迁 Criterion、删 LoopHandler/add_handler、统一 verdict 轮询 [approved]


Completed work lives in the git log (campaigns: frame-block, cpu-breadth,
pw-breadth, se-breadth / RM1, python-on-demand-install).

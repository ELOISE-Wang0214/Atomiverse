# Specs index

- [scf-driver-refactor-06-criterion-unify](scf-driver-refactor-06-criterion-unify.md) — 终止即组件：三内置判据迁 Criterion、删 LoopHandler/add_handler、统一 verdict 轮询 [approved]
- [scf-slater-determinant-01-type](scf-slater-determinant-01-type.md) — 新公共值类型 cpu::scf::SlaterDeterminant：按自旋道存 C/ε/f、density 委托 linalg、validate(S) 不变量、with_occupation ΔSCF 入口；零消费者接线 [approved]
- [scf-slater-determinant-02-rehome](scf-slater-determinant-02-rehome.md) — 疏散 scf_result.h 同住者：三枚举→scf_enums.h、SCFInitialGuess→scf_guess.h、GradientResult→gradient_result.h、IterationRecord 冗余删除；行为零变化 [approved]
- [scf-slater-determinant-03-unbag](scf-slater-determinant-03-unbag.md) — SCFResult 消亡：方法持有 SlaterDeterminant、run() 归 void、消费面全迁（前置 01/02 + driver-06） [approved]
- [scf-slater-determinant-04-methods](scf-slater-determinant-04-methods.md) — cpu/scf 13 个 build_* 收编为成员方法/构造函数（kernel 逐字搬移零性能回归）+ no_factory_scrub build_* 通配门禁（scf 域）与旧名全树黑名单 [approved]


Completed work lives in the git log (campaigns: frame-block, cpu-breadth,
pw-breadth, se-breadth / RM1, python-on-demand-install).

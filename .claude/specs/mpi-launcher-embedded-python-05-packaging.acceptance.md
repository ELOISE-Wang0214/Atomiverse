---
slug: mpi-launcher-embedded-python-05-packaging
criteria:
  - id: ac-001
    summary: pyproject 定义 ATV_BUILD_LAUNCHER env-override，默认 OFF
    type: code
    pass_when: |
      python/pyproject.toml 的 [tool.scikit-build.cmake.define] 含
      ATV_BUILD_LAUNCHER = {env = "ATV_BUILD_LAUNCHER", default = "OFF"}。
    status: pending
  - id: ac-002
    summary: launcher install 规则 SKBUILD/BINDIR 双路且 RPATH 就位
    type: code
    pass_when: |
      app/atomiverse/CMakeLists.txt 中 launcher target 的 install(TARGETS ... RUNTIME
      DESTINATION ...) 在 SKBUILD 下指向 ${SKBUILD_SCRIPTS_DIR}、否则指向
      ${CMAKE_INSTALL_BINDIR}，且 target 带 INSTALL_RPATH "$ORIGIN/../lib"；
      不存在任何 console-script / [project.scripts] 入口指向 launcher。
    status: pending
  - id: ac-003
    summary: 纯 CMake install 覆盖库与头文件且 wheel 构建不携带
    type: code
    pass_when: |
      根 CMakeLists.txt 在 NOT SKBUILD 守卫下 install atomiverse 库与
      include/atomiverse 头目录；ATV_BUILD_LAUNCHER=ON 且 MPI 未找到时
      configure 以 FATAL_ERROR 终止并点名两个开关。
    status: pending
  - id: ac-004
    summary: InterpreterLocator 单元测试通过且零解释器执行
    type: code
    pass_when: |
      test_interpreter_locator 注册进 ctest 并通过（秒级）：override 优先、
      PATH 搜索、缺失 throw、版本不匹配 throw；测试源与
      interpreter_locator.{h,cpp} 均不含 Python.h/popen/exec 调用，
      测试期无任何第三方包读取。
    status: pending
  - id: ac-005
    summary: 无静默回退——发现失败与 ABI 不匹配均清晰致命
    type: code
    pass_when: |
      python_host.cpp / interpreter_locator.cpp 中不存在编译期烘焙的
      解释器/前缀路径作为运行时回退；缺 python3 与 minor 版本不匹配两条
      路径均产生含 ATV_PYTHON 提示（后者含两个版本号）的错误并非零退出。
    status: pending
  - id: ac-006
    summary: fresh venv pip 安装后 PATH 上出现真实 ELF 原生 atomiverse
    type: runtime
    evaluator_hint: 手动验收；需网络与全新 venv；不得进 ctest
    pass_when: |
      全新 venv 内 ATV_BUILD_LAUNCHER=ON pip install . 成功；wheel 解包可见
      *.data/scripts/atomiverse；安装后 command -v atomiverse 落在 venv/bin
      且 file 判定为 ELF 可执行（非文本脚本）。
    status: pending
  - id: ac-007
    summary: atomiverse script.py 单 rank 使用 venv 的 site-packages
    type: runtime
    evaluator_hint: 手动验收；同 ac-006 的 venv
    pass_when: |
      atomiverse script.py（单 rank）执行的脚本打印的 sys.prefix 等于该 venv
      前缀且成功 import atomiverse（模块路径位于 venv 的 site-packages）。
    status: pending
  - id: ac-008
    summary: cmake --install 落出可用布局
    type: runtime
    evaluator_hint: 本机 scratch 前缀即可，无需网络
    pass_when: |
      ATV_BUILD_LAUNCHER=ON 配置后 cmake --install build --prefix <scratch>
      产生 bin/atomiverse（ELF）、lib/ 下 atomiverse 库、include/atomiverse/ 头树。
    status: pending
  - id: ac-009
    summary: user-guide 启动器文档齐备并入索引
    type: docs
    pass_when: |
      docs/user-guide/launcher.md 覆盖 python-vs-mpirun/srun 对照、环境要求
      （python3 minor 匹配、libpython/LD_LIBRARY_PATH、Lmod）、ATV_PYTHON 覆盖、
      两条安装路径与 wheel 不可移植说明，且 docs/user-guide/index.md 链接到它。
    status: pending
---

# Acceptance criteria

ac-001–ac-005 为代码/单元级，`/mol:impl` 本地闭环可验；ac-006/ac-007 需要网络与全新 venv，ac-008 需要一次带前缀的 install 运行——三者为验收时手动 runtime 判据，**不得**以任何形式注册进 ctest（tests/ 单元-秒级-零第三方读取是硬规则）；ac-009 为文档判据。

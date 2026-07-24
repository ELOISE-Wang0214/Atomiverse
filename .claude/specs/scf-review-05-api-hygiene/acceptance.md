---
slug: scf-review-05-api-hygiene
criteria:
  - id: ac-001
    summary: MolecularGrid is constructor-built; no build_molecular_grid symbol
    type: code
    pass_when: |
      grid.h declares `MolecularGrid(const native::Atom* atoms, int natm,
      GridLevel level, int n_radial = 50, int n_angular = 110)` as a member
      constructor; `grep -rn 'build_molecular_grid' include/ src/ tests/`
      returns no declaration/definition/call (only historical prose comments,
      if any).
    status: pending
  - id: ac-002
    summary: no_factory_scrub whitelist removed and gate RESTORED TO GREEN
    type: code
    evaluator_hint: "bash tests/architecture/no_factory_scrub.sh <repo-root> (red on dev today)"
    pass_when: |
      tests/architecture/no_factory_scrub.sh contains no `build_molecular_grid`
      token; running it against the repo root exits 0. Note: the scrub is RED on
      dev pre-chain (build_se_part, adc build_m*_ x3, build_imds x3); this spec
      clears the gf2/adc renames and the grid factory, spec 03 clears build_imds
      — with 03 landed before or alongside, exit must be 0.
    status: pending
  - id: ac-003
    summary: compatibility_surface_scrub refreshed and RESTORED TO GREEN
    type: code
    evaluator_hint: "bash tests/cpu/compatibility_surface_scrub.sh <repo-root> (red on dev today)"
    pass_when: |
      tests/cpu/compatibility_surface_scrub.sh exits 0; its grid guard matches
      the MolecularGrid constructor and fails if `natm` is immediately followed
      by `int n_radial` without an explicit GridLevel (implicit-SG0 intent
      preserved), no longer references build_molecular_grid, and the 6
      pre-existing `using molrs::Element` violations (gpw tests/regressions/docs)
      are fixed to qualified molrs::Element usage.
    status: pending
  - id: ac-004
    summary: build_se_part and adc build_m*_ renamed to non-factory verbs
    type: code
    pass_when: |
      `grep -rn 'build_se_part\|build_mss_\|build_mab_ea_\|build_mij_ip_' src/`
      returns nothing; src/cpu/scf/gf2.cpp defines `compute_se_part(...)` (2 call
      sites) and src/cpu/scf/adc.cpp defines compute_mss_/compute_mab_ea_/
      compute_mij_ip_ with all call sites updated.
    status: pending
  - id: ac-005
    summary: rys_roots throws instead of aborting
    type: code
    pass_when: |
      src/cpu/scf/rys_roots.cpp function rys_roots() contains no std::abort()
      and no std::fprintf; the unsupported-nroots and nonzero-err branches
      `throw std::invalid_argument`/`std::runtime_error` carrying the original
      messages.
    status: pending
  - id: ac-006
    summary: unit test pins the rys_roots throw for nroots>7
    type: code
    evaluator_hint: "marker: RysRoots"
    pass_when: |
      tests/cpu/scf/test_rys_roots.cpp asserts EXPECT_THROW/ASSERT_THROW that
      native::rys_roots(8, x, u, w) throws std::invalid_argument; the
      test_cpu_scf_rys_roots binary passes.
    status: pending
  - id: ac-007
    summary: local pi/M_PI constants swept to units::pi
    type: code
    pass_when: |
      In the enumerated cpu/scf files, `grep -nE 'constexpr double (kPi|PI)\s*='`
      and `grep -n 'M_PI'` return nothing; each edited file includes
      atomiverse/core/units.h and uses units::pi. Numerics bitwise-unchanged
      (units::pi == 3.14159265358979323846 == prior literals).
    status: pending
  - id: ac-008
    summary: 11 internal .cpp relocated under scf/internal and CMake repathed
    type: code
    pass_when: |
      cphf_pack, ucphf_pack, orbital_hessian, stability, ecp_integrals, minao,
      one_electron, shell_pair, integrals, integral_derivatives, gradient exist
      only under src/cpu/scf/internal/ (none at src/cpu/scf/ root); CMakeLists.txt
      lists all 11 with the internal/ prefix; the target configures and compiles.
    status: pending
  - id: ac-009
    summary: architecture-witness regression builds, links, runs, exits 0
    type: runtime
    evaluator_hint: "regressions target: reg_scf_api_hygiene"
    pass_when: |
      regressions/layout/scf_api_hygiene/reg_scf_api_hygiene.cpp (registered via
      atv_add_regression(layout scf_api_hygiene)) compiles as an out-of-library
      consumer, constructs a MolecularGrid via the new constructor, catches the
      std::invalid_argument thrown by native::rys_roots(8, 1.0, u, w), and
      returns 0. No third-party runtime, no numeric golden.
    status: pending
  - id: ac-010
    summary: full unit suite and existing SCF/DFT regressions unchanged
    type: runtime
    pass_when: |
      ctest full unit suite passes; the regressions target passes including
      existing regressions/scf/* cases reproducing their cited literature values
      unchanged (behavior-preserving refactor confirmed).
    status: pending
---

# Acceptance criteria

- ac-001..ac-003 — item 10 网格工厂转构造器 + 门禁 lockstep(白名单删除、隐式-SG0 守卫刷新)。
- ac-004 — item 11 `build_se_part`→`compute_se_part`(同时保证 `no_factory_scrub` 的 `cpu/scf` `build_*` 通配门禁绿色)。
- ac-005/ac-006 — item 14 `abort`→`throw` + 新单元测试(唯一行为改动:更好的错误路径,可捕获、可测)。
- ac-007 — item 9 π 常量归一到 `units::pi`,逐位等价。
- ac-008 — item 13 11 个内部 `.cpp` 归位 + CMake 显式清单 11 行重定向。
- ac-009 — 架构见证型回归(遵 `reg_architecture_gates` 先例):库外见证构造器可链接、异常可捕获。
- ac-010 — 全量单元套件 + 既有 SCF/DFT 回归保持不变,是行为保持的端到端护栏。

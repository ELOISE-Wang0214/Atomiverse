---
spec: arch-refactor-06-frame
created: 2026-07-10
criteria:
  - id: ac-001
    summary: "cpu::Frame no longer declares Bonds/Angles/Dihedrals shims and no code uses them"
    type: code
    pass_when: "include/atomiverse/cpu/frame/frame.h contains no `struct Bonds`, `struct Angles`, or `struct Dihedrals` declaration and no `bonds` / `angles` / `dihedrals` data member; ripgrep for `\\.bonds\\.|\\.angles\\.|\\.dihedrals\\.` over tests/cpu/, benchmarks/, and tools/ returns no match on a cpu::Frame instance (cuda::Frame / torch::Frame member accesses in their own backend tests excluded)."
    status: pending
  - id: ac-002
    summary: "Migrated fixtures build and pass on the per-atom topology API"
    type: runtime
    evaluator_hint: "ctest -R test_cpu_md_no_alloc; bench binary bench_nve_limit_lammps_align"
    pass_when: "With the shims deleted, test_cpu_md_no_alloc passes, and all bench_nve_limit_lammps_align cases pass under -DATV_BUILD_BENCHMARKS=ON — including NveLimitLammpsAlign_VelocityVerletConservesEnergy, whose nonzero initial-energy guard and 1e-4 conservation bound now hold with real bond forces from the per-atom layout."
    status: pending
  - id: ac-003
    summary: "Golden round-trip writes image and tag back per the contract"
    type: runtime
    evaluator_hint: "ctest -R test_cpu_materialize"
    pass_when: "A new test in tests/cpu/test_materialize.cpp does materialize → mutate (x/v/f/charge, image entries, permuted tags, scalars) → write_back → re-read and passes: molrs atoms columns ix/iy/iz equal the mutated int8 image values (widened to i32) and the atoms.tag column equals the permuted TagInt values (narrowed to i32), alongside the previously covered x/v/f/charge/scalars fields."
    status: pending
  - id: ac-004
    summary: "type/mass/topology/box are provably not written back"
    type: runtime
    evaluator_hint: "ctest -R test_cpu_materialize"
    pass_when: "The same golden test mutates frame.atoms.type and frame.atoms.mass on the cpu::Frame before write_back and then asserts the molrs atoms.type / atoms.mass columns, the bonds/angles/dihedrals blocks, and the simbox still hold their pre-materialize values; and a tag value > INT32_MAX makes write_back throw std::runtime_error naming the tag field (no silent truncation)."
    status: pending
  - id: ac-005
    summary: "materialize.h documents the full write_back field contract"
    type: code
    pass_when: "include/atomiverse/cpu/frame/materialize.h's doxygen comments enumerate, for every field in {x/y/z, vx/vy/vz, fx/fy/fz, charge, image, tag, type, mass, topology, box, scalars/stress/stress_kinetic}, whether write_back writes it and why the non-written ones are designated per-run static/immutable; state that materialize never reads the scalars block (write-only channel); and name write_frame_xyz_meta as the per-frame metadata path kept outside write_back."
    status: pending
  - id: ac-006
    summary: "molrs.hpp comment matches the single always-f64 bridge reality"
    type: code
    pass_when: "include/atomiverse/molrs.hpp no longer mentions bridge_f32/bridge_f64 or configure-time precision selection; its comment states molrs-core hardcodes f64 with a single always-double bridge, consistent with CMakeLists.txt:658-659."
    status: pending
  - id: ac-007
    summary: "CUDA dual-topology and contract inheritance are documented at the source"
    type: code
    pass_when: "include/atomiverse/cuda/frame/frame.cuh's header comment describes both topology representations (single-GPU flat BondList/AngleList/DihedralList; multi-GPU per-atom SoA with global-tag partners) and qualifies the 'field set mirrors cpu::Frame' claim accordingly; src/cuda/frame/transfer.cu's header comment states kernels are not yet ported to per-atom topology so cpu::to_cuda re-flattens every transfer; include/atomiverse/cuda/frame/materialize.cuh notes cuda::write_back inherits the cpu::write_back contract via to_cpu and that to_cpu's returned field set (x/v/f/tag/charge/image/scalars/box) covers every contract-written field."
    status: pending
  - id: ac-008
    summary: "Stale frame.bonds doc references corrected tree-wide"
    type: code
    pass_when: "ripgrep for `frame\\.bonds|frame\\.angles|frame\\.dihedrals` over include/atomiverse/cpu/, src/cpu/, and docs/ returns no remaining reference that describes current cpu::Frame storage (halo_exchange.h/.cpp, domain.cpp, pair_lj_cut.h, docs/developer/cpu-internals.md, docs/tutorial/ch10-classical-md.md all reworded to the per-atom frame.atoms.num_bond/bond_atom layout)."
    status: pending
  - id: ac-009
    summary: "Full suite green and PE stage-5 benchmark unperturbed"
    type: runtime
    evaluator_hint: "ctest --test-dir build; benchmark binary polyethylene_tension stage 5"
    pass_when: "The full default ctest suite passes after the frame.h shim deletion, and a polyethylene_tension run's stage-5 throughput (Matom-step/s) lies within ±5% of the most recent baseline CSV under benchmarks/md/pe/performance/runs/."
    status: pending
out_of_scope:
  - "CUDA kernel 迁移到 per-atom 拓扑（仅文档化现状）"
  - "精度 tier 变更（precision-migration spec）"
  - "box/topology 回写支持（钦定 per-run 静态）"
  - "torch::Frame 拓扑表示"
---

# Acceptance — arch-refactor-06-frame

ac-001/ac-002 覆盖 shim 删除与 fixture 迁移；ac-003/ac-004 是 write_back 契约的 golden
round-trip 正反两面（写字段往返、非写字段不动 + tag 溢出 fail-loud）；ac-005–ac-008 固化
文档/注释修正（契约文档、molrs.hpp、CUDA 双拓扑、陈旧引用清扫）；ac-009 是链级热路径门禁
（frame.h 为热路径 include，全量套件 + PE stage-5 无扰动）。CUDA 侧行为由契约继承覆盖
（ac-007 文档 + 既有 test_cuda_materialize 于 GPU 节点回归，见 spec Testing strategy——
GPU 节点执行不设为绑定判据，因登录节点无 CUDA toolchain）。

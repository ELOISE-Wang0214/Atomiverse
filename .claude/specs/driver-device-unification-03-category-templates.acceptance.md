---
slug: driver-device-unification-03-category-templates
criteria:
  - id: ac-001
    summary: Four device-generic category templates exist with one run()/validate each
    type: code
    pass_when: |
      `ATV::EvaluateDriver<DT>`, `ATV::DynamicsDriver<DT>`,
      `ATV::OptimizeDriver<DT>`, `ATV::ResponseDriver<DT>` are defined in
      `facade/driver_categories.h` (+ `_impl.h`), each with exactly one `run(...)`
      body and one `validate_components()` body. Dynamics adds set_integrator/
      set_thermostat/set_barostat; each template carries a `molrs::FrameRef`
      overload.
    status: pending
  - id: ac-002
    summary: The 7 hand-written category headers are deleted
    type: code
    pass_when: |
      None of cpu/{evaluate,dynamics,optimize,response}_driver.h or
      cuda/{evaluate,dynamics,optimize}_driver.cuh exist. `grep -rn` finds no
      `#include` of any of them.
    status: pending
  - id: ac-003
    summary: Concrete drivers reparented onto the templates
    type: code
    pass_when: |
      `cpu::md::MD : ATV::DynamicsDriver<cpu::device_tag>`,
      `cuda::md::MD : ATV::DynamicsDriver<cuda::device_tag>`; Evaluate/Minimize
      reparented on both backends; `cpu::vib::Hessian : ATV::ResponseDriver<cpu::device_tag>`.
    status: pending
  - id: ac-004
    summary: Old Backend public staging names removed; stage_* is the only public surface
    type: code
    pass_when: |
      `set_frame`, `materialize`, `retrieve`, `drain_analysis`, `sync` no longer
      appear as public methods on cpu::Backend / cuda::Backend (private impl may
      keep them). All driver call sites use stage_in/stage_out/drain.
    status: pending
  - id: ac-005
    summary: Category-template test + full regression pass; continuation preserved
    type: runtime
    pass_when: |
      `ctest -R test_category_templates` exits 0 (Dynamics NVE + cumulative
      step counter across run(5);run(5)). Full `ctest --test-dir build
      --output-on-failure` exits 0 — test_md, cuda_md*, minimizer, evaluate,
      vib/hessian, and hot-path no-sync guards all green (unified teardown adds
      no per-step sync; NHC/Langevin extended state persists across run calls).
    status: pending
---

---
slug: driver-api-02-aimd-wrappers
status: in-progress
created: 2026-07-27
revised: 2026-07-27
domain: dynamics / electronic
queue: driver-api
deferred: false
depends_on: [driver-api-01-evaluator-registration]
unblocks: []
grilled: true
worktree: grok-worktrees/driver-api-02-aimd-wrappers
branch: wf/driver-api-02-aimd-wrappers
---

# Spec: `driver-api-02-aimd-wrappers` — CPMD / XL / SGCP are types, not flags

**No backward compatibility.** Worktree: `wf/driver-api-02-aimd-wrappers`.

## Summary

Remove **public** AIMD mode switches on single-point methods. Dynamics strategy is a **concrete evaluator type** registered with `set_evaluator`.

| Public type | Owns | Role |
| ----------- | ---- | ---- |
| bare `GPWEvaluator` / `RHF`/`RKS`/… | — | BOMD / single-point only |
| `cpu::electronic::CarParrinello` | `GPWEvaluator` | fictitious-mass CP + RATTLE |
| `cpu::electronic::XlBomdForce` | `GPWEvaluator` | Niklasson XL-BOMD |
| `cpu::electronic::SecondGenCp` | `GPWEvaluator` | ASPC + partial SCF |
| `cpu::electronic::CarParrinelloMethod<M>` | `M` (RHF/RKS/…) | Gaussian AO-CPMD |
| `cpu::electronic::XlBomdMethod<M>` | `M` | Gaussian XL-BOMD |

```cpp
sim.set_integrator(VelocityVerlet{dt});
sim.set_evaluator(electronic::CarParrinello{std::move(gpw), /*mu=*/2000.0});
sim.set_observer(md::CpmdConserved{});  // or on(end_of_step).add
```

## Design

### GPWEvaluator

- AIMD members (`cpmd`, `cpmd_mu`, `xl_bomd`, `second_gen_cp`, …) move to **private**.
- Friend: `electronic::CarParrinello`, `XlBomdForce`, `SecondGenCp`.
- Bare evaluator: `timestep_required() == false`; no CP without wrapper.
- Wrappers: inherit `cpu::GPW` + `TimestepConsumer` (CP only requires dt), own `GPWEvaluator`, set private flags in **ctor only**, forward `build`/`launch`/accessors.

### MeanFieldMethod

- AIMD members become **protected** (not public).
- `CarParrinelloMethod<M>` / `XlBomdMethod<M>` publicly inherit `M`, set flags in ctor.

### Registration

`set_evaluator(std::move(wrapper))` — kind remains `GPW` or nested `HF`.

### Mutual exclusion

Type-level: one wrapper type per object. Internal guards for illegal method combos (OT+CP, …) stay in `build()`.

## Tasks

1. GPW private flags + three wrapper headers.
2. MFM protected flags + method templates.
3. Migrate unit tests + regressions (`cpmd-*`, `xl-bomd-*`, `second-gen-*`).
4. Build `test_cpu_gpw_cpmd` (+ related) green.

## Out of scope

- QM/MM × CPMD; multi-k UKS CP productization beyond current gates.
- Extracting CP launch logic out of `gpw_evaluator.cpp` into separate TU (flags-out is enough for API).

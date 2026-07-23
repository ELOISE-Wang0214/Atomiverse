# PySCF method loop — staff model

## Why does `agent_budget` exist?

It is a **Grok host / platform** control, not a science design choice.

| | |
|--|--|
| **What** | Cap on how many **LLM sessions** one workflow *run* may start (`agent()` + each `parallel()` slot), including `resume_from` turns |
| **Why** | Each session costs tokens, time, and sandbox resources. Without a cap, a buggy loop can spawn unbounded children |
| **Limits** | Default **128**, hard max **1024** per run |
| **What it is not** | Not “how many employees you have”, not “how many PySCF features exist” |

### Employee metaphor (accurate)

| Metaphor | Reality |
| -------- | ------- |
| **Headcount** | `staff_count` / `parallelism` (default **5**) — fixed crew |
| **Employees stay on staff** | Same `agent_id` continued with **`resume_from`** each wave |
| **This shift’s payroll** | `agent_budget` — how many *turns* the crew may take before the run ends |
| **Next day / next shift** | Scheduler starts a **new** run and hires again (`agent_id` does not survive process/run boundaries) |

So yes: **think of 5 employees repeating the job**.  
Budget is only **how long this shift can pay them**, not how many people you hired.

```text
staff_count = 5          →  five people
each wave                →  each person gets one feature (resume)
agent_budget             →  total person-turns this shift
                           e.g. 256 turns ≈ many waves of 5 before shift ends
```

## Flow

```text
census (once)
hire crew of 5 on first work wave
loop:
  foreman assigns ≤5 pending features
  parallel: staff[i] resume_from → pin+impl+verify for assignment[i]
  close queue/roadmap
until queue empty | payroll (budget) empty | max_waves
```

## Launch (recommended)

```text
agent_budget: 256          # shift payroll — NOT headcount
args:
  staff_count: 5           # or parallelism: 5
  max_waves: 20
  do_perf: false
  mode: exhaust
```

## Queue

`docs/developer/pyscf-method-queue.md` — survives shifts.

## Run log

(autonomous closes append lines here)

- `h2-b3lyp-hybrid-tda` | Hybrid TDDFT-TDA (B3LYP exact-exchange kernel) | S1=25.928364822972 eV (PySCF 2.14.0) | perf_status=n_a | status=closed
- `ump2-o2-sto3g` | UMP2/O2/STO-3G (triplet, cart, conv_tol=1e-12) | e_corr=-0.096393144130 Ha, e_tot=-147.730339929580 Ha (UHF e_hf=-147.633946785450 Ha; PySCF 2.14.0; tol 1e-6 Ha) | perf_status=ok | status=closed

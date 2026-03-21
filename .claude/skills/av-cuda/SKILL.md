---
name: av-cuda
description: CUDA execution model audit -- classifies components, flags host violations, validates sync discipline.
argument-hint: "[path or module]"
user-invocable: true
---

# CUDA Execution Model Audit

Delegate to the `av-cuda` agent with the user's argument (or all CUDA backend files if none given).

The agent will:
1. Read CLAUDE.md "CUDA execution model" section for the three component categories
2. Discover all CUDA components via glob (`src/cuda/**/*.cu`, `include/atomiverse/cuda/**/*.cuh`)
3. Classify each component: compute (device-resident), CPU bridge, or host-dependent
4. Audit all synchronization points for justification
5. Check for undocumented host dependencies not in the CLAUDE.md exceptions table
6. Report findings in severity-sorted `[SEVERITY] file:line` format with classification table

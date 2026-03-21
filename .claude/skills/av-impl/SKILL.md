---
name: av-impl
description: Implementation workflow for Atomiverse -- scales from quick additions to full multi-phase pipelines.
argument-hint: "<feature description or path to spec file>"
user-invocable: true
---

# Atomiverse Implementation

Read CLAUDE.md for architecture rules, conventions, and build commands.

## Step 1 — Assess Scope

Classify the task:

- **SMALL** (new component following existing pattern, < 3 files): Skip to Step 3. This includes standalone CUDA kernel or CPU component development.
- **MEDIUM** (new component with new patterns, 3-8 files): Start at Step 2.
- **LARGE** (new driver type, new Frame fields, architectural change): Full pipeline from Step 2.

State the scope classification and get user confirmation before proceeding.

## Step 2 — Domain & Spec (MEDIUM/LARGE only)

- If the feature involves equations or physics: invoke `/av-litrev` to verify scientific basis. Abort if no credible basis found for domain-critical work.
- If input is a file path, read it as the spec.
- Otherwise, invoke `/av-spec` to produce a structured specification.
- Confirm spec covers: component design, phase bindings, lifecycle needs, kernel design.

## Step 3 — Architecture Check

Scan the codebase to find the closest existing component as a pattern:
```bash
grep -rl "ComponentBase" include/atomiverse/
```

For the new code, verify:
- Which layer? (Facade header, CUDA backend .cu, CPU backend .cpp, or both)
- Facade header CUDA-free and xtensor-free?
- Component lifecycle: needs build()/destroy()?
- Phase bindings: default_bindings() or cross-driver (empty)?

For LARGE scope: delegate to `av-architect` agent for full validation.

## Step 4 — TDD (RED)

Delegate to `av-tester` agent:
- Write failing tests in `tests/test_{name}.cu`
- Required: component basics, lifecycle, integration, domain validation
- Build and confirm tests fail for the right reason

## Step 5 — Implement (GREEN)

1. Header in `include/atomiverse/{module}/` (CUDA-free, xtensor-free)
2. Implementation: CUDA in `src/cuda/{module}/*.cu`, CPU in `src/cpu/{module}/*.cpp`
3. Register in CMakeLists.txt
4. No template parameters for Real/Int — these are compile-time typedefs
5. Build and confirm tests pass

## Step 6 — Verify

Run in parallel where possible:
- Build: `cmake --build build -j$(nproc)`
- Tests: `ctest --test-dir build --output-on-failure`
- For MEDIUM/LARGE: delegate to `av-architect` agent for architecture validation
- For MEDIUM/LARGE: delegate to `av-documenter` agent for documentation

Summarize: files changed, tests passing, remaining TODOs.

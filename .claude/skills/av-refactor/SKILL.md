---
name: av-refactor
description: Restructure code while maintaining all architectural invariants.
argument-hint: "<what to refactor and why>"
user-invocable: true
---

# Refactor

Read CLAUDE.md for architecture rules and conventions.

## Procedure

1. **Snapshot invariants.** Before changing anything:
   ```bash
   cmake --build build -j$(nproc)
   ctest --test-dir build --output-on-failure
   ```
   Record which tests pass and the current API surface.

2. **Understand scope.** What is being restructured and why? Identify:
   - What moves where / gets renamed
   - What interfaces change
   - What must stay the same (invariants to preserve)

3. **Check architecture.** Delegate to `av-architect` agent on the current state to ensure the target area is well-understood.

4. **Execute incrementally.** One logical change at a time. Build + test after each change.

5. **Verify invariants.** After completion:
   - Same tests pass (no regressions)
   - Delegate to `av-architect` agent to confirm architecture compliance
   - No new backend-specific headers in facade (no CUDA, no xtensor), no lifecycle violations

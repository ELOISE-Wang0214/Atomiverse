---
name: av-test
description: Run tests and analyze test quality for Atomiverse.
argument-hint: "[module or test path]"
user-invocable: true
---

# Test & Coverage Analysis

Read CLAUDE.md for build/test commands.

1. Build and run tests:
   ```bash
   cmake --build build -j$(nproc)
   ctest --test-dir build --output-on-failure
   ```
   If a module is specified, filter with `-R "{module}"`.

2. Delegate to `av-tester` agent for test quality and coverage analysis.

3. Report: pass/fail/skip counts, per-file breakdown, missing test categories.

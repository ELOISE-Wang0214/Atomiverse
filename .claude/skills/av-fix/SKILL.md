---
name: av-fix
description: Quick fix workflow -- diagnose, fix, verify in minimal steps.
argument-hint: "<bug description, error message, or failing test>"
user-invocable: true
---

# Quick Fix

Read CLAUDE.md for build commands and conventions.

## Procedure

1. **Reproduce.** Build and run to confirm the failure:
   ```bash
   cmake --build build -j$(nproc) 2>&1 | tail -50
   ctest --test-dir build --output-on-failure
   ```

2. **Diagnose.** Follow the diagnosis procedure from `/av-debug` (classify → gather → diagnose by type). Stop once root cause is identified.

3. **Fix.** Make the minimal change that resolves the issue.
   - Check CLAUDE.md if the fix touches architecture boundaries (lifecycle, backend isolation).

4. **Verify.** Build and run full test suite — confirm no regressions.

5. **Report.** One-line summary of root cause and fix.

---
name: av-arch
description: Validate Atomiverse architecture compliance.
argument-hint: "[path or module]"
user-invocable: true
---

# Architecture Validation

Delegate to the `av-architect` agent with the user's argument (or "all files" if none given).

The agent will:
1. Read CLAUDE.md for current architecture rules
2. Discover current files via glob
3. Run all architecture checks (backend isolation, lifecycle, Frame/Driver separation, naming)
4. Report findings in severity-sorted `[SEVERITY] file:line` format

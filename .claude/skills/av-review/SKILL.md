---
name: av-review
description: Comprehensive code review across all dimensions.
argument-hint: "[path or list of modified files]"
user-invocable: true
---

# Code Review

Read CLAUDE.md for architecture rules and conventions.

Determine the files to review: use the argument, or `git diff --name-only` for recent changes.

## Review Dimensions

Run these in parallel by delegating to agents:

1. **Architecture** — delegate to `av-architect` agent
2. **Performance** — delegate to `av-optimizer` agent
3. **Documentation** — delegate to `av-documenter` agent
4. **Domain correctness** — delegate to `av-scientist` agent

## Output

Aggregate findings into a summary table:

| Dimension          | Critical | High | Medium | Low |
|--------------------|----------|------|--------|-----|
| Architecture       |          |      |        |     |
| Performance        |          |      |        |     |
| Documentation      |          |      |        |     |
| Domain Correctness |          |      |        |     |

Verdict: APPROVE / REQUEST CHANGES / BLOCK

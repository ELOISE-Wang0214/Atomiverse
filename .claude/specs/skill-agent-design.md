---
title: Skill / Agent Design Guideline (Atomiverse Reference Model)
status: approved
created: 2026-04-17
---

# Skill / Agent Design Guideline

This spec is **descriptive** of the current Atomiverse `.claude/{skills,agents}/`
layout and **prescriptive** as a checklist for re-examining other projects. It
captures *what each skill/agent owns* and *the design rules that produced the
orthogonal split*. When auditing another repository, walk Sections 4 and 5
against its `.claude/` tree and produce a delta report.

## 1. Two-Layer Model

```
user ──▶ /skill  (workflow + orchestration, user-invocable)
              │
              ▼
          agent  (single expertise axis, never user-invocable directly)
              │
              ▼
       tools (Read/Grep/Glob/Bash, sometimes Write/Edit)
```

- **Skill** = a *verb* the user types (`/av-impl`, `/av-fix`). It owns the
  *procedure*: scope assessment, ordering, parallel fan-out, when to delegate,
  how to aggregate results. A skill is short — it reads CLAUDE.md, classifies
  the request, calls one or more agents, summarizes.
- **Agent** = a *role* with one expertise. It owns *checks, knowledge, and
  edits inside its domain*. It is stateless across invocations and is reached
  only through a skill (or transitively from another agent only when the
  agent's frontmatter explicitly allows it — see Rule O3).

The two layers must not collapse. A skill must not contain expert knowledge
that belongs in an agent. An agent must not contain workflow choreography that
belongs in a skill.

## 2. Why This Split

Three pressures justify the cost of maintaining two layers:

1. **Reusability.** Multiple skills delegate to the same agent (`av-architect`
   is used by `/av-impl`, `/av-refactor`, `/av-arch`, `/av-review`). Putting
   the architecture rules in one place keeps them coherent.
2. **Parallelism.** A skill like `/av-review` fans out to four agents in one
   message; agents have isolated context windows so the orchestrator stays
   small.
3. **Safety.** Read-only agents (`av-architect`, `av-cuda`, `av-optimizer`,
   `av-scientist`) cannot accidentally edit code; write-capable agents
   (`av-tester`, `av-hpc`, `av-documenter`) have edit rights only inside
   their declared scope.

## 3. Knowledge Hierarchy

```
.claude/NOTES.md   evolving decisions, captured by /av-note
        │                       │
        ▼ promotion when stable │
   CLAUDE.md  ◀──────── single source of truth — read by every agent
        │
        ▼
   skills/agents read from it; never duplicate it
```

Every agent's first line is *"Read CLAUDE.md and .claude/NOTES.md before
running any checks."* Project-specific facts live in CLAUDE.md or NOTES.md, not
in agent prompts. Agent prompts hold only **unique knowledge** the agent needs
that does not belong in CLAUDE.md (algorithm patterns, grep heuristics,
tolerances, etc.).

## 4. Per-Skill Responsibilities

A skill manifest (`SKILL.md`) is at most ~80 lines. It states scope, procedure,
and what to delegate. It does **not** restate CLAUDE.md.

### Plan phase

| Skill | Owns | Delegates to | Writes |
|---|---|---|---|
| `/av-spec` | Convert request → structured spec at `.claude/specs/{slug}.md`; update INDEX | `av-scientist` (if physics) | `specs/` |
| `/av-litrev` | Establish credibility (papers, reference implementations, validation cases) before code is written | `av-scientist` (+ web) | none |

### Build phase

| Skill | Owns | Delegates to | Writes |
|---|---|---|---|
| `/av-impl` | Scope classification (SMALL/MEDIUM/LARGE), pipeline orchestration | `av-architect`, `av-tester`, `av-documenter` | code, tests, docs |
| `/av-fix` | Minimal-diff bug fix loop (reproduce → diagnose → fix → verify) | reuses `/av-debug` procedure inline | code |
| `/av-refactor` | Snapshot invariants → incremental restructure → re-verify | `av-architect` | code |

### Diagnose phase

| Skill | Owns | Delegates to | Writes |
|---|---|---|---|
| `/av-debug` | Classify failure (build / test / runtime), gather, root-cause | none | **never** |

### Verify phase

| Skill | Owns | Delegates to | Writes |
|---|---|---|---|
| `/av-review` | Aggregate multi-axis findings into severity table + verdict | `av-architect` ‖ `av-optimizer` ‖ `av-documenter` ‖ `av-scientist` (parallel) | none |
| `/av-arch` | Layering / lifecycle / backend isolation | `av-architect` | none |
| `/av-cuda` | GPU execution model classification + sync audit | `av-cuda` | none |
| `/av-perf` | Cross-backend perf anti-patterns | `av-optimizer` | none |
| `/av-hpc` | CPU parallelism review **or** implementation | `av-hpc` | code (impl mode) |
| `/av-test` | Run ctest, surface gaps | `av-tester` | tests |

### Document & memory

| Skill | Owns | Delegates to | Writes |
|---|---|---|---|
| `/av-docs` | Narrative tutorial / Doxygen audit | `av-documenter` | docs |
| `/av-note` | Capture decision + conflict / stale sweep + promote stable rules to CLAUDE.md | none | NOTES.md, CLAUDE.md |

## 5. Per-Agent Responsibilities

An agent manifest declares `tools:` (capability surface) and `model: inherit`.
It states **role**, **unique knowledge** (not in CLAUDE.md), **procedure**, and
**output format** (`[SEVERITY] file:line` for review agents).

| Agent | Single axis of expertise | Tools | Read-only? |
|---|---|---|---|
| `av-architect` | Module/layer boundaries, Driver/Component/Frame separation, lifecycle | R/G/Glob/Bash | yes |
| `av-cuda` | GPU execution model — compute / CPU-bridge / host-dependent classification, sync discipline | R/G/Glob/Bash | yes |
| `av-hpc` | CPU parallelism — OpenMP, MPI, SIMD/xsimd, cache, ERI hot path | + Write/Edit | **no** |
| `av-optimizer` | Cross-backend perf anti-patterns and hot-path measurement | R/G/Glob/Bash | yes |
| `av-scientist` | Physics correctness — equations, units, conservation laws, references | + WebSearch/WebFetch | yes |
| `av-tester` | TDD cycle — RED→GREEN, GoogleTest, lifecycle/domain test categories, tolerances | + Write/Edit | **no** |
| `av-documenter` | Doxygen audit + narrative tutorials in the project's textbook voice | + Write/Edit | **no** |

Each agent solves one and only one of: *"is this well-structured?"*,
*"is this fast?"*, *"is this correct in the GPU model?"*, *"is this fast on
CPU?"*, *"is this physically correct?"*, *"is this tested?"*, *"is this
documented?"*.

## 6. Design Rules (the orthogonality contract)

These are the rules to enforce when designing or auditing skill/agent sets.
They are stated as **must / should / never**.

### Orthogonality (O)

- **O1.** Each agent owns exactly one expertise axis. If two agents would both
  answer the same question, merge or split until each has a unique axis.
- **O2.** Agents do not call other agents. The skill is the only orchestrator.
  *Exception:* `/av-hpc`'s implementation mode may write code that the
  `av-architect` agent will later validate — but the validation call is made
  by the parent skill, not by `av-hpc` itself.
- **O3.** A skill names its delegates explicitly. There is no implicit routing.
- **O4.** A skill does not contain expert knowledge. If a check is non-trivial
  enough to need worked examples or reference data, it belongs in an agent.

### Knowledge locality (K)

- **K1.** CLAUDE.md is the single source of truth for stable rules. Every
  agent reads it first.
- **K2.** NOTES.md holds evolving decisions; `/av-note` enforces conflict
  detection and promotes stable entries into CLAUDE.md (then deletes from
  NOTES).
- **K3.** Agent prompts contain only *unique knowledge* not in CLAUDE.md
  (heuristics, grep patterns, tolerances). They never duplicate architecture
  rules.
- **K4.** Skill prompts contain only *workflow* (steps, when to delegate,
  output format). They never restate domain rules.

### Capability hygiene (C)

- **C1.** Read-only agents (`av-architect`, `av-cuda`, `av-optimizer`,
  `av-scientist`) declare only Read/Grep/Glob/Bash. Granting Write/Edit must
  be justified by the agent's role.
- **C2.** A skill that writes (specs, code, docs, tests) must say so in its
  description so the user knows what will change.
- **C3.** Diagnose-only skills (`/av-debug`) state explicitly that they never
  edit. This must be enforced in the skill's procedure, not just in the
  description.

### Workflow shape (W)

- **W1.** Plan → Build → Verify → Document → Memory is the canonical loop.
  Every project should have at least one skill in each phase.
- **W2.** RED before GREEN: the implementation skill must call the test agent
  to write failing tests before the implementation agent runs.
- **W3.** Reviews fan out in parallel. The review skill issues all delegate
  calls in one message.
- **W4.** Architecture validation happens at three points: scope assessment
  (planning), after implementation, after refactor.
- **W5.** Bug fixing is a separate, minimal loop (`/av-debug` → `/av-fix`),
  not a degenerate case of `/av-impl`. Keep them distinct.

### Output format (F)

- **F1.** Review agents output severity-sorted `[SEVERITY] file:line —
  message` lines. Skills aggregating multiple agents render a Critical /
  High / Medium / Low table and a verdict (APPROVE / REQUEST CHANGES /
  BLOCK).
- **F2.** Skills end with a one-line summary suitable for the user (files
  changed, tests passing, remaining TODOs).

## 7. Re-Examination Checklist (apply to any project)

Use this checklist when auditing another project's `.claude/` directory.
Output one finding per row, with `[SEVERITY] file:line — message`.

### Layer presence

- [ ] Are skills (`.claude/skills/*/SKILL.md`) and agents (`.claude/agents/*.md`)
      both present, or is the project conflating the two?
- [ ] Does every user-visible verb live in a skill, not an agent?
- [ ] Is there a CLAUDE.md serving as single source of truth?

### Orthogonality

- [ ] List each agent and the question it answers. Are any two answering the
      same question? (Rule O1)
- [ ] Do agents invoke other agents? (Rule O2 violation)
- [ ] Do skills name delegates explicitly, or do they rely on implicit
      routing? (Rule O3)
- [ ] Does any skill embed expert knowledge that belongs in an agent? (Rule O4)

### Knowledge

- [ ] Do agents read CLAUDE.md first? (Rule K1)
- [ ] Is there a `/note`-like skill enforcing conflict detection between
      NOTES.md and CLAUDE.md, plus promotion of stable rules? (Rule K2)
- [ ] Do agent prompts duplicate CLAUDE.md content? (Rule K3 violation)
- [ ] Do skill prompts restate domain rules instead of pointing to CLAUDE.md?
      (Rule K4 violation)

### Capability

- [ ] Are read-only agents declared with read-only tools? Any unjustified
      Write/Edit? (Rule C1)
- [ ] Does each writing skill announce what it writes? (Rule C2)
- [ ] Are diagnose-only skills enforced as such in their procedure, not just
      labeled? (Rule C3)

### Workflow

- [ ] Plan / Build / Verify / Document / Memory — does the project have at
      least one skill in each phase? (Rule W1)
- [ ] Does the implementation skill enforce RED before GREEN? (Rule W2)
- [ ] Do review skills fan out in parallel? (Rule W3)
- [ ] Are architecture checks scheduled at the three required points?
      (Rule W4)
- [ ] Are debug and fix kept as separate minimal loops? (Rule W5)

### Output

- [ ] Do review agents emit `[SEVERITY] file:line` lines? (Rule F1)
- [ ] Do skills end with a one-line user-facing summary? (Rule F2)

### Per-project domain coverage

For a domain-heavy project, verify there is one agent per *axis of risk* the
project actually has. For Atomiverse those axes are: architecture, GPU
execution, CPU HPC, cross-backend perf, science, tests, docs. For another
project they will differ — but every axis where the project can fail
silently must have its own agent, and every agent must map to a real risk.

## 8. Out of Scope

- This spec does not prescribe the **set** of skills/agents another project
  must have. It prescribes the *shape* of that set: orthogonality, knowledge
  locality, capability hygiene, workflow phases, output format. The actual
  list (`av-cuda`, `av-hpc`, etc.) is Atomiverse-specific.
- This spec does not cover slash-command implementation details (the harness
  side of skills/hooks). It addresses only the design contract.

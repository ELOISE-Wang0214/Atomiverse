---
name: av-spec
description: Convert feature request into a structured technical specification, persisted as a trackable project artifact.
argument-hint: "<feature description>"
user-invocable: true
---

# Specification Generator

Read CLAUDE.md for architecture, component model, and conventions.

## Procedure

1. **Parse the request.** What capability is being added? Which layer is affected?

2. **Search domain basis.** If physics is involved, delegate to `av-scientist` agent.

3. **Scan codebase for interaction points:**
   - Glob `include/atomiverse/**/*.h` to discover existing components as patterns
   - Check if new Frame fields are needed
   - Check if a new Driver type is needed

4. **Output structured spec** with these sections:
   - Summary (1 paragraph)
   - Domain basis (equations, references, units) — if applicable
   - Component design (class, constructor params, bindings, lifecycle)
   - Frame changes (if any)
   - Dependency graph (what blocks what)
   - Files to create/modify
   - Testing strategy (lifecycle, correctness, integration, domain)

5. **Persist the spec.** Write the spec to `.claude/specs/{slug}.md` where `{slug}` is a short kebab-case name (e.g., `amorphous-pe`, `nose-hoover`). Use this frontmatter:

   ```markdown
   ---
   title: <spec title>
   status: draft | approved | in-progress | done
   created: YYYY-MM-DD
   ---
   ```

6. **Update the spec index.** Append a one-line entry to `.claude/specs/INDEX.md`:
   ```
   - [{slug}]({slug}.md) — <one-line summary> [status]
   ```
   Create INDEX.md if it doesn't exist.

7. **Report to user.** Show the spec file path so they can review and approve. Update `status: approved` when the user confirms.

## Lifecycle

Specs are living documents. When implementation begins (`/av-impl`), update status to `in-progress`. When all items are complete and verified, update to `done`. Use `/av-spec` again on the same slug to revise an existing spec instead of creating a new one.

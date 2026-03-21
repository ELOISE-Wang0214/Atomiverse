---
name: av-note
description: Capture an architectural decision or correction into .claude/NOTES.md, with conflict detection and stale-note cleanup.
argument-hint: "<what to remember>"
user-invocable: true
---

# Architecture Note

Capture an important architectural decision, correction, or lesson learned.

## Procedure

### 1. Conflict check (MANDATORY — run before any write)

Read both `.claude/NOTES.md` and `CLAUDE.md`. For the new note, check:

- **Duplicate**: An existing note or CLAUDE.md section already says the same thing → tell the user, don't add.
- **Contradiction**: An existing note or CLAUDE.md rule says the **opposite** → flag the conflict explicitly:
  ```
  CONFLICT: new note "<X>" contradicts existing "<Y>" in {file}:{line}.
  ```
  Ask the user which is correct. Then:
  - Update or delete the wrong entry.
  - If the wrong entry is in CLAUDE.md, fix CLAUDE.md.
  - If the wrong entry is in NOTES.md, delete or rewrite it.
- **Supersede**: The new note refines or replaces an older note → update the old note in-place (don't keep both).

### 2. Scan for stale notes

Every time this skill runs, also scan ALL existing NOTES.md entries:

- **Outdated**: Does the note reference files, functions, or patterns that no longer exist? Grep the codebase to verify. Delete notes that reference dead code.
- **Already promoted**: Is the note's rule now fully covered in CLAUDE.md? If yes, delete it from NOTES.md.
- **Contradicts current code**: Does the note say X but the code does Y? Flag to the user for resolution.

Report any cleanup actions taken.

### 3. Decide where the new note belongs

- Stable, proven rule already verified in code → add to CLAUDE.md directly.
- New decision that may still evolve → add to .claude/NOTES.md.
- Already covered → tell the user, don't duplicate.

### 4. Write the entry

In NOTES.md, use this format:

```markdown
---

## [YYYY-MM-DD] Short title

Description of the decision or correction.

**Rule**: The concrete rule to follow.
```

### 5. Promotion check

Check if any existing NOTES.md entries are now stable enough to promote to CLAUDE.md:
- The rule is verified in current code (grep to confirm).
- The rule has not changed in recent conversations.
- If promoting: add to the appropriate CLAUDE.md section, then delete from NOTES.md.

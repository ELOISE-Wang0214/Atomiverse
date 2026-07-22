---
spec: hf-deferred-backlog
created: 2026-07-21
criteria:
  - id: ac-001
    summary: Inventory is the single residual HF parking lot after campaign close
    type: docs
    pass_when: >
      .claude/specs/hf-deferred-backlog.md exists with P0–P5 checklists and
      permanent non-goals; no hf-capability-gaps-* or
      scf-driver-refactor-06-* files remain under .claude/specs/; INDEX.md
      points at this inventory as the open HF residual entry.
    status: verified
    last_checked: 2026-07-21
    verified_by: human
  - id: ac-002
    summary: Umbrella is not an impl target
    type: docs
    pass_when: >
      Spec frontmatter and Summary state that /mol:impl must not run on
      hf-deferred-backlog; work is carved into child specs; Tasks do not
      describe implementing P0–P5 under this slug.
    status: verified
    last_checked: 2026-07-21
    verified_by: human
  - id: ac-003
    summary: Closed campaign is recorded as superseded, not reopened
    type: docs
    pass_when: >
      Spec notes supersedes hf-capability-gaps-01…21 and
      scf-driver-refactor-06-criterion-unify; residual unfinished hard gates
      (CFMM, x2c, C-PCM, embed) appear under P1 rather than as resurrected
      old slugs.
    status: verified
    last_checked: 2026-07-21
    verified_by: human
out_of_scope:
  - Implementing any P0–P5 checklist item under this slug
  - Restoring deleted hf-capability-gaps-* acceptance ledgers
  - PySCF feature parity as a global success metric
---

# Acceptance — hf-deferred-backlog

This slug **records** residual HF work; it does not deliver physics. "Done"
for the parking lot itself means: campaign specs closed, inventory written,
and process rules (child-spec only) are explicit. Individual checklist rows
are closed only by **child** specs later.

## AC-001 — Inventory replaces the open HF campaign

After the 2026-07-21 manual close, `.claude/specs/` must not retain the
old capability-gap chain files. The deferred inventory holds P0–P5 and
permanent non-goals.

## AC-002 — No mega-impl

Operators and automation must not run `/mol:impl hf-deferred-backlog`.
Work starts with a new `/mol:spec` on one theme.

## AC-003 — History lives in git, residuals live in P1

Unfinished oracle/docs/hard-gates from the closed chain are listed under
P1 so they can be re-specified tightly without resurrecting 21 files.

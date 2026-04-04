---
name: av-docs
description: Write or rewrite Atomiverse documentation in a narrative, textbook-style voice.
argument-hint: "<feature name or spec slug>"
user-invocable: true
---

# Tutorial Writer

Delegate to the `av-documenter` agent with full context, but do not hand it a
rigid outline and do not optimize for length.

Before delegating, collect and verify:
1. The relevant spec or design note, if one exists.
2. The implementation headers and source files that define the current behavior.
3. The existing user-facing documentation for the same topic.
4. The real example program, build target, and output artifacts that the text
   will mention.

If a path, target, plot, or API does not exist, do not mention it as if it
does. Either correct the tutorial to the current codebase or add the missing
artifact as part of the task.

## What a good tutorial looks like

Write for a technically literate reader who is new to the domain. The result
should read like a short chapter, not like a generated checklist.

Prefer a continuous narrative:
- Start with the physical or computational problem in plain language.
- Introduce notation only after the reader knows why it matters.
- Let equations support the explanation instead of replacing it.
- Braid theory and code together. The matching example should appear inside the
  teaching flow, not as an appendix.
- Move from intuition to code and back again in the order that best fits the
  topic. Do not force a fixed section template onto every tutorial.
- Keep each tutorial anchored to one concrete numbered example under
  `examples/tutorial/`. A CUDA companion may share the same tutorial number as
  a secondary translation, but the main reading path should still stand on its
  own.
- Put usage, setup, backend integration, and framework interop material in the
  user guide rather than in the computational-chemistry tutorial track.

## Style rules

- Prefer prose over lists. Use bullets only when sequence or enumeration is
  genuinely clearer than paragraphs.
- Do not add a table of contents unless the user explicitly asks for one.
- Avoid large comparison tables. If a fact can be stated in a sentence, write a
  sentence.
- Avoid numbered heading scaffolds like "1. What problem are we solving?" unless
  the material is truly procedural.
- Keep sections compact. The tutorial should be as long as necessary and no
  longer. There is no target line count.
- Use concrete file names, real commands, and real output from the repository.
- Do not fabricate end-to-end examples, plots, or build targets.
- Write mathematics in standard Markdown LaTeX: `$...$` for inline math and
  `$$...$$` for display math. Never use Doxygen `\f` markers in Markdown.
- Every symbol must be defined in prose.
- Do not foreground library names, container types, backend plumbing, header
  paths, or internal execution phases unless they are the subject of the page.

## Technical bar

- Verify that every referenced file and API name exists in the current tree.
- When describing the algorithm, follow the current implementation, not an
  older version of the design.
- If the code has evolved, explain the current structure plainly instead of
  preserving stale historical wording.
- Every tutorial should have a nearby numbered example program that a reader can
  build and run while reading. Prefer CPU examples unless the topic is
  inherently backend-specific and belongs in the user guide.

---
name: av-docs
description: Write or rewrite a textbook-style feature tutorial with E2E example and plots.
argument-hint: "<feature name or spec slug>"
user-invocable: true
---

# Tutorial Writer

Delegate to the `av-documenter` agent with full context.

Before delegating, collect:
1. The spec file at `.claude/specs/{slug}.md` (if it exists)
2. The implementation headers and source files
3. Any existing tutorial at `docs/tutorial/{slug}.md`
4. The closest existing tutorial as a style reference (e.g., `docs/tutorial/h2-rhf.md`)

Pass all of this to `av-documenter` with the instruction to write a **tutorial** (not a Doxygen audit).

## What a correct tutorial looks like

### Audience
Assume **zero domain background**. The reader may be a software engineer, physicist,
or student with no computational chemistry experience. Never assume familiarity with
quantum mechanics, Hartree-Fock, or electronic structure.

### Required structure (in this order)

1. **Motivation** — What physical problem are we solving? Why does it matter?
   Use an analogy that requires no domain knowledge. One page max.

2. **Theory from first principles** — Build up to the method step by step.
   Each equation must be preceded by a paragraph explaining what it means physically.
   No equation without prose. No symbol without definition.

3. **The algorithm** — Walk through the SCF loop or solver steps using
   pseudocode-style English before showing any C++. Annotate each step with
   its physical meaning.

4. **Atomiverse API walkthrough** — Show exactly how to call the solver.
   Every field in the input/output struct must be explained.
   Include the complete include list so the reader can copy-paste to compile.

5. **Key path analysis** — Trace through the most important code paths
   (`build_uhf_G`, `ortho_diag`, etc.) with line-level annotations.
   Explain *why* each step is done, not just *what* it does.

6. **E2E example** — A complete, compilable Atomiverse example in `examples/`.
   Must produce **at least two plots** using matplotplusplus (guarded by
   `#ifdef ATV_HAS_MATPLOT`) that verify the physics:
   - One plot showing the key physical observable (e.g., potential energy curve)
   - One plot showing a diagnostic (e.g., convergence, spin contamination)
   Reference the build command (`cmake -B build -DATV_USE_MATPLOT=ON ...`).

7. **Interpreting the output** — Walk through the expected plot shapes and
   numerical values. Explain what "correct" looks like and what failure looks like.

8. **Common pitfalls** — What goes wrong if the user sets N_α wrong?
   What does non-convergence look like?

### Style rules
- Prose paragraphs, not bullet lists, for theory sections
- **All math must use LaTeX notation.** Inline: `$...$`. Display: `$$...$$` on its own line.
  Never use Unicode math symbols (α β γ Σ ∑ √) or ASCII pseudo-equations in prose or
  code blocks. Every equation must be typeset.
- Concrete numbers over abstract statements ("$-0.4666$ Hartree" not "negative")
- No jargon without definition ("Hartree-Fock is a method for..." not just "HF")
- No forward references ("as we saw above" requires the thing was actually above)

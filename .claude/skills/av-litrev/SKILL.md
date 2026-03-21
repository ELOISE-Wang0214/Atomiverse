---
name: av-litrev
description: Literature and reference-implementation review for domain-critical features.
argument-hint: "<method or algorithm name>"
user-invocable: true
---

# Literature Review

Delegate to `av-scientist` agent for domain expertise, then extend with:

## Procedure

1. **Literature search.** Use WebSearch for foundational papers, reviews, corrections.
2. **Document equations.** Governing equations with all terms defined, in project units (eV, Angstrom, fs, amu).
3. **Reference implementations.** Search LAMMPS, VASP, ASE, etc. for API patterns and validation cases.
4. **Validation targets.** Analytical solutions, published benchmarks, conservation laws, limiting cases.

## Output

- Credibility assessment: YES / NO / CONDITIONAL
- Key references (numbered list)
- Implementation guidance (algorithm variant, numerical details, pitfalls)
- Validation plan (ordered tests with expected values)
- Open questions

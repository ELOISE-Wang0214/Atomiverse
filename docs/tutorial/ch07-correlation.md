# Chapter 7: Beyond Mean Field — Correlation Methods

Every electronic-structure method in the previous chapters — Hartree-Fock
([Tutorial 2](02-rhf.md), [Tutorial 3](03-uhf.md)) and Kohn-Sham DFT
([Tutorial 4](04-uks.md), [Chapter 6](ch06-hybrid-dft.md)) — is a **mean-field**
theory. Each electron moves in the *average* field of all the others. That
averaging is what makes the SCF loop tractable, and it is also its central
approximation: real electrons avoid each other instantaneously, and mean-field
theory does not fully capture that avoidance.

This is a **status chapter**. It explains what mean-field theories miss, sketches
the hierarchy of methods built to recover it, and states plainly where
Atomiverse stands: these wavefunction correlation methods are **not implemented**
— they are out of scope for the current engine, documented here so the manual's
map of computational chemistry is complete and honest, not so you can run one
today.

You will learn:

1. **The concept** of electron correlation and the correlation energy.
2. **The hierarchy** — MP2, configuration interaction, coupled cluster — and the
   cost-versus-accuracy trade that orders them.
3. **What Atomiverse provides today**: DFT as the practical, affordable route to
   correlation, and why that is the engine's answer for now.
4. **Where a correlation method would live** in the Driver/Component
   architecture if it were added — which turns out to reinforce the same
   pattern the rest of the book teaches.

## The Concept: What Mean-Field Theory Misses

For a given basis set, Hartree-Fock gives the best possible single Slater
determinant — the best wavefunction of the form "each electron in its own
orbital." The exact non-relativistic energy is lower. The difference is, by
definition, the **correlation energy**:

$$ E_\text{corr} = E_\text{exact} - E_\text{HF}. $$

It is a small fraction of the total energy — typically under 1% — but it is
comparable to, or larger than, the **chemical** energies we care about: reaction
enthalpies, barrier heights, binding energies. Getting total energies to 99%
and still being useless for chemistry is the recurring hazard of this field, and
correlation is why.

Two flavours are worth naming:

- **Dynamic correlation** — the short-range, instantaneous avoidance of
  electrons (the "Coulomb hole"). Captured by mixing in many determinants that
  each promote electrons into virtual orbitals.
- **Static (non-dynamic) correlation** — the near-degeneracy that appears when a
  single determinant is qualitatively wrong (bond breaking, diradicals,
  transition metals). This one breaks single-reference methods and needs
  multi-reference treatments.

Hartree-Fock captures **neither**. Kohn-Sham DFT ([Chapter 6](ch06-hybrid-dft.md))
folds dynamic correlation into the exchange-correlation functional — an
*approximate, implicit* treatment that is astonishingly cost-effective but not
systematically improvable. The wavefunction hierarchy below is the *explicit,
systematically improvable* alternative, bought at steeply rising cost.

## The Hierarchy: MP2, CI, Coupled Cluster

Post-Hartree-Fock methods start from the converged HF orbitals and add
correlation by allowing electrons to occupy the empty (virtual) orbitals. They
form a rough ladder of accuracy and cost, all scaling far worse than the
$\mathcal{O}(N^3)$–$\mathcal{O}(N^4)$ of SCF:

| method | idea | formal cost | role |
|---|---|---|---|
| **MP2** | 2nd-order perturbation theory in the fluctuation potential | $\mathcal{O}(N^5)$ | cheapest correlation; good for dispersion, weak on metals/π |
| **CISD / CCSD** | include single + double excitations (variational / exponential) | $\mathcal{O}(N^6)$ | CCSD is a reliable workhorse for closed-shell molecules |
| **CCSD(T)** | CCSD plus a perturbative triples correction | $\mathcal{O}(N^7)$ | the **"gold standard"** of single-reference chemistry |
| **Multi-reference (CASSCF, MRCI)** | several reference determinants | worse still | for static correlation: bond breaking, diradicals |

The pattern is a treadmill: each rung buys roughly an order of magnitude in the
error at the price of a higher power of the system size. CCSD(T) at the complete
basis-set limit is accurate enough to compete with experiment for small
molecules, and completely intractable for a protein. **Coupled cluster is not
variational** (its energy can dip below the exact energy), which is the price of
its favourable size-consistency. These are the trade-offs a practitioner weighs;
this book's job here is to place them on the map, not to teach how to run them.

## What Atomiverse Provides Today

Atomiverse's answer to correlation is **density-functional theory**. A hybrid
functional such as B3LYP or PBE0 ([Chapter 6](ch06-hybrid-dft.md)) delivers much
of the accuracy of an MP2/CCSD calculation for a fraction of the cost, by
approximating correlation inside the exchange-correlation functional rather than
constructing it explicitly from excited determinants. For the great majority of
ground-state structures, energies, and forces — and for every gradient-driven
workflow in this manual (AIMD, geometry optimization, vibrations) — DFT is the
method of choice, and it is what the engine is built around.

What Atomiverse does **not** provide:

- No MP2, no configuration interaction, no coupled cluster, no multi-reference
  methods. There is no `mp2_solve`, no CC amplitude solver, no post-SCF
  correlation module anywhere in the codebase.
- Consequently, no systematically-improvable path to the exact energy, and no
  first-principles treatment of static correlation.

If your science *requires* CCSD(T)-quality correlation or a genuine
multi-reference treatment, Atomiverse is not the right tool today, and this
chapter is the honest place to say so.

## Where Correlation Would Live in the Architecture

It is worth sketching where a correlation method *would* fit, because the answer
reinforces the architectural through-line of the whole manual.

A post-Hartree-Fock method is not a new kind of loop. It is a **post-SCF
step**: it consumes the *converged* output of an SCF solve — the occupied and
virtual molecular orbitals, their orbital energies, and the two-electron
integrals transformed into the MO basis — and returns a correlation energy (and,
in principle, a correlated density and gradient).

In Atomiverse's vocabulary that is not a new `Driver` category. The Gaussian SCF
engine already exposes exactly the inputs such a method needs:

- `cpu::scf::build_integrals(frame, basis, charge, /*compute_eri=*/true)`
  produces the two-electron integrals (`integrals.ERI`) that an MO-basis
  transform would consume — the same integrals the CUDA `RHFFock` kernel reads.
- `cpu::scf::rhf_solve(...)` returns the converged orbitals and their energies.

A correlation method would therefore be a component that **reads the SCF result
and writes a corrected energy onto `frame.scalars`** — architecturally an
*evaluator that runs after the SCF evaluator*, in the same slot pattern that
[Chapter 13](ch13-observables.md)'s observers occupy after the integrator. It
would reuse the existing integral engine rather than introduce a parallel one,
exactly as the hybrid-gradient routine in [Chapter 6](ch06-hybrid-dft.md) reuses
the ERI machinery. "One loop, many physics" extends cleanly to "one SCF result,
many post-processings" — the door is architecturally open even though nothing
walks through it yet.

## Status and Further Reading

**Status:** planned / not implemented. This chapter documents intent and
context only; there is no runnable example and no validation case, because there
is no code to validate. When the manual's method-status matrix
([`.claude/specs/manual-discovery.md`](../../.claude/specs/manual-discovery.md) §5)
lists a correlation method as `absent (planned)`, this chapter is what that row
points to.

For the methods Atomiverse *does* provide as its route to correlation, return to
[Chapter 6: Hybrid Functionals and Analytic Gradients](ch06-hybrid-dft.md). For
the periodic counterpart of the electronic-structure problem, continue to
[Tutorial 5: Plane-Wave DFT for Periodic Systems](05-pw.md).

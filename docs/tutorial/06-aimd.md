# Tutorial 6: Ab Initio Molecular Dynamics (AIMD)

The previous five chapters showed how to compute, in turn, a one-shot SCF
energy, a closed-shell RHF wave function, an open-shell UHF wave function, a
Kohn-Sham DFT calculation, and an L-BFGS geometry optimization. Each chapter
asked the same question — *what does the molecule look like at one
geometry?* — and answered it with a different level of theory.

This chapter changes the question. We no longer want one geometry; we want
a *trajectory*. Atoms have masses; given a force on each atom we can
integrate Newton's equations and watch the molecule evolve in time. When
the forces come from a self-consistent electronic-structure calculation
performed at every time step, the procedure is called **ab initio
molecular dynamics**, and the most common variant — running a fully
converged SCF at every nuclear configuration — is **Born-Oppenheimer MD
(BOMD)**.

The example to read alongside this chapter is
[`examples/06_aimd.cpp`](../../examples/06_aimd.cpp).
It propagates $\mathrm{H}_2$ for 50 fs in the canonical (NVT) ensemble at
100 K, using RHF/STO-3G forces and a Nose-Hoover thermostat.

## The Driver / Component Recipe

Atomiverse does not ship an `AIMDDriver` class. AIMD is built from
components you have already met:

- The **MD driver** `cpu::md::MD` owns the time-stepping loop. It is the
  same driver that runs classical molecular dynamics with Lennard-Jones
  forces.
- A **velocity-Verlet integrator** `cpu::md::VelocityVerlet` advances
  positions and velocities each step.
- An **SCF evaluator** `cpu::scf::SCFEvaluator` plays the role of the force
  field. It is bound at `md::Phase::force` — the same slot a classical
  pair potential would occupy. At each MD step it rebuilds molecular
  integrals at the current geometry, runs SCF, computes the analytic
  nuclear gradient, and writes forces into `frame.atoms.fx/fy/fz`.

The three pieces compose:

```cpp
cpu::scf::SCFEvaluator scf{config};
scf.predictor(/*order=*/3);                       // ASPC density warm-start

cpu::md::MD sim;
sim.add(cpu::md::InitVelocity{T_target, seed});   // Maxwell-Boltzmann init
sim.add(cpu::md::VelocityVerlet{dt});             // 0.5 fs per step
sim.add(cpu::md::NoseHoover{T_target, tau, dt});  // NVT
sim.on(md::Phase::force).add(std::move(scf));     // SCF as the "force field"
sim.build();
sim.run(frame, n_steps);
```

The conceptual recipe is: **classical MD wiring, ab initio force
evaluator**. Nothing in the loop "knows" the forces are quantum
mechanical; the SCF evaluator simply advertises `Capability::Force` and
the driver treats it like any other force-providing component.

## Why ASPC Matters

A naive BOMD implementation runs a full SCF every step starting from a
generic initial guess (say, a Hartree core diagonalization). For a small
molecule like $\mathrm{H}_2$ this converges in two or three iterations
and the cost is dominated by integral construction. For a hundred-atom
system in a polarized basis, however, full SCF can take 15–30 iterations
per step, and a 10 ps trajectory ($2\times10^4$ MD steps) becomes
prohibitive.

The Always Stable Predictor-Corrector (ASPC) of Kolafa exploits the
smoothness of nuclear motion. Because the density matrix changes slowly
between adjacent MD steps, a Lagrange extrapolation of the density from
the last $K$ steps gives a much better SCF guess than a generic one. A
short corrector step refines the extrapolated density, after which the
calculation either accepts the corrected density (cheap) or falls back to
a full SCF (rare).

In Atomiverse this is one method call:

```cpp
scf.predictor(/*order=*/3, /*corrector_steps=*/2,
              /*corrector_tol=*/1e-4, /*fallback=*/true);
```

For production-scale AIMD (water, peptides, transition-metal complexes)
ASPC commonly cuts the SCF iteration count by 30–50%, bringing the per-MD-step
cost within roughly $1.2\times$ a single fresh SCF. For very small systems
the picture flips: each SCF converges so quickly that ASPC's bookkeeping
overhead exceeds the savings. The
an ASPC-on vs ASPC-off comparison run prints both costs and the
ratio so you can see the trade-off on your machine for whatever system
you are studying.

## Choosing an Ensemble

Atomiverse ships every standard MD ensemble component, and each one
composes with the SCF evaluator with a single `add()` call. The energy
balance changes; the wiring does not.

**NVE — microcanonical.** Drop the thermostat entirely. Total energy
$E = E_{\mathrm{pot}} + E_{\mathrm{kin}}$ is conserved up to integrator
error. Energy drift is the canary for integration quality:

```cpp
sim.add(cpu::md::VelocityVerlet{dt});
sim.on(md::Phase::force).add(std::move(scf));
```

**NVT (Nose-Hoover) — canonical, deterministic.** A single thermostat
chain coupled to the velocities. Best for equilibrium-property runs
where you want a smooth time series.

```cpp
sim.add(cpu::md::NoseHoover{T_target, tau, dt});
```

**NVT (Langevin) — canonical, stochastic.** Adds friction and a random
kick at each step, satisfying the fluctuation-dissipation theorem.
Better for getting out of metastable wells; trajectories from different
seeds are independent samples.

```cpp
sim.add(cpu::md::FixLangevin{T_target, damp, dt, seed});
```

**NPT (Berendsen barostat).** Couples the box edge length to a target
pressure. Combine with a thermostat for proper NPT.

```cpp
sim.add(cpu::md::NoseHoover{T_target, tau_T, dt});
sim.add(cpu::md::Barostat{P_target, tau_P, dt});
```

Every one of these is independently tested against the SCF evaluator
in `tests/cpu/scf/test_scf_evaluator.cpp` (the `BOMDTest` cases).

## Reading The Output

The tutorial example prints two summary lines:

```text
[aimd] step=1   PE_eV=  -30.3887  KE_eV=  0.0129  T_K=  99.61
[aimd] step=100 PE_eV=  -30.1052  KE_eV=  0.0242  T_K= 186.94
```

`PE_eV` is the SCF total electronic energy plus nuclear-nuclear
repulsion, in eV (the same number you would read from `frame.energy`).
`KE_eV` is the classical kinetic energy of the nuclei, computed by the
auto-added `KineticEnergy` component at `end_of_step`. `T_K` is the
instantaneous kinetic temperature derived from $T = 2 E_{\mathrm{kin}} /
(N_{\mathrm{dof}} k_B)$ with $N_{\mathrm{dof}} = 3 N - 3$ (centre-of-mass
removed).

One pitfall is worth flagging up front:

**Statistical fluctuations are huge for small systems.** $\mathrm{H}_2$
has only three vibrational/rotational degrees of freedom after
centre-of-mass removal. The instantaneous temperature can swing
between $50$ K and $700$ K even when the ensemble average is
$100$ K. This is correct equilibrium thermodynamics for a small
system; it is *not* a bug in the thermostat. Trust ensemble
averages over many trajectories, not a single short one.

A historical pitfall — **wrap × non-minimum-image SCF** — is worth
mentioning because it shaped the design of `SCFEvaluator::launch`:
`VelocityVerlet` wraps atoms through periodic boundaries, but the
SCF integral builder is not minimum-image aware. If an atom on the
box boundary wrapped across, the next SCF would see the two atoms
~one box length apart and collapse onto the unphysical RHF ionic
limit. To eliminate this trap, `SCFEvaluator::launch` reassembles
the molecule via minimum-image relative to atom 0 before each
integral build, then restores the wrapped positions before
returning. The user does not have to centre atoms in the box; any
geometry inside the cell is acceptable. (This anchor-on-atom-0
scheme assumes a single isolated fragment per SCF call, which
matches the current scope of CPU SCF.)

## Why Not RHF At High Temperature

RHF cannot describe homolytic bond cleavage. The closed-shell ansatz
forces both electrons to share one spatial orbital, which at large
internuclear separation gives an unphysical mixture of $\mathrm{H}^+ +
\mathrm{H}^-$ and the ordinary covalent dissociation product. The
energy along the dissociation curve is therefore wrong — too high near
the asymptote — and the forces have the wrong sign in some intermediate
range.

In a hot AIMD trajectory this manifests as the molecule appearing to
"fly apart" once thermal motion samples the unphysical region. The
example uses $T = 100$ K, well below the energy where the RHF curve
goes wrong, so the molecule oscillates around its equilibrium bond
length and stays bound. UKS with a GGA functional (the recipe of
Tutorial 4) handles dissociation properly and would let you run hot
AIMD on the same code path — just swap the SCF configuration.

## What Has Been Validated

The Phase 5 spec
([`.claude/specs/phase5-aimd.md`](../../.claude/specs/phase5-aimd.md))
lists the AIMD-specific regressions covered by the test suite:

- **NVE conservation laws** — energy, momentum, and bond oscillation on
  $\mathrm{H}_2$, exercised in `test_scf_evaluator.cpp::BOMDTest`.
- **NVT / NPT composition** — the SCF evaluator coexists with
  Nose-Hoover, Langevin, and Berendsen barostat without conflict.
- **Long-trajectory drift** — 200-step NVE on $\mathrm{H}_2$ keeps total
  energy bounded; the threshold is much tighter than the production-scale
  $10^{-4}$ eV/atom/ps target. See
  `tests/cpu/scf/test_aimd_conservation.cpp`.
- **State-leakage check** — at any point in an ASPC-warmstarted run, a
  fresh stand-alone single-point on the snapshot geometry agrees with
  the propagation force to $\sim 10^{-5}$ eV/Å. This catches regressions
  where the predictor accidentally caches geometry-dependent state and
  forgets to invalidate it.

Long water-box trajectories and pair-correlation comparisons against
CPMD or CP2K references are out of scope for the in-tree tests; they
belong in a separate benchmark that runs offline. The composition tests
above are sufficient to certify that AIMD wiring is correct.

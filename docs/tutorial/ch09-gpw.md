# Chapter 9: Gaussian and Plane Waves (GPW)

The manual has now presented the two great basis families of electronic
structure. The molecular chapters ([Chapter 2](ch02-rhf.md) through
[Chapter 6](ch06-hybrid-dft.md)) used **atom-centered Gaussians** — compact,
few functions per atom, but the two-electron Coulomb term costs formally
$\mathcal{O}(N^4)$ integrals. [Chapter 8](ch08-plane-wave.md) used **plane waves** —
a uniform basis that makes the Hartree potential an $\mathcal{O}(N \log N)$ FFT,
but needs hundreds of functions per atom.

**GPW — Gaussian and Plane Waves** — is the method that takes the best of each.
It is the electronic-structure engine behind CP2K's *Quickstep*, and Atomiverse
reproduces its numbers for the same input.

> **The GPW idea.** Keep the *wavefunction* in a small Gaussian basis, but
> represent the *density* on a uniform real-space grid. The Kohn-Sham matrix is
> built in the Gaussian basis (like [Chapter 5](ch05-kohn-sham-dft.md)); the Hartree and
> exchange-correlation potentials are evaluated on the grid by FFT (like
> [Chapter 8](ch08-plane-wave.md)). A dual representation, bridged by mapping the Gaussian
> density onto the grid and back.

You will learn:

1. **The concept** — the dual Gaussian/grid representation and why it scales
   better than pure Gaussians for the electrostatics.
2. **The representation** — a `gpw_solver` inner engine behind a
   `GPWEvaluator` component, the third member of the evaluator family.
3. **The C++ API** — `GpwElement`, `GpwEvaluatorConfig`, and `single_point`,
   walked through `examples/09_gpw.cpp` (an isolated water molecule).
4. **The expected result** — CP2K Quickstep parity, to the sub-µHartree.
5. **The validation** — against CP2K, and the honest experimental scope.
6. **The recurring pattern** — a third evaluator that drops into the same
   Evaluate / MD / Minimize drivers as the Gaussian and plane-wave ones.

GPW needs FFTW3 and is off by default; build with
`-DATV_BUILD_GPW=ON` (alongside `-DATV_BUILD_RHF=ON -DATV_BUILD_PW=ON`).

## The Concept: Two Representations, One Density

The Kohn-Sham energy has terms that live most naturally in *different* bases.

The **kinetic energy** and the **short-range pseudopotential** are cheap and
accurate in a small Gaussian basis — a handful of contracted functions per atom,
exactly as in the molecular DFT of [Chapter 5](ch05-kohn-sham-dft.md). But the **Hartree
term** — the classical electron-electron Coulomb repulsion — is the expensive
one for Gaussians, because it couples every basis-function pair to every other.

Plane waves invert this. The Hartree potential is trivial in reciprocal space:
solve Poisson's equation by dividing the density's Fourier coefficients by
$G^2$, an FFT away. GPW therefore does the following each SCF iteration:

1. Build the density matrix in the **Gaussian** basis.
2. **Collocate** it onto a uniform real-space **grid** — map the compact
   Gaussian density $\rho(\mathbf r)$ onto grid points.
3. FFT to reciprocal space, solve Poisson for the Hartree potential, evaluate
   the XC potential on the grid — all $\mathcal{O}(N \log N)$.
4. **Integrate** the resulting potential back against the Gaussian basis
   functions to form the potential matrix contribution.
5. Add kinetic + pseudopotential (Gaussian side), diagonalize, rebuild the
   density — the familiar SCF loop.

Core electrons are removed by **GTH (Goedecker-Teter-Hutter) pseudopotentials**,
the same analytic norm-conserving form used in the plane-wave chapter, which
keeps the required grid cutoff modest. The "plane-wave cutoff" here is the
resolution of that real-space grid, set in Rydberg.

The result: Gaussian compactness for the wavefunction, plane-wave-grid scaling
for the electrostatics. This is why Quickstep can run condensed-phase DFT and
AIMD on thousands of atoms.

## The Representation: An Inner Solver Behind an Evaluator

GPW follows the exact architecture of the other electronic-structure methods in
this manual. The SCF machinery is an **inner solver** — `cpu::gpw::gpw_scf`,
driven by `cpu::gpw::GPWEvaluator` — and it is *not* a `Driver` subclass, just
as [Gaussian SCF is an inner solver](ch02-rhf.md) behind `SCFEvaluator` and
[plane-wave DFT](ch08-plane-wave.md) is one behind `PWEvaluator`.

`GPWEvaluator` reports `role() == Evaluator` and `capabilities() == Energy`
(plus `Force` when finite-difference forces are enabled). Its `build()` records
the per-element basis and pseudopotential templates; each `launch()` assembles
the atoms at the frame's current geometry (converting Å → Bohr), runs the GPW
SCF, and writes the total energy to `frame.scalars.energy` in eV — the same
contract every evaluator honours. Because it is an ordinary Energy evaluator, it
plugs into the `Evaluate`, `MD`, and `Minimize` drivers with no special support.

## The C++ API: Elements, Config, `single_point`

The runnable companion is [`examples/09_gpw.cpp`](https://github.com/MolCrafts/Atomiverse/blob/dev/examples/09_gpw.cpp) —
an isolated water molecule, the electronic structure of CP2K's getting-started
`h2o.inp`. Keep it open beside this text.

**Step 1 — build the Frame.** GPW reads `frame.atoms.type` as the atomic number
$Z$ and needs a cell (`cpu::Box`); an isolated molecule sits in a box large
enough to avoid its images:

```cpp
cpu::Frame f(3);
const int    z[3]  = {8, 1, 1};                    // O, H, H
const double x[3]  = {5.0000, 5.7586, 4.2414};     // Å
const double y[3]  = {5.0000, 5.0000, 5.0000};
const double zc[3] = {5.0000, 5.5043, 5.5043};
for (int i = 0; i < 3; ++i) {
    f.atoms.type[i] = z[i];
    f.atoms.x[i] = x[i]; f.atoms.y[i] = y[i]; f.atoms.z[i] = zc[i];
}
f.box = cpu::Box::cubic(10.0);                      // 10 Å box
```

**Step 2 — load the basis and pseudopotential per element.** GPW consumes
CP2K-format data files (`BASIS_MOLOPT`, `GTH_POTENTIALS`), so a `GpwElement`
pairs a contracted Gaussian basis with a GTH potential, keyed by $Z$:

```cpp
const std::string data = /* path to a CP2K data/ directory */;

cpu::gpw::GpwElement O, H;
O.z = 8;
O.basis = cpu::gpw::parse_gth_basis(data + "/BASIS_MOLOPT", "O", "DZVP-MOLOPT-GTH");
O.pot   = cpu::gpw::parse_gth_potential(data + "/GTH_POTENTIALS", "O", "GTH-PBE-q6");
H.z = 1;
H.basis = cpu::gpw::parse_gth_basis(data + "/BASIS_MOLOPT", "H", "DZVP-MOLOPT-GTH");
H.pot   = cpu::gpw::parse_gth_potential(data + "/GTH_POTENTIALS", "H", "GTH-PBE-q1");
```

The `q6`/`q1` suffixes are the valence electron counts the pseudopotential keeps
(6 for oxygen, 1 for hydrogen) — the core is replaced by the GTH form.

**Step 3 — configure and run.** `GpwEvaluatorConfig` carries the grid cutoff, XC
choice, electrostatic boundary condition, and the element list:

```cpp
cpu::gpw::GpwEvaluatorConfig cfg;
cfg.cutoff_ry  = 400.0;                     // real-space grid resolution [Ry]
cfg.xc         = cpu::gpw::XcKind::PBE;      // PADE (LDA) or PBE (GGA)
cfg.nonperiodic = true;                      // isolated: Martyna-Tuckerman electrostatics
cfg.elements   = {O, H};

cpu::gpw::GPWEvaluator eval(cfg);
cpu::single_point(frame, eval);              // frame.scalars.energy now set [eV]

const double e_ha = eval.last_energy_hartree();          // total energy [Hartree]
const bool   ok   = eval.last_result().converged;
const int    n    = eval.last_result().iters;
```

Set `cfg.nonperiodic = false` for a genuinely periodic system (a crystal or
slab) — then the grid electrostatics are periodic, as in
[`examples/10_gpw_znoedl.cpp`](https://github.com/MolCrafts/Atomiverse/blob/dev/examples/10_gpw_znoedl.cpp), which runs
CP2K's periodic ZnO electric-double-layer slab. Optional knobs: `cfg.mix` (density
mixing), `cfg.max_iter`, `cfg.tol`, `cfg.d3` (Grimme DFT-D3 dispersion, requires
`cfg.d3_param_file`), and `cfg.compute_forces` (finite-difference forces — see
the scope note below).

## The Expected Result: CP2K Quickstep Parity

Because GPW is *the same method* as Quickstep — same basis, same
pseudopotentials, same grid electrostatics — the target is not "close to CP2K"
but "CP2K's number." For the water input above (PBE, DZVP-MOLOPT-GTH,
Martyna-Tuckerman) the CP2K reference total energy is

$$ E_\text{total} = -17.21233714\ \text{Hartree}, $$

and `examples/09_gpw.cpp` prints `|dE|` against it. The `GpwResult` also exposes
the energy decomposition — `e_kinetic`, `e_loc_sr` (short-range local PP),
`e_nonlocal`, `e_es` (total electrostatic), `e_xc` — so you can check parity
term by term, which is how you localise a discrepancy if the total drifts.

This is a **numerical-tolerance** result against a reference *engine* (see the
[benchmark validation policy](https://github.com/MolCrafts/Atomiverse/blob/dev/benchmarks/README.md)): the two codes
implement the same equations, so agreement should be to the convergence
threshold, not merely to chemical accuracy.

## The Validation: Against CP2K

CP2K Quickstep is the reference oracle for GPW, exactly as PySCF is for the
Gaussian SCF chapters and ABINIT/CP2K for the plane-wave chapter. The validation
strategy is CP2K parity on:

- **Isolated molecules** (Martyna-Tuckerman) — the water single point above.
- **Periodic systems** — the ZnO slab of `examples/10_gpw_znoedl.cpp`, matching
  CP2K's `ZnO_MD.inp` electronic structure (DZVP-MOLOPT-SR-GTH, GTH-PBE, 320 Ry).

A dedicated `benchmarks/gpw/` suite that pins these parities under CI is on the
roadmap ([`.claude/specs/gpw-gaussian-plane-wave.md`](https://github.com/MolCrafts/Atomiverse/blob/dev/.claude/specs/gpw-gaussian-plane-wave.md));
until then the two examples are the executable checks, each printing its
deviation from the CP2K reference.

## The Recurring Pattern: A Third Evaluator, Same Drivers

Stand back and count the electronic-structure evaluators this manual now has:

| evaluator | basis | electrostatics | chapter |
|---|---|---|---|
| `SCFEvaluator` | Gaussian | analytic ERIs | [2](ch02-rhf.md)–[6](ch06-hybrid-dft.md) |
| `PWEvaluator` | plane wave | reciprocal-space FFT | [8](ch08-plane-wave.md) |
| `GPWEvaluator` | Gaussian + grid | grid FFT (Poisson) | 9 (this) |

All three are `role() == Evaluator`, `capabilities() == Energy`, and all three
write `frame.scalars.energy`. So all three are interchangeable at the *driver*
level: `cpu::single_point(frame, eval)` for an energy, `on(md::Phase::force)
.add(eval)` for AIMD, `on(minimizer::Phase::force).add(eval)` for geometry
optimization. The physics inside each is completely different; the seam the
driver sees — an Energy-capable component at the force phase — is identical.
This is "swap the evaluator, keep the driver" ([Chapter 15](ch15-cross-method.md))
paying off for a *third* electronic-structure family.

## Scope and Status

GPW is **experimental**. Be precise about what works today:

- **Closed-shell only** — a restricted (spin-paired) formulation. No open-shell
  GPW.
- **Γ-point only** — no k-point sampling; correct for large periodic cells where
  a single k-point suffices, not for small unit cells needing a k-mesh.
- **XC:** PADE (LDA) and PBE (GGA). No hybrid functionals on the grid.
- **Forces are finite-difference**, not analytic: `cfg.compute_forces = true`
  displaces each atom by `cfg.fd_step` Bohr and differences the energy — correct
  but expensive ($\mathcal{O}(N)$ extra SCF solves). Analytic GPW forces are a
  later capability.
- The CP2K-format I/O helpers (`parse_gth_basis`, `parse_gth_potential`) live
  under `cpu/gpw/internal/` and are not yet a frozen public interface.

For the periodic plane-wave method GPW builds on, see
[Chapter 8](ch08-plane-wave.md); for the molecular Gaussian DFT it shares a basis with,
see [Chapter 5](ch05-kohn-sham-dft.md).

## Validation Reference

The two runnable GPW calculations each validate against a committed CP2K
reference energy printed at the end of the run:

- Isolated H₂O (PBE, DZVP-MOLOPT-GTH, Martyna-Tuckerman) vs CP2K
  `-17.21233714` Ha: [`examples/09_gpw.cpp`](https://github.com/MolCrafts/Atomiverse/blob/dev/examples/09_gpw.cpp).
- Periodic ZnO slab (CP2K `ZnO_MD.inp` electronic structure):
  [`examples/10_gpw_znoedl.cpp`](https://github.com/MolCrafts/Atomiverse/blob/dev/examples/10_gpw_znoedl.cpp).

The design, CP2K-parity targets, and planned `benchmarks/gpw/` migration are
tracked in
[`.claude/specs/gpw-gaussian-plane-wave.md`](https://github.com/MolCrafts/Atomiverse/blob/dev/.claude/specs/gpw-gaussian-plane-wave.md)
and its production extension
[`.claude/specs/gpw-edl-aimd.md`](https://github.com/MolCrafts/Atomiverse/blob/dev/.claude/specs/gpw-edl-aimd.md).

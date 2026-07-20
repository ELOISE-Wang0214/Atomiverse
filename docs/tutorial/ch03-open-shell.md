# Chapter 3: Unrestricted Hartree-Fock

Restricted Hartree-Fock is a good first model because it explains molecular
orbitals, density matrices, and self-consistency with the least amount of
machinery. It also fails in exactly the situations that matter most for open
shells and bond breaking. This chapter is about that failure and about the
smallest change that repairs it: allow the $\alpha$ and $\beta$ electrons
to occupy different spatial orbitals.

| | |
|---|---|
| **Example** | [`examples/03_uhf.cpp`](https://github.com/MolCrafts/Atomiverse/blob/dev/examples/03_uhf.cpp) |
| **Results** | [`examples/results/03_uhf_pes.json`](https://github.com/MolCrafts/Atomiverse/blob/dev/examples/results/03_uhf_pes.json) |
| **Reference limit** | \(2\times E(\mathrm{H})=-0.9332\) Ha |

```bash
cmake --build build --target 03_uhf
./build/examples/03_uhf
```

## Live result: H₂ dissociation (STO-3G)

RHF cannot dissociate to two neutral H atoms; UHF can. The scan below is the
actual table printed by `03_uhf` (columns: \(R\), \(E_\mathrm{RHF}\),
\(E_\mathrm{UHF}\), \(\langle S^2\rangle\)).

```molplot preset="molplot" theme="auto" width="100%" aspect="16:9"
mark: {type: line, point: true, tooltip: true, strokeWidth: 2}
data:
  values:
    - {R: 0.50, method: RHF, E: -1.042996}
    - {R: 0.75, method: RHF, E: -1.116151}
    - {R: 1.00, method: RHF, E: -1.066109}
    - {R: 1.25, method: RHF, E: -0.989114}
    - {R: 1.50, method: RHF, E: -0.910874}
    - {R: 1.75, method: RHF, E: -0.841349}
    - {R: 2.00, method: RHF, E: -0.783793}
    - {R: 2.50, method: RHF, E: -0.702944}
    - {R: 3.00, method: RHF, E: -0.656048}
    - {R: 3.50, method: RHF, E: -0.629820}
    - {R: 4.00, method: RHF, E: -0.614870}
    - {R: 0.50, method: UHF, E: -1.042996}
    - {R: 0.75, method: UHF, E: -1.116151}
    - {R: 1.00, method: UHF, E: -1.066109}
    - {R: 1.25, method: UHF, E: -0.994190}
    - {R: 1.50, method: UHF, E: -0.957707}
    - {R: 1.75, method: UHF, E: -0.943080}
    - {R: 2.00, method: UHF, E: -0.937213}
    - {R: 2.50, method: UHF, E: -0.933867}
    - {R: 3.00, method: UHF, E: -0.933285}
    - {R: 3.50, method: UHF, E: -0.933182}
    - {R: 4.00, method: UHF, E: -0.933166}
encoding:
  x: {field: R, type: quantitative, title: "R (Å)", scale: {zero: false}}
  y:
    field: E
    type: quantitative
    title: "E (Hartree)"
    scale: {zero: false, nice: true}
  color: {field: method, type: nominal, title: null}
  tooltip:
    - {field: method, type: nominal}
    - {field: R, type: quantitative, format: ".2f"}
    - {field: E, type: quantitative, format: ".6f"}
```

At \(R=4\) Å, UHF sits at **−0.933166 Ha** — within \(3\times 10^{-5}\) Ha of
the exact one-electron limit \(2\times E(\mathrm{H})=-0.9332\) Ha. RHF is still
**−0.61487 Ha**, more than 0.3 Ha too high. Spin contamination
\(\langle S^2\rangle\) rises from 0 near equilibrium toward 1 as the UHF
solution becomes two separate radicals.

## Why RHF Stops Being Enough

The central RHF assumption is that one spatial orbital is occupied by two
electrons of opposite spin. Near equilibrium, that is often a reasonable
picture. For $\mathrm{H}_2$ at $0.74$ Angstrom, it is exactly the right
place to begin.

The problem appears when the bond is stretched. Far apart, the physically
correct dissociation limit is two separate hydrogen atoms, each carrying one
electron. RHF refuses to describe that situation honestly because it insists
that both electrons remain in the same delocalized spatial orbital. The result
is not a small quantitative error. It is a qualitatively wrong picture of the
electrons.

RHF also has no way to describe a genuine odd-electron system such as the
hydrogen atom. There is simply no second electron available to occupy the same
spatial orbital with opposite spin. That is why the UHF example in this
repository begins with the hydrogen atom: it is the smallest case that forces
you to abandon the restricted ansatz.

## What UHF Changes And What It Keeps

UHF keeps the mean-field Hartree-Fock idea. The electrons still move in an
average field built from the current density, and the calculation is still an
SCF fixed-point problem. What changes is the wavefunction ansatz. Instead of one
orbital coefficient matrix $C$, we now have two:

$$
\varphi_i^\alpha(\mathbf{r}) = \sum_\mu C^\alpha_{\mu i}\chi_\mu(\mathbf{r}),
\qquad
\varphi_j^\beta(\mathbf{r}) = \sum_\mu C^\beta_{\mu j}\chi_\mu(\mathbf{r}).
$$

That one change forces the density and Fock matrices to split by spin:

$$
P^\alpha = C^\alpha_{\mathrm{occ}}\left(C^\alpha_{\mathrm{occ}}\right)^T,
\qquad
P^\beta = C^\beta_{\mathrm{occ}}\left(C^\beta_{\mathrm{occ}}\right)^T,
\qquad
P = P^\alpha + P^\beta.
$$

The Coulomb field depends on the total density because every electron repels
every other electron. Exchange depends only on same-spin density because the
exchange term comes from antisymmetry of the wavefunction. That gives the two
UHF Fock matrices:

$$
F^\alpha = H_{\mathrm{core}} + J(P) - K(P^\alpha),
\qquad
F^\beta = H_{\mathrm{core}} + J(P) - K(P^\beta).
$$

Once those matrices are built, the logic is completely familiar. Diagonalize
each spin channel, occupy the lowest orbitals in that channel, rebuild the two
densities, and iterate until they stop moving.

## Reading The Example

The key change in the code is as small as the key change in the theory:
instead of one occupation pattern, you now choose one for each spin channel.
In the example, that choice appears directly in calls such as

```cpp
auto ints = /* IntegralPack is private — use RHF */;
cpu::scf::UHF uhf(n_alpha, n_beta); /* uhf.tol = tol; */
config.max_iter = max_iter;
cpu::scf::UHF uhf(n_alpha, n_beta);
solver.build(ints);
const auto& result = solver.run(ints);
```

The input you must choose carefully is the spin occupation. `n_alpha` and
`n_beta` are not bookkeeping details. They define the spin state. Their sum
must equal the total number of electrons, and their difference fixes
$M_S = (N_\alpha - N_\beta)/2$.

This is why the first example in the file is the hydrogen atom. It forces the
issue immediately: one electron means one occupied spin channel and one empty
spin channel. The second and third examples then show what happens when the
same unrestricted machinery is applied to $\mathrm{H}_2$ near equilibrium and
far from equilibrium.

An equilibrium closed-shell system such as $\mathrm{H}_2$ at $0.74$
Angstrom should come back with nearly identical $\alpha$ and $\beta$
channels and a spin contamination near zero. A hydrogen atom should come back
with `n_alpha = 1`, `n_beta = 0`, and a total energy near the STO-3G reference
of $-0.4666$ Hartree. A stretched $\mathrm{H}_2$ scan should show the
$\alpha$ and $\beta$ channels separating as the bond length grows.

## Charge and Multiplicity Conventions

**Concept.** An open-shell calculation is specified by two user-visible
integers: the molecular charge $q$ and the spin multiplicity $M$. Everything
the solver needs — how many $\alpha$ and how many $\beta$ electrons to occupy
— follows deterministically from those two numbers. Getting this translation
wrong is the single most common way to converge cleanly to the wrong state,
so Atomiverse makes it an explicit, validated step rather than a convention
you must remember.

**Math.** The electron count folds the charge in at integral-build time,

$$ N_e = \sum_A Z_A - q, $$

and the multiplicity $M = 2S + 1$ fixes the spin split through
$n_\alpha - n_\beta = M - 1$:

$$ n_\alpha = \frac{N_e + M - 1}{2}, \qquad
   n_\beta = \frac{N_e - M + 1}{2}. $$

These are integers only when $N_e$ and $M - 1$ have the same parity (even
$N_e$ pairs with odd $M$, odd $N_e$ with even $M$), and $n_\beta \ge 0$
requires $M \le N_e + 1$. Any violation means the requested state does not
exist.

**Algorithm.** The resolver is a struct you instantiate directly —
construction computes and validates in one step, and the resolved counts feed
the existing config builders:

```cpp
/* integrals private inside RHF/UHF/RKS/UKS */       // H2O+, N_e = 9
cpu::scf::SpinResolution occ(ints.n_electrons, /*multiplicity=*/2);
cpu::scf::UHF uhf(occ.n_alpha, occ.n_beta);  // e.g. (5, 4)
```

Worked table (all rows are validated against PySCF 2.13.1 and
regression-locked — see below):

| System | charge | $M$ | $N_e$ | $n_\alpha$ | $n_\beta$ |
|---|---|---|---|---|---|
| $\mathrm{H_2O}$ | 0 | 1 | 10 | 5 | 5 |
| $\mathrm{OH^-}$ | −1 | 1 | 10 | 5 | 5 |
| $\mathrm{H_2O^+}$ | +1 | 2 | 9 | 5 | 4 |
| O atom | 0 | 3 | 8 | 5 | 3 |

Rejected combinations (the constructor throws `std::invalid_argument`):
$\mathrm{H_2O}(0, 2)$ — even $N_e = 10$ with even $M$ breaks parity;
$\mathrm{H_2O^+}(+1, 1)$ — odd $N_e = 9$ with odd $M$ breaks parity;
$\mathrm{H}(0, 4)$ — $M = 4 > N_e + 1 = 2$ would need $n_\beta = -1$.

One convention contrast worth memorizing: PySCF's `mol.spin` is **not** the
multiplicity — it is $2S = n_\alpha - n_\beta = M - 1$. A PySCF doublet is
`spin=1`; the same state here is `SpinResolution(N_e, /*multiplicity=*/2)`.

**Limitations.** The resolver fixes bookkeeping only. It cannot tell you
*which* multiplicity is the ground state (the O atom is a triplet because of
Hund's rule, not because $M=3$ is the only valid input), and it does not
protect against spin contamination — that diagnostic is the next section.

<!-- Every (charge, M, n_alpha, n_beta) tuple and every rejected combination
     in this section is regression-locked by
     tests/cpu/scf/test_spin_resolution.cpp (test_cpu_scf_spin_resolution):
     WorkedTableMatchesTutorialCh03, OhMinusChargeFoldsIntoElectronCount,
     RejectsInvalidCombinations, FeedsUhfConfigRawCounts,
     FeedsUksConfigRawCounts. Change one side only in lockstep with the
     other. -->

## Spin Contamination: What It Means Here

The price of unrestricted flexibility is that the UHF determinant is usually
not an exact eigenfunction of $\hat S^2$. That is what the spin
contamination number measures. In the usual notation, the contaminated amount is

$$
\Delta\langle \hat S^2 \rangle
= \langle \hat S^2 \rangle_{\mathrm{UHF}} - S(S+1).
$$

Near equilibrium in a closed-shell singlet, this number should be close to
zero. When $\mathrm{H}_2$ is stretched and the unrestricted solution breaks
spin symmetry, the contamination rises. That is not automatically a bug. In the
bond-breaking regime it is part of the mechanism by which UHF escapes the
artificially restricted RHF picture.

The repository example makes this tangible. Near equilibrium, RHF and UHF stay
close because the electrons still behave like a closed-shell pair. As the bond
is stretched, the unrestricted solution lowers the energy by letting the two
spins localize differently. The growth of spin contamination is the visible
trace of that symmetry breaking.

## ΔSCF: Ionization Energies Beyond Koopmans

The previous chapter's Koopmans section (see the Koopmans discussion in
[Chapter 2: Restricted Hartree-Fock](ch02-rhf.md)) estimated the ionization
energy of a closed-shell molecule as $I_K = -\varepsilon_\text{HOMO}$: the
negative of the highest occupied orbital energy. That estimate is free — you
already have the orbital energies once the neutral SCF converges — but it rests
on a physical fiction. It assumes that when you remove an electron, the
remaining electrons do not move. They do. ΔSCF is the honest alternative: solve
the cation self-consistently and take an energy difference,

$$ I_{\Delta\text{SCF}} = E^{+}(N-1) - E(N). $$

Both energies are evaluated at the **same** geometry, basis, and convergence
threshold. That makes this a *vertical* ionization energy — the nuclei are
frozen at the neutral equilibrium while an electron is removed.

### Why ΔSCF Is Always Lower Than Koopmans

The relationship between the two estimates is not empirical; it is a two-line
variational argument. Koopmans' $-\varepsilon_\text{HOMO}$ is exactly the energy
of the $(N-1)$-electron determinant you get by deleting one electron from the
neutral orbitals and **freezing everything else**. The ΔSCF cation energy
$E^{+}$ is the result of variationally optimizing that same $(N-1)$-electron
wavefunction over the same Hilbert space. A variational optimum can only lower
the energy, so $E^{+} \le E_\text{frozen}$, and therefore

$$ I_{\Delta\text{SCF}} \le I_K = -\varepsilon_\text{HOMO}. $$

The gap between them, $I_K - I_{\Delta\text{SCF}} \ge 0$, is the **orbital
relaxation energy**: how much the remaining electrons lower their energy by
rearranging around the newly created hole. For a real molecule it is a fraction
of an eV to a couple of eV, and the ΔSCF value is the better one.

There is one clean teaching case where the two estimates are *exactly* equal.
For He and H₂ in a minimal (STO-3G) basis the cation has a single electron in a
spatial orbital that is fixed by symmetry — there is no second electron to
relax and no other orbital of the right symmetry to mix in. The relaxation
energy is identically zero, and the example prints $I_K = I_{\Delta\text{SCF}}$
to solver precision. Watch that row: it is the boundary condition that makes the
inequality meaningful everywhere else.

### Spin Bookkeeping for the Cation

Removing one electron from a closed-shell neutral turns an $N$-electron singlet
into an $(N-1)$-electron doublet ($M = 2$). Using the charge/multiplicity rules
from the section above, the electron count folds the $+1$ charge in at
integral-build time ($N_e = N - 1$), and the doublet split is

$$ n_\alpha = \frac{N}{2}, \qquad n_\beta = \frac{N}{2} - 1. $$

Concretely: the neutral has $N/2$ doubly-occupied levels; the cation keeps all
$N/2$ $\alpha$ orbitals and drops one $\beta$ electron. For H₂O ($N = 10$) that
is $(n_\alpha, n_\beta) = (5, 4)$; for He and H₂ ($N = 2$) it is $(1, 0)$.

### Warm-Starting the Cation

An open-shell SCF can converge cleanly to the *wrong* state — a hole in a lower
orbital rather than the HOMO, or an excited configuration. The defense is to
warm-start the cation from the neutral density instead of from a bare Hcore
guess. Atomiverse exposes this through `SCFInitialGuess`: seed both spin
channels with half of the converged neutral density,

```cpp
cpu::scf::RHF neutral;
neutral.basis = basis;
neutral.charge = 0;
neutral.tol = tol;
neutral_cfg.max_iter = max_iter;
cpu::scf::SCFSolver neutral_solver{neutral_cfg};
neutral_solver.build(ints_n);
const auto& neutral = neutral_solver.run(ints_n);

cpu::scf::SCFInitialGuess guess;
guess.P_alpha = 0.5 * neutral.P_alpha; // split total restricted density evenly
guess.P_beta  = 0.5 * neutral.P_alpha;
guess.has_alpha = guess.has_beta = true;

cpu::scf::UHF cation(N/2, N/2 - 1);
cation.basis = basis;
cation.charge = 1;
cation.tol = tol;
cpu::scf::SCFSolver cation_solver{cfg};
cation_solver.build(ints_c);
const auto& cation = cation_solver.run(ints_c, &guess);
```

A bare Hcore cold start converges to the same energy on these small systems —
the example and its tests check that the two agree — but density warm-starting
is the habit to build: on larger or near-degenerate cations it is the
difference between the ground-state ion and a plausible-looking excited
solution.

### Monitoring Spin Contamination

Because the cation is unrestricted, check its $\langle \hat S^2 \rangle$ the
same way the previous section did. A clean doublet has

$$ \langle \hat S^2 \rangle = S(S+1) = \tfrac{1}{2}\cdot\tfrac{3}{2} = 0.75, $$

so `SCFResult::spin_contamination` (which reports $\langle \hat S^2 \rangle -
S(S+1)$) should sit at ≈ 0. The example prints $\langle \hat S^2 \rangle =
\text{spin\_contamination} + 0.75$ for every cation; values noticeably above
0.75 warn that the unrestricted determinant has picked up higher-spin character
and the ΔSCF number should be read with caution.

### Vertical Versus Experiment

These are *vertical* ionization energies at a fixed neutral geometry. A measured
*adiabatic* ionization energy also lets the cation relax its geometry, which
lowers $E^{+}$ further and makes the adiabatic value smaller than the vertical
one. So even a well-converged ΔSCF number is only a semi-quantitative match to
an experimental ionization energy in a small basis — useful for trends and for
the pedagogy of relaxation, not for spectroscopic accuracy. Enlarging the basis
(STO-3G → 6-31G) and, eventually, adding correlation (a later chapter) moves the
numbers toward experiment; relaxing the cation geometry would convert vertical
to adiabatic. The worked table for all five systems in both bases is
[`examples/12_delta_scf.cpp`](https://github.com/MolCrafts/Atomiverse/blob/dev/examples/12_delta_scf.cpp).

## Running The Example

The build path is straightforward:

```bash
cmake -B build -DATV_BUILD_CUDA=OFF -DATV_BUILD_HF=ON -DATV_BUILD_EXAMPLES=ON
cmake --build build --target 03_uhf
./build/examples/03_uhf
```

The example is a good companion to this chapter because it places RHF and UHF
side by side instead of treating UHF as a separate world. You can watch the
method collapse back to RHF when the restricted approximation is safe, then
watch it diverge only when that extra freedom becomes physically necessary.

## Where Readers Usually Go Wrong

The most common mistake is giving the wrong spin occupation. If the total
electron count is wrong, everything downstream is wrong. If the spin split is
wrong, the solver may converge cleanly to the wrong state. The second common
mistake is to read nonzero spin contamination as a sign that the calculation
failed. In the bond-breaking example it is the opposite: contamination is the
visible trace of unrestricted symmetry breaking doing the job RHF could not do.

The more subtle mistake is to compare only total energies and ignore the rest of
the result. UHF is a spin-resolved method. If the $\alpha$ and $\beta$
densities, orbital energies, and contamination all look like the restricted
solution, then the unrestricted machinery has not actually changed the physics,
even if the solver was configured for UHF.

The next chapter keeps the spin-resolved structure of UHF but changes the mean
field itself. Instead of explicit Hartree-Fock exchange, it uses an
exchange-correlation functional evaluated on a real-space grid:
[Chapter 5: Unrestricted Kohn-Sham DFT (UKS)](ch05-kohn-sham-dft.md).

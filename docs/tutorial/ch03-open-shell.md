# Chapter 3: Unrestricted Hartree-Fock

Restricted Hartree-Fock is a good first model because it explains molecular
orbitals, density matrices, and self-consistency with the least amount of
machinery. It also fails in exactly the situations that matter most for open
shells and bond breaking. This chapter is about that failure and about the
smallest change that repairs it: allow the $\alpha$ and $\beta$ electrons
to occupy different spatial orbitals.

The example for this chapter is
[`examples/03_uhf.cpp`](../../examples/03_uhf.cpp). That
example is worth reading because it does not talk about UHF in the abstract. It
runs three concrete cases: a hydrogen atom, equilibrium $\mathrm{H}_2$, and a
stretched $\mathrm{H}_2$ bond scan where the restricted approximation finally
breaks.

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
auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
auto result = cpu::scf::uhf_solve(ints, n_alpha, n_beta, tol, max_iter);
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

## Running The Example

The build path is straightforward:

```bash
cmake -B build -DATV_BUILD_CUDA=OFF -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
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
even if the function call name says `uhf_solve`.

The next chapter keeps the spin-resolved structure of UHF but changes the mean
field itself. Instead of explicit Hartree-Fock exchange, it uses an
exchange-correlation functional evaluated on a real-space grid:
[Chapter 5: Unrestricted Kohn-Sham DFT (UKS)](ch05-kohn-sham-dft.md).

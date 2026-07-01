# Chapter 2: Atomic Orbitals to Single-Point RHF

This chapter builds the first real quantum-chemical result in the manual — the
total energy of a molecule — and it does so in two deliberate steps. First it
isolates the one idea that makes electronic structure hard: the matrix you want
to diagonalize depends on the eigenvectors that come out of that same
diagonalization, so one solve is never enough. Then it replaces the abstract
version of that loop with the smallest genuine molecular problem in the
repository, closed-shell Restricted Hartree-Fock (RHF) on H₂.

Keep two examples open as you read: [`examples/01_scf.cpp`](../../examples/01_scf.cpp)
for the bare self-consistent loop, and [`examples/02_rhf.cpp`](../../examples/02_rhf.cpp)
for the real RHF calculation the loop turns into.

## Why SCF Has to Iterate

Suppose there is a fixed symmetric matrix $H_0$, a current density-like matrix
$P$, and a scalar coupling strength $\alpha$. A stripped-down model of the
self-consistent field (SCF) problem says the matrix we solve at iteration $k$ is

$$ H^{(k)} = H_0 + \alpha P^{(k)}. $$

Once $H^{(k)}$ is available, we diagonalize it,

$$ H^{(k)} C^{(k)} = C^{(k)} \varepsilon^{(k)}, $$

take the occupied columns of $C^{(k)}$, and rebuild the next density as

$$ P^{(k+1)} = C_{\mathrm{occ}}^{(k)} \left(C_{\mathrm{occ}}^{(k)}\right)^T. $$

Nothing in this loop is chemically realistic yet, and that is intentional — it
makes the fixed-point structure obvious. The current density changes the next
matrix, the next matrix changes the occupied subspace, and the occupied subspace
changes the next density. The calculation stops only when the density hardly
moves anymore:

$$ \left\|P^{(k+1)} - P^{(k)}\right\|_F < \mathrm{tol}. $$

If you keep only one picture from the first half of this chapter, keep this one:
**SCF is an attempt to make the matrix used for the solve consistent with the
orbitals that the solve itself produces.**

## Reading the Bare Loop

[`examples/01_scf.cpp`](../../examples/01_scf.cpp) is deliberately stripped down
— a realistic molecular example would hide the SCF pattern under too many new
ideas at once. It builds a fixed matrix $H_0$, sets the first density guess to
zero, and then runs exactly the equations above: form $H^{(k)}$, diagonalize,
keep the occupied columns, rebuild the density. The most important object is not
$H_0$ but the projector $P = C_{\mathrm{occ}} C_{\mathrm{occ}}^T$ — the
mathematical record of which states are occupied, and the thing fed back into
the next iteration.

The example prints an energy-like scalar $E = \mathrm{Tr}(H P)$ and the density
residual at every step; here those are only diagnostics that let you watch the
iteration settle. Build and run it:

```bash
cmake -B build -DATV_BUILD_EXAMPLES=ON
cmake --build build --target 01_scf
./build/examples/01_scf
```

The residual should shrink until the loop stops changing. That is the whole
lesson of the first half: SCF is not a single diagonalization, it is a search for
a density that reproduces itself.

## What Changes When the Problem Becomes Molecular

Now keep the loop and replace the introductory matrix with real physics: H₂ in
the minimal STO-3G basis — tiny enough to inspect by eye, rich enough to
introduce every central object in closed-shell Hartree-Fock.

We are no longer solving for eigenvectors of an arbitrary symmetric matrix. We
solve for molecular orbitals expanded in a finite atomic-orbital basis
$\{\chi_\mu\}$. Because those basis functions overlap, the eigenproblem is not
the ordinary $FC = C\varepsilon$ form. It is the **generalized** eigenproblem

$$ F C = S C \varepsilon, $$

where $S$ is the overlap matrix — a measure of how far the basis is from being
orthonormal. The matrix $H_{\mathrm{core}}$ collects the one-electron physics
(kinetic energy plus electron-nucleus attraction). The four-index tensor
$(\mu\nu|\lambda\sigma)$ stores the electron-electron repulsion integrals. These
fixed quantities come entirely from the geometry and the basis set; once they
are available, the SCF loop can begin.

Closed-shell RHF assumes electrons come in pairs sharing the same spatial
orbital. For H₂ that is exactly right: two electrons, one occupied spatial
orbital, and a density matrix

$$ P = 2\, C_{\mathrm{occ}} C_{\mathrm{occ}}^T. $$

The factor of $2$ is the signature of the closed-shell model — each occupied
spatial orbital carries one $\alpha$ and one $\beta$ electron. The effective
matrix is now the Fock matrix,

$$ F_{\mu\nu} = H_{\mathrm{core},\mu\nu} +
\sum_{\lambda\sigma} P_{\lambda\sigma}
\left[(\mu\nu|\lambda\sigma) - \tfrac{1}{2}(\mu\lambda|\nu\sigma)\right]. $$

This is the RHF version of the introductory update $H = H_0 + \alpha P$. The
form is more elaborate, but the logic is identical: the current density
determines the next Fock matrix, which determines the next occupied orbitals,
which determine the next density.

## Reading the RHF Example

The setup in [`examples/02_rhf.cpp`](../../examples/02_rhf.cpp) is intentionally
short:

```cpp
cpu::Frame frame(2);
frame.atoms.x[0] = Real(0);     frame.atoms.y[0] = Real(0); frame.atoms.z[0] = Real(0);
frame.atoms.x[1] = Real(0.74);  frame.atoms.y[1] = Real(0); frame.atoms.z[1] = Real(0);
frame.atoms.type[0] = 1;
frame.atoms.type[1] = 1;
frame.box = Mat3::cubic(100.0);
```

Those lines carry more chemistry than they appear to: two nuclei fixed at the
equilibrium bond length, each of charge $+1$, so two electrons to place into
molecular orbitals. The rest of the file follows the RHF logic directly:

```cpp
auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
auto result    = cpu::scf::rhf_solve(integrals, tol, 100);   // result.energy [Hartree]
```

> **Element coverage:** the STO-3G integral path supports elements **H through O
> (Z = 1–8)**. Atoms with higher atomic number throw at `build_integrals` rather
> than returning a wrong answer (there are "no fallback values" by design).
> Heavier-element basis sets are a planned extension.

The first line prepares the fixed molecular quantities (overlap, one-electron
terms, electron-electron repulsion, nuclear repulsion). The second starts the
SCF cycle — exactly the loop from the first half of this chapter, except the
effective matrix is now the Fock matrix and the energy is a real observable:

$$ E_{\mathrm{RHF}} =
\tfrac{1}{2}\mathrm{Tr}\!\left[P(H_{\mathrm{core}} + F)\right] + V_{\mathrm{nuc}}. $$

If you want the same lesson expressed on CUDA, the companion file is
[`examples/02_rhf.cu`](../../examples/02_rhf.cu), but the CPU example is the
place to learn the method itself.

## The Result That Convinces You

At $R = 0.74$ Å in STO-3G, RHF should converge to a total energy near
$\mathbf{-1.1167}$ **Hartree**. Small differences in the last digits from
precision settings are normal; the point is that the example connects the
equations you have just read to a real molecular energy that matches
expectation. This is a **numerical-tolerance** result — the reference is the
textbook/PySCF value, and the SCF benchmark suite pins Atomiverse's RHF energies
against PySCF to `1e-7` (see the [validation policy](../../.claude/specs/manual-system.md)).

Build and run:

```bash
cmake -B build -DATV_BUILD_CUDA=OFF -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
cmake --build build --target 02_rhf
./build/examples/02_rhf
```

A useful check after convergence is that the reported electron count stays
consistent with the molecule you asked for: the density must still represent two
electrons in a closed-shell state.

## What to Take Into Chapter 3

The mental shift is now complete. The iteration still has the same shape it had
in the bare loop, but every object in it has become physical: the density is an
electronic density matrix, the effective matrix is a Fock matrix, and the final
scalar is a real molecular energy. The next chapter keeps this machinery but
drops the closed-shell restriction — where RHF stops being enough and
unrestricted Hartree-Fock begins:
[Chapter 3: Unrestricted Hartree-Fock (UHF)](ch03-open-shell.md).

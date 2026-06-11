# Tutorial 2: Restricted Hartree-Fock (RHF)

The first chapter taught the shape of the SCF loop without tying it to any real
physics. This chapter keeps the loop and replaces the introductory matrix with the
smallest genuine molecular problem in the repository: $\mathrm{H}_2$ in the
minimal STO-3G basis. That problem is tiny enough to inspect by eye but rich
enough to introduce every central object in closed-shell Hartree-Fock.

The example to read while you work through this chapter is
[`examples/02_rhf.cpp`](../../examples/02_rhf.cpp). It sets
up $\mathrm{H}_2$ at a bond length of $0.74$ Angstrom and should converge to a
total energy near $-1.1167$ Hartree. If you later want to see the same lesson
expressed on CUDA, the companion file is
[`examples/02_rhf.cu`](../../examples/02_rhf.cu), but the CPU
example is the place to learn the method itself.

## What Changes When The Problem Becomes Molecular

The conceptual jump from the toy chapter is that we are no longer solving for
eigenvectors of an arbitrary symmetric matrix. We are solving for molecular
orbitals expanded in a finite atomic-orbital basis
$\{\chi_\mu\}$. Because those basis functions overlap, the eigenproblem is
not the ordinary $F C = C \varepsilon$ form. It is

$$
F C = S C \varepsilon,
$$

where $S$ is the overlap matrix. In prose, $S$ tells us how far the
basis is from being orthonormal. The matrix $H_{\mathrm{core}}$ collects the
one-electron physics: kinetic energy plus attraction between electrons and
nuclei. The four-index tensor $(\mu\nu|\lambda\sigma)$ stores the
electron-electron repulsion integrals. These fixed quantities come entirely from
the geometry and the basis set. Once they are available, the SCF loop can begin.

Closed-shell RHF assumes that electrons come in pairs and that the two spins in
each pair share the same spatial orbital. For $\mathrm{H}_2$ that is exactly
the right first approximation: two electrons, one occupied spatial orbital, and
a density matrix

$$
P = 2 C_{\mathrm{occ}} C_{\mathrm{occ}}^T.
$$

The factor of $2$ is the mathematical signature of the closed-shell model.
It says that each occupied spatial orbital carries one $\alpha$ electron and
one $\beta$ electron.

The effective matrix is now the Fock matrix

$$
F_{\mu\nu} = H_{\mathrm{core},\mu\nu} +
\sum_{\lambda\sigma}
P_{\lambda\sigma}
\left[
(\mu\nu|\lambda\sigma) - \frac{1}{2}(\mu\lambda|\nu\sigma)
\right].
$$

This is the RHF version of the introductory update $H = H_0 + \alpha P$. The form is
more complicated, but the logic is the same. The current density determines the
next effective matrix; the next effective matrix determines the next occupied
orbitals; the occupied orbitals determine the next density.

## Reading The Example

The setup in [`examples/02_rhf.cpp`](../../examples/02_rhf.cpp)
is intentionally short:

```cpp
cpu::Frame frame(2);
frame.x[0] = Real(0);     frame.y[0] = Real(0); frame.z[0] = Real(0);
frame.x[1] = Real(0.74);  frame.y[1] = Real(0); frame.z[1] = Real(0);
frame.type_id[0] = 1;
frame.type_id[1] = 1;
frame.box = Mat3::cubic(100.0);
```

Those few lines say more chemistry than they first appear to. The two nuclei
are fixed at the equilibrium bond length, each carries charge $+1$, and the
calculation therefore has two electrons to place into molecular orbitals. Once
the geometry is defined, the rest of the file follows the RHF logic directly:

```cpp
auto integrals = cpu::scf::build_integrals_sto3g(frame);
auto result = cpu::scf::rhf_solve(integrals, tol, 100);
```

The first line prepares the fixed molecular quantities: overlap, one-electron
terms, electron-electron repulsion, and nuclear repulsion. The second line
starts the SCF cycle itself. Conceptually, that cycle is exactly the one from
Tutorial 1, except that the effective matrix is now the Fock matrix and the
energy has become a real physical observable. The total energy is

$$
E_{\mathrm{RHF}} =
\frac{1}{2}\mathrm{Tr}\left[P(H_{\mathrm{core}} + F)\right] + V_{\mathrm{nuc}}.
$$

By this point the energy is no longer a pedagogical placeholder. It is the
quantity you would compare against a textbook value.

## What Result Should Convince You The Chapter Worked

At $R = 0.74$ Angstrom in STO-3G, RHF should converge to a total energy near
$-1.1167$ Hartree. Small differences in precision are normal, but the
important thing in this chapter is whether the example connects the equations
you have just read to a real molecular energy that matches expectation.

Build and run it with

```bash
cmake -B build -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
cmake --build build --target 02_rhf
./build/examples/02_rhf
```

A useful check after convergence is that the reported electron count stays
consistent with the molecule you asked for. In prose, the density must still
represent two electrons in a closed-shell state.

The important mental shift from Tutorial 1 is now complete. The iteration still
has the same shape, but every object in it has become physical: the density is
an electronic density matrix, the effective matrix is a Fock matrix, and the
final scalar is a real molecular energy.

The next chapter keeps this machinery but drops the closed-shell restriction.
That is where RHF stops being enough and where unrestricted Hartree-Fock begins:
[Tutorial 3: Unrestricted Hartree-Fock (UHF)](03-uhf.md).

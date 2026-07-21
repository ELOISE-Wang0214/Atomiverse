# Chapter 20. Molecular properties and environment

A converged single-determinant wave function already provides the energy and
orbitals. This chapter extracts observables and places the molecule in an
external field, a dielectric continuum, or a scalar-relativistic one-electron
Hamiltonian. The order follows the physics: one-electron expectation values,
population and electrostatic-potential charges, response to nuclear
displacements and optical fields, electrostatic environments, and scalar
relativity.

---

## 20.1 Multipole moments

For a single determinant, the expectation value of any one-electron operator
\(\hat O\) is

\[
\langle \Phi | \hat O | \Phi \rangle = \operatorname{Tr}(P O).
\]

The electric dipole and traceless quadrupole (Buckingham convention) are
therefore

\[
\boldsymbol\mu = \sum_A Z_A \mathbf R_A - \operatorname{Tr}(P\,\mathbf r),
\qquad
\Theta_{ab} = \tfrac12 \sum q\,(3 r_a r_b - r^2 \delta_{ab}).
\]

Under a translation of the origin by \(\mathbf d\),
\(\boldsymbol\mu' = \boldsymbol\mu - q_{\mathrm{tot}}\mathbf d\). The dipole
of a neutral system is origin-independent; for charged systems the origin must
be reported with the value. This engine uses the **centre of nuclear charge**
and stores it as `kMultipoleOrigin`, so a rigid translation of the molecule
does not change the reported numbers.

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Pople631Gs;
rhf.run(frame);

cpu::MultipoleObserver obs{rhf.result_handle()};
// Or attach at end_of_step of a host driver; results in kDipole, kQuadrupole,
// kMultipoleOrigin
```

One-electron multipole integrals use Obara–Saika shifted-centre formulae
(*J. Chem. Phys.* **84**, 3963 (1986)). The translation identity
\(M^1(\mathbf C') = M^1(\mathbf C) + (\mathbf C-\mathbf C')S\) is checked to
\(10^{-13}\) in unit tests.

Water RHF/6-31G* \(\mu_z\) agrees with PySCF (table at the end of the chapter).
The HF overestimate of the water dipole (\(\sim 2.23\,\mathrm{D}\) versus
experiment \(\sim 1.85\,\mathrm{D}\)) is a methodological bias, not an
implementation defect.

---

## 20.2 Populations and ESP charges

Atomic charges are not observables: every partition is a convention. The three
routes below optimise different goals and should not be treated as
interchangeable numbers.

### Mulliken and Löwdin

Mulliken population analysis splits the overlap population equally between
centres and is notoriously basis-set sensitive. Löwdin first symmetrically
orthogonalises with \(S^{1/2}\) and then sums over atomic blocks:

\[
q_A = Z_A - \sum_{\mu\in A} (S^{1/2} P S^{1/2})_{\mu\mu}.
\]

Symmetric orthogonalisation is closest to the original basis in a least-squares
sense, but the partition remains orbital-space based; with diffuse bases
“which atom owns a function” is poorly defined. Implementation uses a
unit-diagonal AO metric (\(S_{ii}=1\)) so Cartesian d/f engineering norms do
not enter the charges.

```cpp
auto q  = cpu::lowdin_charge(rhf.result_handle(), frame); // [e]
auto sz = cpu::lowdin_spin(rhf.result_handle(), frame);
```

### ESP fitting: MK and CHELPG

Force-field charges are often chosen to **reproduce the electrostatic potential
outside the molecule**. Sampling outside the van der Waals surface,

\[
V(\mathbf r_k)
  = \sum_A \frac{Z_A}{|\mathbf R_A-\mathbf r_k|}
  - \sum_{\mu\nu} P_{\mu\nu}
    \langle\mu|\,1/|\mathbf r-\mathbf r'|\,|\nu\rangle,
\]

one solves a constrained least-squares problem with
\(\sum q = q_{\mathrm{mol}}\) (Lagrange multipliers / KKT).

- **MK** (Besler–Merz–Kollman): nested shells scaled from Bondi radii.
- **CHELPG** (Breneman–Wiberg): cubic grid points in a band outside the vdW
  surface.

```cpp
cpu::MKCharges mk;
auto q_mk = mk.fit(rhf.result_handle(), frame);

cpu::ChelpgCharges chelpg;
auto q_chelpg = chelpg.fit(rhf.result_handle(), frame);
```

Buried atoms make \(1/r\) columns nearly linearly dependent. When singular
values collapse, an SVD minimum-norm solution is used: the total charge remains
reliable, but the split among collinear atoms is not unique. Fitted charges on
symmetry-equivalent atoms agree only to the grid discretisation
(\(\sim 10^{-3}\,e\)).

For water RHF/6-31G*, Löwdin oxygen charge is about \(-0.73\,e\); MK and CHELPG
are about \(-0.81\,e\). The first is an orbital partition; the second fits the
external potential.

---

## 20.3 Analytic second derivatives and vibrations

Hellmann–Feynman theory removes the density response from **first** energy
derivatives. **Second** derivatives do not enjoy that cancellation: the orbital
response \(U^x\) must be solved. Coupled-perturbed Hartree–Fock (CPHF; Pople
et al., 1979) gives, for each nuclear coordinate \(x\),

\[
{}^1(A+B)\,U^x = -B^x.
\]

The left-hand operator is the same \({}^1(A+B)\) that appears in orbital
stability analysis and in SOSCF. The right-hand side contains skeleton
derivative Fock and overlap contributions; each piece is checked by finite
differences. The assembled analytic Hessian is first validated against finite
differences of analytic gradients, then against external codes. The finite-
difference Hessian pipeline remains available; both feed the same
`vib::HessianResult`.

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Pople631Gs;
rhf.tol = 1e-12;
cpu::vib::HessianResult H = rhf.hessian(frame); // eV/Å²

// Mass-weight, project translations/rotations, diagonalise → frequencies [cm⁻¹]
// See regressions/scf/h2o_rhf_freq
```

`hessian()` is a method output, parallel to `forces()`; vibrational analysis
belongs to vib components. Analytic Hessians are currently limited to
closed-shell RHF; systems without virtual orbitals (e.g. He/STO-3G) are
rejected. Water RHF/6-31G* harmonic frequencies agree with PySCF at the
\(10^{-4}\,\mathrm{cm}^{-1}\) level (table below).

---

## 20.4 Time-dependent HF and excitation energies

Ground-state expectation values answer “what state the molecule is in.”
Excitation energies come from another response: for an oscillating field of
frequency \(\omega\), the polarisability \(\alpha(\omega)\) has poles at the
excitation energies. Time-dependent Hartree–Fock (TDHF; historically RPA)
writes the self-sustained density oscillation as

\[
\begin{pmatrix} A & B \\ -B & -A \end{pmatrix}
\begin{pmatrix} X \\ Y \end{pmatrix}
= \omega
\begin{pmatrix} X \\ Y \end{pmatrix}.
\]

\(X\) are occupied→virtual excitation amplitudes and \(Y\) de-excitation
amplitudes. The blocks \(A\) and \(B\) share the orbital Hessian kernel with
stability analysis, SOSCF, and CPHF.

For real orbitals the problem reduces to a symmetric eigenproblem

\[
(A-B)^{1/2}(A+B)(A-B)^{1/2} Z = \omega^2 Z.
\]

The **Tamm–Dancoff approximation (TDA)** drops \(B\) (hence \(Y\)), giving
\(A X = \omega X\), which is CIS on an HF reference. Variational structure
implies \(\omega_{\mathrm{RPA}} \le \omega_{\mathrm{TDA}}\) for corresponding
roots.

If the reference is triplet-unstable, \({}^3(A+B)\) has a negative eigenvalue
and RPA \(\omega^2\) can be negative—not a numerical bug, but a signal that
perturbation theory about that reference fails. TDA’s \(A\) remains Hermitian,
so its spectrum stays real. Accordingly `tdhf()` runs a stability check first;
on instability it throws with guidance to use `tda()` or re-converge along the
unstable mode.

Transition dipoles (length gauge) and oscillator strengths are

\[
\boldsymbol\mu_{0n}
  = \sqrt2 \sum_{ia}(X+Y)_{ia}\langle\phi_i|\mathbf r|\phi_a\rangle,
\qquad
f = \tfrac23\,\omega[\mathrm{Ha}]\,|\boldsymbol\mu_{0n}|^2[\mathrm{a.u.}].
\]

Triplet channels are spin-forbidden: \(f=0\) exactly.

Formaldehyde is a standard teaching case: the lowest singlet \(n\to\pi^*\)
(\(^1A_2\)) is dipole-forbidden (\(f=0\)); \(\pi\to\pi^*\) is a strong
absorption polarised along C=O. Single-excitation theories typically overshoot
excitation energies by a few tenths of an eV; they are tools for assignment and
semi-quantitative trends, not for spectroscopic accuracy.

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Pople631Gs;
rhf.tol = 1e-12;
rhf.run(frame); // formaldehyde geometry in regressions/scf/ch2o_tdhf

auto tda  = rhf.tda(3);
auto tdhf = rhf.tdhf(3);
```

---

## 20.5 External fields and point-charge embedding

Electrostatic environments do not introduce a new interaction class; they
modify the **one-electron** Hamiltonian. Point-charge embedding is the QM side
of electrostatic QM/MM; a uniform field couples through the dipole. Two-
electron integrals and the Frame layout are unchanged.

Point charges (positions in Å, charges in \(e\)):

\[
h_{\mu\nu} \leftarrow h_{\mu\nu}
  - \sum_C q_C
    \langle\mu|\,1/|\mathbf r-\mathbf R_C|\,|\nu\rangle,
\]

with classical nucleus–charge and charge–charge terms added to the nuclear
repulsion. The integrals reuse the nuclear-attraction primitives.

Uniform field \(\mathbf F\) (atomic units), origin \(\mathbf R_0\) at the
centre of nuclear charge (same convention as §20.1):

\[
h_{\mu\nu} \leftarrow h_{\mu\nu}
  + \mathbf F\cdot\langle\mu|(\mathbf r-\mathbf R_0)|\nu\rangle,
\qquad
E_{\mathrm{nuc}} \leftarrow E_{\mathrm{nuc}}
  - \sum_A Z_A \mathbf F\cdot(\mathbf R_A-\mathbf R_0).
\]

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::STO3G;
rhf.point_charges = {
    {0.0, 0.0,  3.0,  0.5},
    {0.0, 0.0, -3.0, -0.5},
};
rhf.external_field = {0.0, 0.0, 1e-4};
rhf.run(frame);
// Analytic forces with a non-empty environment are not implemented
```

The first-order identity \(E(\mathbf F)=E(0)-\boldsymbol\mu\cdot\mathbf F\)
provides a cross-check: central differences of the energy with respect to a
small field should recover the analytic dipole. For water/STO-3G and
\(F_z=10^{-4}\) a.u. the two agree to \(10^{-5}\) a.u.

In a static field the energy is unbounded below. With strong fields or very
diffuse bases the SCF can “ionise” along \(-\mathbf F\). Finite-field work
should keep \(|\mathbf F|\sim 10^{-4}\) a.u.

When comparing to PySCF `qmmm.mm_charge`, note that PySCF omits the MM–MM
self-energy by default, while this engine includes it in the nuclear term;
add that constant explicitly if total energies are compared
(`regressions/scf/h2o_embedding_field`).

---

## 20.6 Continuum solvation: C-PCM

The solvent may be treated as a continuum of dielectric constant
\(\varepsilon\), with the solute in a cavity of atom-centred spheres. Apparent
surface charges on the cavity respond to the solute potential and are updated
with the density each SCF iteration—unlike the fixed charges of §20.5.

On a tessera discretisation,

\[
S\mathbf q = -f(\varepsilon)\,\mathbf V,
\qquad
f(\varepsilon)=\frac{\varepsilon-1}{\varepsilon+x},
\qquad
\Delta G_{\mathrm{el}}=\tfrac12\sum_k q_k V_k.
\]

\(x=0\) is C-PCM (Barone–Cossi; PySCF default); \(x=1/2\) is original COSMO
(Klamt). At low \(\varepsilon\) the two can differ by a few per cent; for water
(\(\varepsilon=78.4\)) they nearly coincide. Sphere radii are Bondi radii times
a scale (default 1.2); Lebedev points inside other spheres are discarded. The
surface is discontinuous in nuclear coordinates, so geometry scans may show
kinks; surface smoothing (SWIG) is out of scope.

Only the **electrostatic** \(\Delta G_{\mathrm{el}}\) is included. Cavitation,
dispersion, and repulsion are absent, so the number is not a free energy of
solvation comparable to experiment. Analytic forces with solvent on are not
implemented.

```cpp
cpu::scf::RHF sol;
sol.basis = cpu::scf::BasisSet::Pople631Gs;
sol.solvent_eps = 78.4;
sol.solvent_radii_scale = 1.2;
sol.solvent_x = 0.0; // C-PCM
sol.run(frame);
double dg_kcal = sol.e_solvent() * 627.509474;
```

`solvent_eps = 0` (default) is identical to the gas-phase path. Water
RHF/6-31G* gives \(\Delta G_{\mathrm{el}}\approx -8.0\,\mathrm{kcal\,mol}^{-1}\),
within \(0.1\,\mathrm{kcal\,mol}^{-1}\) of PySCF C-PCM (table;
`regressions/scf/h2o_rhf_631gs_cpcm`).

---

## 20.7 Scalar relativity: x2c-1e

Near heavy nuclei, electrons have large local momentum. Dirac theory contracts
s/p shells and lowers their energy, while d/f shells expand through screening.
All-electron calculations beyond the second period without scalar relativity
systematically raise core energies. Spin–orbit coupling is a separate channel
and is not treated here.

Exact two-component (x2c) methods reduce the four-component Dirac problem to a
positive-energy two-component Hamiltonian (Kutzelnigg–Liu; Peng–Reiher).
**x2c-1e** modifies only the one-electron Hamiltonian: in the contracted AO
basis one builds \(\{S,T,V,W\}\) with

\[
W_{\mu\nu}=\langle\nabla\mu|V|\nabla\nu\rangle,
\]

assembles the kinetically balanced Dirac matrix, extracts the decoupling
matrix \(X\) and renormalisation \(R\), and replaces the SCF core Hamiltonian
by \(h^{\mathrm{X2C}}\). Overlap, two-electron integrals, and the form of the
Fock iteration are unchanged (two-electron relativistic corrections are outside
x2c-1e). The speed of light is \(c=137.03599967994\) a.u., matching PySCF.

| | Effective core potential (ECP) | x2c-1e |
|---|---|---|
| Use | Heavy valence (e.g. def2-ECP) | All-electron bases |
| Mechanism | Relativity parameterised into the potential | Explicit one-electron Hamiltonian |
| Mutual exclusion | Must not be combined with x2c | Must not be combined with ECP |

Relativistic ECPs already include scalar effects; stacking x2c double-counts
them. The engine rejects that combination per element. Analytic forces with
x2c on are not implemented.

Expectation values of untransformed non-relativistic property operators
(dipole, ESP, multipoles) computed with x2c orbitals contain **picture-change
error**. This implementation does not transform property operators; such
numbers are only rough estimates.

```cpp
cpu::scf::RHF x2c;
x2c.basis = cpu::scf::BasisSet::Def2TZVP;
x2c.x2c = true;
x2c.run(ar_frame);
```

For Ar RHF/def2-TZVP the relativistic energy lowering is about
\(1.84\,\mathrm{Ha}\), within \(10^{-5}\,\mathrm{Ha}\) of PySCF `x2c1e()` with
`xuncontract=False`.

---

## External numerical comparison

Values below come from in-tree regressions and precomputed PySCF 2.14.0
references (methodology in the [appendix](appendix-pyscf-oracle.md)). The
build and test suite do not depend on Python.

| System | Quantity | Engine | PySCF | Tol. |
|---|---|---|---|---|
| H₂O RHF/6-31G* | \(\mu_z\) / a.u. | −0.87578587 | −0.87578587 | \(10^{-6}\) |
| H₂O RHF/6-31G* | Löwdin \(q_\mathrm{O}\) / e | −0.7284 | −0.7284 | \(10^{-6}\) |
| H₂O RHF/6-31G* | MK \(q_\mathrm{O}\) / e | −0.8071 | −0.8071 | \(10^{-3}\) |
| H₂O RHF/6-31G* | Harmonic \(\nu\) / cm⁻¹ | 1843.7, 4036.2, 4146.4 | same | 1 |
| CH₂O RHF/6-31G* | TDA \(n\to\pi^*\) / eV | 4.6690 | 4.6690 | \(10^{-4}\) |
| CH₂O RHF/6-31G* | TDHF \(n\to\pi^*\) / eV | 4.4999 | 4.4999 | \(10^{-4}\) |
| CH₂O RHF/6-31G* | TDHF \(\pi\to\pi^*\) \(f\) | 0.1796 | 0.1796 | \(10^{-3}\) |
| H₂O STO-3G + point charges | \(E\) (with MM–MM) / Ha | see reg. | see REFERENCE | \(10^{-8}\) |
| H₂O 6-31G* C-PCM | \(\Delta G_\mathrm{el}\) / kcal mol⁻¹ | −8.02 | −7.96 | 0.1 |
| Ar def2-TZVP x2c-1e | \(E\) / Ha | see reg. | −528.64422317 | \(10^{-5}\) |

### Further reading

- Jensen, *Introduction to Computational Chemistry*, molecular properties, solvation, and relativity.
- Szabo & Ostlund, *Modern Quantum Chemistry*, §3 (expectation values); §7 (time-dependent methods).
- Dreuw & Head-Gordon, *Chem. Rev.* **105**, 4009 (2005) (excitation energies).
- Barone & Cossi, *J. Phys. Chem. A* **102**, 1995 (1998) (C-PCM).
- Peng & Reiher, *Theor. Chem. Acc.* **131**, 1081 (2012) (x2c).

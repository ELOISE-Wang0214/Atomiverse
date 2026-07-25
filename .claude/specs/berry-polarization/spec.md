# berry-polarization — Berry-phase (King–Smith–Vanderbilt) polarization on PW

**Loop:** cp2k-method-loop, wave 3. **Domain:** response. **Backend:** CPU PW.
**Depends on:** `bands-dos-si` (NSCF-on-frozen-Veff surface, `rebuild_k_tables`) — done.
**Held sibling:** `wannier-centers` (shares string/localization surface; do not absorb).

## What

Formal polarization of a periodic insulator via the discretized Berry-phase
(King–Smith & Vanderbilt, PRB 47, 1651 (1993); Resta, RMP 66, 899 (1994)):
for each lattice direction \(i\), strings of \(J\) k-points along \(\mathbf b_i\)
over an \(N_\perp \times N_\perp\) transverse mesh,

\[
\varphi_i = \bigl\langle \operatorname{Im}\ln \textstyle\prod_{j=0}^{J-1}
  \det S^{(j)} \bigr\rangle_\mathrm{strings},\qquad
S^{(j)}_{mn} = \langle u_{m,k_j} | u_{n,k_{j+1}} \rangle
\]

with the closing link \(k_{J-1} \to k_0 + \mathbf b_i\) evaluated through the
periodic-gauge coefficient shift \(c_{n,k_0+\mathbf b}(\mathbf G) =
c_{n,k_0}(\mathbf G + \mathbf b)\). In the plane-wave basis the overlap is an
integer-triple (h,k,l) match between the two k-spheres.

Electronic fractional polarization (spin factor \(f=2\), closed shell):
\(p^\mathrm{el}_i = f\,\varphi_i/2\pi\) (units of \(e\), defined mod \(f\)).
Ionic: \(p^\mathrm{ion}_i = \sum_a Z^\mathrm{val}_a \tau_{a,i}\) (fractional
coords). Total \(p_i = p^\mathrm{ion}_i + p^\mathrm{el}_i\) modulo the quantum.
Cartesian cell dipole \(\boldsymbol\mu = \sum_i p_i \mathbf a_i\) [e·Å].

CP2K twin: Quickstep Berry-phase dipole (`&DFT / PRINT MOMENTS PERIODIC T`),
the same modern-theory-of-polarization observable (Γ-limit in CP2K).

## Public surface (architecture-conform)

`include/atomiverse/cpu/planewave/berry_polarization.h`, mirroring the
`Bands`/`Dos` precedent — concrete class, members as parameters, Frame in,
result struct out; **no** `*Config` bag, **no** factories:

```cpp
cpu::pw::BerryPolarization berry;   // members: k_perp, k_string, directions
auto r = berry.compute(method, frame);
// r.phase_el / r.p_el_frac / r.p_ion_frac / r.p_frac / r.dipole_e_ang
```

Pure free helpers (same header, `gap_from_spectrum` precedent):
`phase_from_link_dets` (arg of normalized det product) and
`wrap_centered(p, quantum)`. String geometry is a method
(`string_kpoints(dir)`) so unit tests need no SCF.

## Internal work

- Extract `run_restricted_scf` (currently file-local in `bands_dos.cpp`) into
  `internal/restricted_scf.{h,cpp}` — second real use (rule: extract on the
  second caller). `bands_dos.cpp` switches to the shared one; behaviour
  identical.
- `berry_polarization.cpp`: NSCF Davidson at string k-points on frozen
  \(V_\mathrm{eff}\) (exactly `prob.n_bands` = n_occ bands), occupied
  coefficient copies, per-link overlap determinants (complex LU, partial
  pivoting), string phase average with 2π branch alignment to string 0.
- Valence charges: `valence_from_z` per atom (the same count the SCF uses).

## Root-cause fixes (iron law — discovered while validating quantization)

The Berry phase is first-order sensitive to eigenvector quality and exposed
two latent bugs in the shared PW solver (both fixed here, both verified by
the restored exact phase quantization + full PW suite rerun):

1. `internal/g_vector.cpp` — per-axis enumeration bounds ignored the
   k-offset of the sphere centre and the obliquity of the reciprocal basis,
   clipping sphere-edge G-vectors (spheres at fractional labels 0.125 vs
   0.875 had different G counts). Exact index windows
   \(-f_i \pm g_\mathrm{max}|\mathbf a_i|/2\pi\); FFT sizing unchanged.
2. `internal/davidson.cpp` — the subspace Hamiltonian was truncated to
   `sum.real()` and diagonalized with the real `jacobi_eigh`; at non-TRIM
   k this corrupts Ritz vectors at first order (eigenvalues only at second
   order). New `internal::jacobi_eigh_hermitian` (complex Jacobi) +
   complex subspace matrix. CUDA band solver has the analogous truncation —
   out of scope (CUDA loop backlog).

## Out of scope

Wannier centers / MLWF localization (held sibling), spin-polarized strings,
metallic smearing (throws, as `Bands`), GPW Berry phase, CUDA, finite-field
response. BaTiO₃ (no Ba/Ti/O UPF in tree; Si carries the pins).

## Case

Bulk Si (diamond, a = 5.43 Å, 2-atom fcc cell, ONCVPSP LDA UPF in-tree,
ecut 15 Ry, SCF k-mesh 2):

1. **Centrosymmetric formal polarization** — total \(p_i\) must sit on the
   symmetry-allowed lattice \(\{0, \tfrac12 Q\}\) (Spaldin, J. Solid State
   Chem. 195, 2 (2012)).
2. **Born effective charge ≈ 0** — displace one sublattice atom by
   \(\delta = 0.02\) fractional along \(\mathbf a_3\); ASR + site symmetry give
   \(Z^* = 0\) exactly for an elemental crystal (Pick–Cohen–Martin, PRB 1,
   910 (1970)). The electronic Berry phase must cancel the bare ionic
   \(\Delta p^\mathrm{ion}_3 = Z^\mathrm{val}\delta = 0.08\) — a sharp sign +
   magnitude test of the string machinery.

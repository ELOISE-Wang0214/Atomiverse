# dftb2-multi-k-forces — analytic Bloch-sum nuclear forces for multi-k DFTB2

status: in_progress
branch: wf/dftb2-multi-k-forces
oracle: dftb+ (family); pins are literature-formula / internal-consistency
        (no DFTB+ binary on this host)

## Summary

`cpu::semiempirical::DFTB2` samples a Monkhorst–Pack k-mesh for the **energy**
(restricted `MioSketch`, landed as `dftb2-kpoints`) but throws on
`with_forces` / `with_stress` as soon as any `k_mesh[d] > 1`
(`src/cpu/semiempirical/dftb2.cpp:1595`). This spec adds the **analytic
multi-k nuclear force** for the restricted MioSketch path:

`F = −∇E` of exactly the lattice-summed, Bloch-phased functional the multi-k
SCC converged, so a central finite difference of the multi-k energy reproduces
it. Multi-k **stress** stays throwing (see *Out of scope*).

## Domain basis

SCC-DFTB2 gradient — Elstner et al., *Phys. Rev. B* **58**, 7260 (1998),
Eq. (15); periodic/Bloch force assembly as in DFTB+ (Hourahine et al.,
*J. Chem. Phys.* **152**, 124101 (2020)).

At Γ the total energy is
\(E = \mathrm{Tr}(\rho H^0) + \tfrac12\sum_{AB}\gamma_{AB}\Delta q_A\Delta q_B
+ E_{\mathrm{rep}}\)
and the fixed-Δq gradient is Hellmann–Feynman + Pulay + ∂γ + ∂E_rep.

For a k-mesh the band term becomes
\(E_{\mathrm{band}}=\sum_k w_k\,\mathrm{Re}\,\mathrm{Tr}[P(k)H^0(k)]\),
with \(H^0(k)_{\mu a,\nu b}=\sum_{\mathbf T}e^{2\pi i\,\kappa\cdot n_{\mathbf T}}
h_{\mu\nu}(\mathbf R_b+\mathbf T-\mathbf R_a)\).

Substituting the Bloch sum and using Hermiticity of \(P(k)\) collapses the
band energy onto **real** image-resolved density blocks

\[
\tilde P_{ij}(a,b,\mathbf T)=\sum_k w_k\,
   \mathrm{Re}\!\left[P(k)_{i_a,j_b}\,e^{-2\pi i\,\kappa\cdot n_{\mathbf T}}\right],
\qquad
\tilde W_{ij}(a,b,\mathbf T)=\sum_k w_k\,
   \mathrm{Re}\!\left[W(k)_{i_a,j_b}\,e^{-2\pi i\,\kappa\cdot n_{\mathbf T}}\right],
\]

with \(P(k)=2\sum_o C_{\cdot o}C^{\dagger}_{o\cdot}\) and
\(W(k)=2\sum_o\varepsilon_o C_{\cdot o}C^{\dagger}_{o\cdot}\). Then, per ordered
pair \(a<b\) and image \(\mathbf T\) with bond \(\mathbf d=\mathbf R_b+\mathbf T-\mathbf R_a\),

\[
\frac{\partial E}{\partial d_\alpha}=
2\sum_{ij}\Big[\tilde P_{ij}\big(\partial_\alpha h_{ij}
+\tfrac12(V_a+V_b)\,\partial_\alpha s_{ij}\big)-\tilde W_{ij}\,\partial_\alpha s_{ij}\Big]
+\Delta q_a\Delta q_b\,\partial_\alpha\gamma_{ab}
+\partial_\alpha E^{\mathrm{rep}}_{ab},
\]

i.e. **exactly the Γ pair form** with \(2P\to\tilde P\), \(2W\to\tilde W\).
At Γ (\(\kappa=0\), one k, \(w=1\)) it reduces to the shipped Γ expression, so
the 1×1×1 limit is an identity, not an approximation.

**Reality of the k-sum.** For a real real-space Hamiltonian \(P(-k)=\overline{P(k)}\)
and \(e^{-i\theta(-k)}=\overline{e^{-i\theta(k)}}\), so a ±k pair contributes
\(w\,\mathrm{Re}(\cdot)\) exactly. The MP mesh in the shipped convention
\(\kappa_d=(2r-N_d-1)/(2N_d)\) is always ±-symmetric, so:

- `k_time_reversal = false` (full mesh): the accumulated imaginary part must
  vanish → **checked, throws**.
- `k_time_reversal = true` (folded): each non-self-paired entry stands for its
  own conjugate partner, so `Re()` is the exact fold and there is nothing left
  to cancel; **self-paired** k (κ_d ∈ {0, ±½}, i.e. Γ from odd meshes) have no
  partner and must be real alone → **checked, throws**.

Self-images (a = a + T, T ≠ 0) have \(\mathbf d=\mathbf T\), independent of any
nuclear coordinate → zero nuclear force (same as the Γ path). They are *not*
zero under strain, which is why stress stays out of scope.

The SCC γ term and the repulsive term are k-independent lattice sums already
implemented for Γ and are reused verbatim; Δq is now the k-integrated Mulliken
population produced by the multi-k SCC.

## Design

- No new public type, no `*Config`, no factory. `DFTB2::with_forces` simply
  stops throwing for multi-k restricted MioSketch. Frame in / Frame out.
- All new code is file-private in the existing anonymous namespace of
  `src/cpu/semiempirical/dftb2.cpp`.
- `SccResult` (file-private) gains the converged Bloch data needed by the
  gradient: `kpts`, `Pk`, `Wk`. It already carries the Γ analogues (`C`, `eps`)
  and the open-shell pair; this is the same pattern, not a new one.
- `scc_energy_restricted_kmesh` already re-diagonalizes once at the converged
  Δq to rebuild `Pk`; that block additionally stores `Wk` from the same
  eigenpairs (no extra diagonalization).
- New file-private
  `analytic_forces_restricted_kmesh(...)` mirrors `analytic_forces_restricted`
  (same pair/image walk, same γ/repulsive terms) with the \(\tilde P/\tilde W\)
  contraction and the imaginary-residue gate.
- Throw messages that stay: they must say **what** is unsupported and **why**.
  - multi-k + `with_stress`
  - multi-k + open shell
  - `ThreeOb` + multi-k

## Files

| File | Change |
| --- | --- |
| `src/cpu/semiempirical/dftb2.cpp` | move `Cplx`/`CplxMatrix`/`KPoint` above `SccResult`; extend `SccResult`; store `Wk`; new `analytic_forces_restricted_kmesh`; re-gate `launch_mio`; sharpen throw messages |
| `include/atomiverse/cpu/semiempirical/dftb2.h` | Doxygen: multi-k is now energy **+ analytic forces**; stress still throws |
| `tests/cpu/semiempirical/test_dftb2.cpp` | replace `KMeshRejectsForces` with the force gates below |
| `regressions/semiempirical/dftb2-multi-k-forces/` | new `reg_dftb2-multi-k-forces.cpp` + `REFERENCE.md` |
| `regressions/CMakeLists.txt` | `atv_add_regression(semiempirical dftb2-multi-k-forces)` |
| `docs/developer/tb-method-queue.md` | tick the slug |

## Tasks

1. Reorder the k-point types above `SccResult`; extend `SccResult` with
   `kpts` / `Pk` / `Wk`; populate them in the converged rebuild.
2. Implement `analytic_forces_restricted_kmesh` incl. the imaginary-residue
   gate.
3. Re-gate `launch_mio`: forces allowed for multi-k restricted; stress and
   open shell keep throwing with explicit messages.
4. Unit tests.
5. Regression + `REFERENCE.md` + CMake registration.
6. Header/queue doc updates.

## Testing

Unit (`tests/cpu/semiempirical/test_dftb2.cpp`, one behaviour per test):

- `KMeshForcesMatchFdChain211` — 2×1×1 mesh, off-eq 1-D C-chain cell,
  analytic F_x vs central FD of the multi-k energy, ≤ 1e-6 eV/Å.
- `KMeshForcesMatchFdChain221` — 2×2×1 mesh, same gate.
- `KMeshForcesSumToZero` — ΣF ≤ 1e-10 eV/Å per Cartesian.
- `KMesh111ForcesMatchGamma` — {1,1,1} reproduces the Γ analytic forces to
  ≤ 1e-10 eV/Å (routing identity).
- `KMeshForcesTimeReversalInvariant` — folded vs unfolded mesh agree to
  ≤ 1e-12 eV/Å (this is the imaginary-cancellation gate on the full mesh).
- `KMeshRejectsStress` / `KMeshRejectsOpenShell` — throws survive.

Regression `regressions/semiempirical/dftb2-multi-k-forces/`: out-of-library
consumer, public API only, off-eq C-chain, 2×1×1 mesh; gates = finite E,
non-trivial |F|, ΣF ≈ 0, analytic-vs-FD identity, and multi-k F ≠ Γ F.
`REFERENCE.md` with oracle block, geometry, pins, citation, `Gaps`.

## Landed beyond the plan — multi-k stress + the Γ self-image virial

The stretch goal shipped because the missing piece turned out to be a real
defect in the surface this slug depends on, not extra scope:

`analytic_forces_restricted` walked ordered pairs `a < b` for the force **and**
the virial. Correct for the force (a self-image bond is **T**, independent of
every nuclear coordinate); **wrong** for the virial (|**T**| scales with the
cell). Invisible to `dftb2-analytic-stress` (molecule in a 20 Å box, every
self-image beyond `r_cut`), and a **14% error** on the periodic C-chain:
σ_xx analytic −0.31716963 vs strain-FD −0.37019871 eV/Å³ — confirmed by a RED
test before the fix.

Both virials (Γ and multi-k) now walk the self-image family — electronic
H⁰/S block (factor 1 per ordered image, whole block), γ_AA (factor ½), and
E_rep (factor ½) — and both meet strain-FD to ≤1e-6 eV/Å³. `DFTB2::with_stress`
no longer throws for multi-k. `dftb2-analytic-stress` is bit-unchanged and its
Gap 2 is marked closed.

## Out of scope

- Open-shell (collinear spin) multi-k — needs an α/β Bloch Fock pair.
- `ThreeOb` multi-k — the 3ob radial path has no Bloch-phased builder.
- Periodic D3 (already an open gap, unrelated).
- Fractional / Fermi occupation (multi-k stays insulator-filled, n_occ per k).

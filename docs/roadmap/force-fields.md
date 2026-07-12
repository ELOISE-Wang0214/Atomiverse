# Force Fields

Legend and caveats: [Roadmap index](index.md).

This is the thinnest area of the engine relative to its infrastructure. The
MD machinery around the potentials (neighbor lists, domain decomposition, halo
exchange, PME, special-bonds exclusions, virial accumulation) is mature and
LAMMPS-shaped; the *catalogue of styles plugged into it* is small — one real
pair style, three bonded styles, no impropers, no many-body potentials.

The Lennard–Jones 12-6 pair style is the only ⚡ in the roadmap: it is
SIMD-vectorized (xsimd, `-march=native`) on CPU, has a device kernel on CUDA,
and both are benchmarked against LAMMPS/Kokkos in
`benchmarks/md/lj/performance/`.

## 1. Classical Bonded Potentials

| Method                             | CPU | CUDA |
| ---------------------------------- | :-: | :--: |
| Harmonic bond potential            |  ✓  |   ✓  |
| Quartic bond potential             |  ○  |   ○  |
| Morse bond potential               |  ○  |   ○  |
| FENE bond potential                |  ○  |   ○  |
| Tabulated bond potential           |  ○  |   ○  |
| Harmonic angle potential           |  ✓  |   ✓  |
| Cosine angle potential             |  ○  |   ○  |
| Restricted-bending angle potential |  ○  |   ○  |
| Urey–Bradley potential             |  ○  |   ○  |
| Tabulated angle potential          |  ○  |   ○  |
| Periodic Fourier torsion           |  ○  |   ○  |
| OPLS torsion                       |  ○  |   ○  |
| Ryckaert–Bellemans torsion         |  ✓  |   ✓  |
| CHARMM torsion                     |  ○  |   ○  |
| Harmonic improper torsion          |  ○  |   ○  |
| Periodic improper torsion          |  ○  |   ○  |
| CVFF improper torsion              |  ○  |   ○  |
| CMAP correction                    |  ○  |   ○  |
| Bond–bond cross term               |  ○  |   ○  |
| Bond–angle cross term              |  ○  |   ○  |
| Angle–angle cross term             |  ○  |   ○  |
| Bond–torsion cross term            |  ○  |   ○  |
| Class-II bonded potential          |  ○  |   ○  |
| Tabulated bonded potential         |  ○  |   ○  |

The torsion style is `DihedralMultiHarmonic`, LAMMPS `dihedral_style
multi/harmonic`: `V = Σₙ Cₙ cosⁿφ`. That is the Ryckaert–Bellemans functional
form, which is why the RB row is ✓ — but it is *not* exposed under the OPLS or
CHARMM parameterizations, and there is **no improper style of any kind**.
Bonded forces are cross-validated against LAMMPS in
`benchmarks/md/pe/validation/`.

## 2. Classical Nonbonded Potentials

| Method                              | CPU | CUDA |
| ----------------------------------- | :-: | :--: |
| Lennard–Jones 12-6 potential        |  ⚡  |   ⚡  |
| Lennard–Jones 9-6 potential|  ⊘  |   ⊘  |
| Mie n-m potential                   |  ○  |   ○  |
| Buckingham potential                |  ○  |   ○  |
| Born–Mayer potential|  ⊘  |   ⊘  |
| Weeks–Chandler–Andersen potential|  ⊘  |   ⊘  |
| Morse nonbonded potential           |  ✓  |   ○  |
| Tabulated pair potential            |  ○  |   ○  |
| Coulomb potential                   |  ○  |   ○  |
| Reaction-field electrostatics       |  ○  |   ○  |
| Shifted-potential electrostatics|  ⊘  |   ⊘  |
| Shifted-force electrostatics|  ⊘  |   ⊘  |
| Damped-shifted-force electrostatics |  ○  |   ○  |
| Wolf summation                      |  ○  |   ○  |
| Direct Ewald summation              |  ✓  |   ✓  |
| Particle Mesh Ewald                 |  ✓  |   ✓  |
| Smooth Particle Mesh Ewald          |  ✓  |   ✓  |
| Particle–particle particle–mesh|  ⊘  |   ⊘  |
| Lennard–Jones PME                   |  ○  |   ○  |
| Fast multipole electrostatics       |  ○  |   ○  |
| Multilevel summation method         |  ○  |   ○  |
| Soft-core Lennard–Jones potential   |  ○  |   ○  |
| Soft-core Coulomb potential         |  ○  |   ○  |

Notes:

- **LJ 12-6** is energy-shifted, honours LAMMPS `special_bonds` 1-2/1-3/1-4
  exclusions, and optionally accumulates the per-pair virial. It carries a
  single ε/σ, **not a per-type mixing table** — multi-component fluids are not
  yet expressible.
- **Morse** exists as a demo style: O(N²) all-pairs, deliberately opting out of
  the neighbor list. Correct, but not a production style — and CPU-only.
- **"Coulomb potential"** is ○ because there is no plain cut-off Coulomb style.
  What exists is `CoulLong`, the real-space Ewald-screened `erfc(αr)/r` half of
  the PME/Ewald split — it is not usable standalone.
- The PME implementation is Essmann **smooth** PME (B-spline charge spreading →
  FFT → influence function), so the "Particle Mesh Ewald" and "Smooth Particle
  Mesh Ewald" rows are the same code. Despite an in-code reference to LAMMPS
  `kspace_style pppm`, this is SPME, not PPPM.
- Reciprocal space is **single-GPU only** on CUDA; the CPU PME replicates the
  grid across MPI ranks rather than distributing the FFT.

## 3. Polarizable Force Fields

| Method                              | CPU | CUDA |
| ----------------------------------- | :-: | :--: |
| Drude oscillator model              |  ○  |   ○  |
| Self-consistent induced dipoles     |  ○  |   ○  |
| Extended-Lagrangian induced dipoles |  ○  |   ○  |
| Thole-damped polarization           |  ○  |   ○  |
| Atomic-multipole electrostatics     |  ○  |   ○  |
| AMOEBA electrostatics               |  ○  |   ○  |
| AMOEBA polarization                 |  ○  |   ○  |
| Fluctuating-charge model            |  ○  |   ○  |
| QEq polarization                    |  ○  |   ○  |
| Charge-equilibration dynamics       |  ○  |   ○  |

## 4. Many-Body and Reactive Potentials

| Method                         | CPU | CUDA |
| ------------------------------ | :-: | :--: |
| Embedded-atom method           |  ✓  |   ○  |
| Modified embedded-atom method  |  ○  |   ○  |
| Finnis–Sinclair potential      |  ○  |   ○  |
| Stillinger–Weber potential     |  ✓  |   ○  |
| Tersoff potential              |  ○  |   ○  |
| Brenner REBO|  ⊘  |   ⊘  |
| AIREBO                         |  ○  |   ○  |
| AIREBO-M                       |  ○  |   ○  |
| EDIP                           |  ○  |   ○  |
| COMB                           |  ○  |   ○  |
| ReaxFF                         |  ○  |   ○  |
| Reactive charge equilibration  |  ○  |   ○  |
| Generic bond-order potential   |  ○  |   ○  |
| Three-body tabulated potential |  ○  |   ○  |
| Many-body tabulated potential  |  ○  |   ○  |

**EAM** (analytic Sutton-Chen, `cpu::EAM`) and **Stillinger-Weber** (`cpu::StillingerWeber`)
are implemented on CPU, single-element and single-rank; both are FD-force-gated and
carry literature regressions (`regressions/md/eam_cu`, `regressions/md/sw_si`). See
[Chapter 16](../tutorial/ch16-many-body-potentials.md).

They **throw under MPI**: a many-body force needs a derived per-atom quantity (EAM's
embedding density) at ghost atoms, and the engine has exactly one forward halo
(`pre_neighbor`) and one reverse halo (`post_force`, forces only) — no mid-force
communication point. Fixing that is a `comm/` change, not a potentials change.

Alloy cross-terms, MEAM, Tersoff and the reactive potentials remain ○ — for those,
the quantum-mechanical route (GPW, plane waves) is still the only option.

## 5. Coarse-Grained Methods

| Method                                     | CPU | CUDA |
| ------------------------------------------ | :-: | :--: |
| MARTINI force field                        |  ○  |   ○  |
| SDK coarse-grained force field             |  ○  |   ○  |
| Dissipative particle dynamics              |  ○  |   ○  |
| Multiparticle collision dynamics           |  ○  |   ○  |
| Brownian dynamics                          |  ○  |   ○  |
| Generic tabulated coarse-grained potential |  ○  |   ○  |
| Iterative Boltzmann inversion potential    |  ○  |   ○  |
| Force-matching coarse-grained potential    |  ○  |   ○  |
| Relative-entropy coarse graining           |  ○  |   ○  |

## 6. Machine-Learning Potentials

| Method                                     | CPU | CUDA |
| ------------------------------------------ | :-: | :--: |
| Behler–Parrinello neural-network potential |  ○  |   ○  |
| ANI potential                              |  ○  |   ○  |
| SchNet|  ⊘  |   ⊘  |
| DimeNet|  ⊘  |   ⊘  |
| DimeNet++                                  |  ○  |   ○  |
| PaiNN                                      |  ○  |   ○  |
| GemNet                                     |  ○  |   ○  |
| NequIP                                     |  ○  |   ○  |
| Allegro                                    |  ○  |   ○  |
| MACE                                       |  ○  |   ○  |
| Deep Potential / DeepMD                    |  ○  |   ○  |
| SOAP-GAP                                   |  ○  |   ○  |
| Atomic cluster expansion                   |  ○  |   ○  |
| Moment tensor potential                    |  ○  |   ○  |
| SNAP                                       |  ○  |   ○  |
| Long-range electrostatic ML potential      |  ○  |   ○  |
| Charge-aware ML potential                  |  ○  |   ○  |
| Spin-dependent ML potential                |  ○  |   ○  |
| Delta-learning potential                   |  ○  |   ○  |
| Ensemble ML potential                      |  ○  |   ○  |

**The LibTorch backend is not a machine-learning potential**, and this table
should not be read as "the torch backend is unfinished ML". `ATV_BUILD_TORCH`
builds a *differentiable classical* force field: LJ and harmonic bond as
`autograd::Function` operators with analytic backward passes w.r.t. positions
**and** the force-field parameters (ε, σ, K, r₀). It is a parameter-fitting
surface — there is no network architecture, no descriptor, no integrator, and
no driver on that backend.

## Superseded (⊘)

| Method | Dominated by |
| ------ | ------------ |
| Lennard–Jones 9-6 | Mie n-m, which contains every (n, m) pair including 9-6 and 12-6 as special cases. |
| Born–Mayer | Buckingham — Born–Mayer *is* Buckingham without the r⁻⁶ attraction. |
| Weeks–Chandler–Andersen | The implemented LJ 12-6: WCA is exactly LJ truncated and shifted at 2^(1/6)σ, expressible today. |
| Shifted-potential, shifted-force electrostatics | Damped-shifted-force (DSF), which is the same idea done correctly — both energy and force continuous at the cutoff. |
| Particle–particle particle–mesh (PPPM) | The implemented smooth PME. Same O(N log N) mesh-Ewald family, same accuracy; no reason to carry a second one. |
| Brenner REBO | AIREBO — REBO plus the LJ and torsion terms it is missing. |
| SchNet, DimeNet | The equivariant generation (PaiNN, NequIP, MACE) and DimeNet++ respectively — strictly better accuracy per parameter and per FLOP. |

Kept ○ deliberately: **Wolf summation** and **reaction-field** are not ⊘ despite
being less accurate than PME — they are O(N), which PME is not, so they win in a
regime PME cannot reach.

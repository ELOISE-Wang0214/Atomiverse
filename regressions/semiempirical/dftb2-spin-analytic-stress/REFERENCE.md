# Reference — collinear spin-DFTB2 analytic cell stress (O₂ triplet, NO doublet)

```text
oracle: dftb+
oracle_version: DFTB+ 24.1 (family oracle) — formula-level only; no DFTB+
  binary exists on this pin host (`module spider dftb` → nothing), so no
  absolute DFTB+ number is claimed and none is invented.
family_oracle: dftb+
```

## Method

Collinear unrestricted SCC-DFTB2 (spin channel of Köhler et al.; second-order
SCC of Elstner et al.):

\[
E = \sum_\sigma \mathrm{Tr}\!\left(P^\sigma H^0\right)
  + \tfrac12\sum_{AB}\gamma_{AB}\,\Delta q_A\Delta q_B
  + \tfrac12\sum_A W_A m_A^2
  + E_{\mathrm{rep}}
\]

with the spin-resolved Fock
\(H^\sigma_{\mu\nu}=H^0_{\mu\nu}+\tfrac12 S_{\mu\nu}
(\varepsilon^\sigma_{A(\mu)}+\varepsilon^\sigma_{A(\nu)})\),
\(\varepsilon^\sigma_A=\sum_B\gamma_{AB}\Delta q_B+\sigma W_A m_A\)
(\(\sigma=+1\) α, \(-1\) β), \(\Delta q_A\) the Mulliken charge fluctuation and
\(m_A=N^\alpha_A-N^\beta_A\) the atomic spin population.

Intensive cell stress (mechanical / DFTB+ / PW–GPW convention):

\[
\sigma_{\alpha\beta}=-\frac{1}{\Omega}\frac{\partial E}{\partial\varepsilon_{\alpha\beta}},
\qquad \text{Voigt } (xx,yy,zz,xy,xz,yz)\ \ [\mathrm{eV/\AA^3}]
\]

evaluated **closed form** at the SCC density as the pair virial
\(\partial E/\partial\varepsilon_{\alpha\beta}
=\sum_{\text{pairs}}(\partial E/\partial R_\alpha)R_\beta\) of the same bond
terms as the analytic nuclear force:

| term | contribution |
| ---- | ------------ |
| Hellmann–Feynman | \(\sum_\sigma P^\sigma\,\partial H^0/\partial R\) |
| Pulay | \(\sum_\sigma\left[P^\sigma\tfrac12(\varepsilon^\sigma_a+\varepsilon^\sigma_b)-W^\sigma\right]\partial S/\partial R\) |
| SCC \(\gamma\) | \(\Delta q_a\Delta q_b\,\partial\gamma_{ab}/\partial R\) on the **total** Mulliken \(\Delta q\) |
| repulsive | \(\partial E_{\mathrm{rep}}/\partial R\) |

with \(W^\sigma_{\mu\nu}=\sum_{o\in\mathrm{occ}\,\sigma}
\varepsilon^\sigma_o C^\sigma_{\mu o}C^\sigma_{\nu o}\).

**The spin-polarization virial does not vanish.**
\(E_{\mathrm{spin}}=\tfrac12\sum_A W_A m_A^2\) has no *explicit* strain
dependence — \(W_A\) is a geometry-independent atomic constant and \(m_A\) is
held at its SCC value — but \(m_A\) depends on the cell through \(S\). Because
the SCC functional is variational (\(\partial E/\partial P^\sigma=H^\sigma\),
\(C^{\sigma\mathsf T}SC^\sigma=I\Rightarrow
\mathrm{Tr}(\partial P^\sigma/\partial\varepsilon\;H^\sigma)
=-\mathrm{Tr}(W^\sigma\,\partial S/\partial\varepsilon)\)), that implicit
response reduces **exactly** to the \(\pm W_A m_A\) half of
\(\varepsilon^\sigma_A\) multiplying \(\partial S/\partial\varepsilon\) in the
Pulay row above. It is included; it is not assumed to be zero.

This is the DFTB+ collinear spin + stress module boundary. It replaces the
multi-SCC affine strain-FD product path (12 extra SCC solves) that the
open-shell branch previously used; strain FD survives only as the **test**
oracle.

## Hamiltonian / tables

- ATV class: `cpu::semiempirical::DFTB2` (`sk_set = MioSketch`, the default)
- SK / Hubbard / spin_W / repulsive: static
  `basis_data/dftb_mio_tables.h` (mio-**style** sketch H/C/N/O; \(W_\mathrm{O}
  = -1.00\) eV, \(W_\mathrm{N} = -0.90\) eV). **Not** bit-identical to
  DFTB.org mio-1-1, so no absolute DFTB+ energy/stress identity is claimed.

## Geometry (Å)

Both molecules use the same tilted bond direction
\(\hat u = (0.6,\,0.7,\,\sqrt{0.15})\) — a unit vector, so \(|R| = r\) and all
six Voigt components are non-trivial (no accidental zeros to hide a sign
error).

**Case 1 — O₂ triplet (headline).** `charge = 0`, `multiplicity = 3`
(\(n_\alpha=7,\ n_\beta=5\)), \(r = 1.30\) Å:

| atom | Z | x | y | z |
| ---- | - | - | - | - |
| O | 8 | 0.00 | 0.00 | 0.00 |
| O | 8 | 0.78 | 0.91 | 0.503487835 |

**Case 2 — NO doublet (discriminator).** `charge = 0`, `multiplicity = 2`
(\(n_\alpha=6,\ n_\beta=5\)), \(r = 1.20\) Å:

| atom | Z | x | y | z |
| ---- | - | - | - | - |
| O | 8 | 0.00 | 0.00 | 0.00 |
| N | 7 | 0.72 | 0.84 | 0.464758002 |

Both isolated in a **cubic 20 Å** box (\(\Omega = 8000\) Å³; every lattice
image is beyond the 4.5 Å element cutoff, so this is the molecular limit).
`scc_tol = 1e-12`, `max_scc = 200`, `with_forces = true`, `with_stress = true`.

## Pins (hard-coded gates)

Per case:

1. **Energy:** finite, \(|E| > 10^{-3}\) eV.
2. **Translation:** \(|\sum_i F_i|\le 10^{-6}\) eV/Å per Cartesian.
3. **Stress:** all six Voigt components finite.
4. **Strain-FD identity:** analytic \(\sigma_v\) vs central affine strain-FD of
   the same converged \(E\) (\(h=2\times10^{-4}\)), **all six** components,
   within **1×10⁻⁸ eV/Å³**. A converged SCC-DFTB2 gradient is stationary and
   the virial is the pair contraction of the same bond terms, so this is
   FD-exact to the central-difference truncation error.

Measured 2026-07-29 (`wf/dftb2-spin-analytic-stress`, 1 thread):

```text
O2 triplet   E = -141.69354377 eV
  sigma[0] analytic=-3.81157017e-04  fd=-3.81157010e-04  |d|=7.149e-12
  sigma[1] analytic=-5.18797051e-04  fd=-5.18797024e-04  |d|=2.682e-11
  sigma[2] analytic=-1.58815424e-04  fd=-1.58815425e-04  |d|=8.292e-13
  sigma[3] analytic=-8.89366373e-04  fd=-8.89366187e-04  |d|=1.854e-10
  sigma[4] analytic=-4.92071593e-04  fd=-4.92071573e-04  |d|=1.961e-11
  sigma[5] analytic=-5.74083525e-04  fd=-5.74083491e-04  |d|=3.414e-11
NO doublet   E = -129.78035712 eV
  sigma[0] analytic=-2.65560301e-03  fd=-2.65560297e-03  |d|=3.630e-11
  sigma[1] analytic=-3.61457076e-03  fd=-3.61457072e-03  |d|=4.482e-11
  sigma[2] analytic=-1.10650125e-03  fd=-1.10650124e-03  |d|=9.707e-12
  sigma[3] analytic=-6.19640702e-03  fd=-6.19640706e-03  |d|=4.134e-11
  sigma[4] analytic=-3.42836874e-03  fd=-3.42836869e-03  |d|=5.554e-11
  sigma[5] analytic=-3.99976353e-03  fd=-3.99976346e-03  |d|=7.240e-11
```

Two orders of magnitude inside the gate on every component of both cases.

### Why NO is in the same case, not a second case

O₂ is a **degenerate** probe of the open-shell virial, on two counts, and a
gate built on it alone would pass an implementation that never wrote the spin
term at all:

1. Homonuclear ⇒ \(\Delta q\equiv 0\) (`reg_spin-dftb2-o2` prints
   `dq_Mulliken = [0.000000, 0.000000]`), so the SCC \(\gamma\) virial is
   multiplied by zero.
2. \(m_\mathrm{O}=m_\mathrm{O}=1\) and \(W\) is one element ⇒ the spin shift
   \(\pm W_A m_A\) is **uniform over atoms**. A uniform \(cS\) added to the
   generalized eigenproblem \(HC=SC\varepsilon\) leaves the eigenvectors
   invariant and shifts \(\varepsilon\to\varepsilon\pm c\), so the O₂ triplet
   density, forces and virial come out numerically identical to the restricted
   ones (measured: analytic σ agrees with the `multiplicity = 1` σ of the same
   geometry to \(\sim10^{-19}\) eV/Å³, and \(E_{\text{triplet}} -
   E_{\text{singlet}} = -1.000000000000\) eV \(= \tfrac12\sum_A W_A m_A^2\)
   exactly, geometry-independent).

NO restores both: \(\Delta q\neq0\), and \(W_\mathrm{N}m_\mathrm{N}\neq
W_\mathrm{O}m_\mathrm{O}\). Discriminating power was verified by construction —
zeroing the \(\pm W_A m_A\) term in the Pulay shift while leaving \(W^\sigma\)
consistent moves the residual from 1e-10/1e-11 to **1.6e-05** (O₂) and
**5.1e-06** (NO) eV/Å³, i.e. three orders **outside** the 1e-8 gate on both.
This also closes the "open-shell \(\Delta q\neq0\) case is owed" item recorded
in `regressions/semiempirical/dftb2-analytic-forces/REFERENCE.md` and
`.../dftb2-spin-analytic-forces/REFERENCE.md`.

## Implementation note — one assembly, not three

`analytic_forces_restricted` and `analytic_forces_unrestricted` were folded
into one spin-resolved `analytic_forces_stress` in
`src/cpu/semiempirical/dftb2.cpp`; the restricted path is its degenerate limit
(\(C^\alpha=C^\beta\), \(\varepsilon^\alpha=\varepsilon^\beta\),
\(n_\alpha=n_\beta=n_{\mathrm{occ}}\), \(m\equiv0\Rightarrow
P^\alpha=P^\beta=P/2\), \(W^\alpha=W^\beta=W/2\)). This is deliberate: the
wrong-sign \(\gamma\) gradient repaired on 2026-07-28 had been written out
three times and was fixed three times. The closed-shell equivalence is pinned
by `DFTB2.SpinAnalyticStressClosedShellLimitMatchesRestricted` (baseline
measured on the parent commit `448dca40`; agreement \(\sim10^{-19}\), gate
1e-10 eV/Å³) and by `reg_dftb2-analytic-stress`, whose numbers are unchanged.

## Gaps (still open)

1. **D3 dispersion is absent from the stress.** `DFTB2::d3` contributes to the
   energy and forces but not to `kStress` on any path. Setting `d3 = true`
   with `with_stress = true` silently yields an electronic+repulsive-only
   stress. A D3 virial is owed (shared gap with
   `regressions/semiempirical/dftb2-analytic-stress`).
2. **No periodic stress case, and the self-image virial is missing.** Both
   molecules sit in a 20 Å box, so only the \(a<b\) pair virial is exercised.
   The pair walk in `analytic_forces_stress`
   (`src/cpu/semiempirical/dftb2.cpp:1325`, `for (int b = a + 1; ...)`) skips
   self-images (\(a=a\), \(T\neq0\)). Those exert **zero nuclear force** (so
   forces are unaffected, and `reg_dftb2-periodic-forces` stays green) but they
   **do** carry a virial once a lattice image falls inside the SK/γ/repulsive
   cutoff. Under dense Γ-PBC (e.g. the 2.8 Å C-chain cell) the analytic stress
   is therefore incomplete. This is **inherited**, not introduced: the shipped
   restricted path had exactly the same pair walk. A periodic strain-FD gate +
   the self-image virial are owed together.
3. **ThreeOb stress is still multi-SCC strain-FD** (`launch_3ob`); only the
   MioSketch path is closed form.
4. **Multi-k stress still throws** (`k_mesh` beyond Γ is energy-only).
5. **No absolute DFTB+ number.** In-tree tables are mio-**style** sketches and
   no DFTB+ binary exists on this host, so every gate here is an
   internal-consistency (analytic vs strain-FD) identity, not an oracle pin.

## Citations

- M. Elstner, D. Porezag, G. Jungnickel, J. Elsner, M. Haugk, Th. Frauenheim,
  S. Suhai, G. Seifert, *Self-consistent-charge density-functional
  tight-binding method for simulations of complex materials properties*,
  Phys. Rev. B **58**, 7260 (1998). DOI 10.1103/PhysRevB.58.7260.
- C. Köhler, G. Seifert, U. Gerstmann, M. Elstner, H. Overhof, Th. Frauenheim,
  *Approximate density-functional calculations of spin densities in large
  molecular systems and complex solids*, Phys. Chem. Chem. Phys. **3**, 5109
  (2001). DOI 10.1039/b105782k. — the collinear spin channel
  \(E_{\mathrm{spin}}=\tfrac12\sum_A W_A m_A^2\) and its \(\pm W_A m_A\)
  Hamiltonian shift.
- C. Köhler, Th. Frauenheim, B. Hourahine, G. Seifert, M. Sternberg,
  *Treatment of collinear and noncollinear electron spin within an approximate
  density functional based method*, J. Phys. Chem. A **111**, 5622 (2007).
  DOI 10.1021/jp068802p.
- B. Hourahine et al., *DFTB+, a software package for efficient approximate
  density functional theory based atomistic simulations*, J. Chem. Phys.
  **152**, 124101 (2020). DOI 10.1063/1.5143190. — DFTB+ 24.1 is the family
  oracle for the force/stress module boundary.

## Perf

```text
perf_exempt: true
perf_exempt_reason: no DFTB+ binary / fair wall on pin host. The analytic path
  replaces 12 multi-SCC energy evaluations (6 Voigt x +/-h) with one post-SCC
  pair contraction.
```

## Units

Public ATV: energy **eV**, positions **Å**, forces **eV/Å**, stress
**eV/Å³** Voigt \((xx,yy,zz,xy,xz,yz)\).

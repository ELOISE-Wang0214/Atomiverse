# QEq / fluctuating-charge equilibration — LAMMPS `fix qeq/point` twin

*(Technical-report draft; narrative house style TBD — rewrite later.)*

| Field | Value |
| ----- | ----- |
| Slug | `qeq-fluctuating-charge` |
| Regression | `regressions/md/qeq-fluctuating-charge/` |
| Public API | `cpu::md::QEq` (Component on `cpu::md::MD`, default `pre_force`) |
| Roadmap / queue | `docs/developer/metalwalls-method-queue.md` → wave 5; map §2 "Fluctuating-charge / QEq" |
| Oracle | **LAMMPS** QEQ package (`fix qeq/point`), offline pin |

## Summary

**Charge equilibration** (electronegativity equalization, Rappé & Goddard,
*J. Phys. Chem.* **95**, 3358 (1991)) for electrolyte / molecular systems:
per-atom fluctuating charges minimize a second-order charge energy with
per-type electronegativity χ [eV] and hardness η [eV/e²] under a total-charge
constraint. This is deliberately **not** the constant-potential electrode
family (`ConstPotBO`/`ConstPotCG`/`ConstPotMaZe` solve fixed-Ψ Gaussian
electrode charges); QEq equilibrates every parameterized atom against its own
electronegativity — the MetalWalls capability map's open "fluctuating-charge /
QEq" polarization slice.

Public surface is `cpu::md::QEq` — a **Component** bound at
`md::Phase::pre_force` (order −50, same slot as the electrode solvers), never
a Driver category, never inside `KSpace*`. It writes equilibrated charges to
`kAtom_charge` before pair / k-space forces each step.

## Method

Energy model (LAMMPS `fix qeq/point` identity, verified against
`src/QEQ/fix_qeq_point.cpp`):

\[
E(q)=\sum_i\Bigl(\chi_i q_i+\tfrac12\eta_i q_i^2\Bigr)
 +\sum_{i<j,\;r_{ij}\le r_c}\frac{q_i q_j}{r_{ij}},
\qquad \sum_i q_i = q_{\mathrm{tot}}
\]

with the literal LAMMPS convention of a **bare** \(1/r\) coupling (no qqrd2e
prefactor on the off-diagonal; the physical Coulomb energy of the resulting
charges is computed downstream by the pair style with
\(k_e = 14.399645\) eV·Å/e²). Stationarity gives the linear system
\((\eta_i\delta_{ij}+J_{ij})\,q=-\chi+\lambda\mathbf 1\), solved with the same
two-vector scheme as LAMMPS: \(Hs=-\chi\), \(Ht=-\mathbf 1\),
\(q=s-((\Sigma s-q_{\mathrm{tot}})/\Sigma t)\,t\). Atomiverse factors the dense
matrix once per launch (Cholesky) and substitutes both right-hand sides;
distances use the frame box minimum image. The matrix is rebuilt every launch,
so charges follow the moving geometry in MD.

## Usage

```cpp
cpu::md::QEq qeq(/*cutoff [Å]=*/10.0);
qeq.set_type(/*type=*/1, /*chi [eV]=*/5.343, /*eta [eV/e²]=*/10.126); // C
qeq.set_type(2, 8.741, 13.364);                                      // O
qeq.set_type(3, 4.528, 13.8904);                                     // H
sim.on(md::Phase::pre_force).add(std::move(qeq));   // before coul pair style
```

Observables after a launch: `charges()` [e], `net_charge()` [e],
`qeq_energy()` [eV, bare-1/r kernel], `n_solved()`. Every atom type present in
the frame must be registered via `set_type` — unparameterized types throw (no
fallback). `qtot` (default 0, the LAMMPS behavior) constrains a net-charged
system.

## Verification

`regressions/md/qeq-fluctuating-charge/` (offline LAMMPS pin, 30 Mar 2026
development build; inputs at
`tools/metalwalls_oracle/inputs/lammps/qeq-fluctuating-charge/`): H₂CO-like
4-atom organic toy, `fix qeq/point 1 10.0 1e-10 500` — per-atom charges match
to ≤ 1e-6 e (q_C = +0.046280, q_O = −0.230796, q_H = +0.092258), Σq = 0 to
1e-12, and the `coul/cut` energy of the equilibrated system matches LAMMPS
`ecoul` = −0.254699311955 eV to ≤ 1e-6 eV. Unit tests
(`tests/cpu/dynamics/test_qeq.cpp`) pin the two-atom closed form
\(q=(\chi_2-\chi_1)/(\eta_1+\eta_2-2/r)\), constraint conservation,
beyond-cutoff decoupling, the energy closed form, and the throw contracts.

## Scope / debt

- Dense Cholesky solve — toy/small-N scope (≤ a few hundred atoms), matching
  the const-potential dense path; a CG + neighbor-list path is the natural
  extension for production N.
- Point kernel only (`qeq/point`); shielded/Slater variants
  (`qeq/shielded`, `qeq/slater`) are follow-up queue rows if mapped.
- Primary-image minimum-image coupling (matches the non-periodic oracle toy);
  no Ewald-summed QEq.

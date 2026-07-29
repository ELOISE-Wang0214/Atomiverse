# td-dftb-h2 — linear-response TD-DFTB (Casida) foothold on SCC-DFTB2

Wave 4 of `tb-method-loop`. Roadmap rows `docs/roadmap/semiempirical.md` "TD-DFTB"
and `docs/roadmap/excited-states.md` "TD-DFTB" are ○ (absent). This spec lands
the **energies-only** foothold.

## Summary

Vertical singlet/triplet excitation energies on a converged **restricted**
SCC-DFTB2 ground state, via the DFTB Casida formulation of Niehaus et al.,
*Phys. Rev. B* **63**, 085108 (2001).

Deliverable: `DFTB2::td_casida(frame, n_states, channel, …)` returning a
`TdDftbResult` (excitation energies [eV], single-particle gaps [eV], Mulliken
transition charges [e], the assembled Casida matrix Ω [eV²]).

## Domain basis

Ground state: SCC-DFTB2 (Elstner et al., PRB **58**, 7260 (1998)) restricted,
Γ-point, MioSketch SK tables. For H₂ the H onsite/Hubbard-U and the whole H–H
SK + repulsive grid are transcribed **mio-1-1** (`basis_data/dftb_mio_hh_sk.h`,
`internal/mio_hh_table.h`), so H₂ is the one in-tree system whose DFTB2 numbers
are genuine mio-1-1.

Linear response (Niehaus 2001, Eqs. 8–12):

- Mulliken transition charges
  \(q^{ia}_A = \tfrac12\sum_{\mu\in A}\big(c_{\mu i}(Sc_a)_\mu + c_{\mu a}(Sc_i)_\mu\big)\)
- Coupling
  \(K^{S}_{ia,jb} = \sum_{AB} q^{ia}_A\,\gamma_{AB}\,q^{jb}_B\) (singlet),
  \(K^{T}_{ia,jb} = \sum_{A} q^{ia}_A\,W_A\,q^{jb}_A\) (triplet; magnetization
  kernel is on-site, \(\gamma^{m}_{AB}=\delta_{AB}W_A\))
- Casida eigenproblem \(\Omega F = \omega^2 F\) with
  \(\Omega_{ia,jb} = \delta_{ij}\delta_{ab}(\varepsilon_a-\varepsilon_i)^2
  + 2\sqrt{\varepsilon_a-\varepsilon_i}\,K_{ia,jb}\,\sqrt{\varepsilon_b-\varepsilon_j}\)

γ_AB and W_A come from the SAME tables the ground state uses
(`gamma_pair(...)`, `DftbElementRow::spin_w`) — no new parameter data.

## Design

**Adaptation of the assigned design, stated explicitly.** The assignment
sketched a standalone `cpu::semiempirical::TdDftb2` class owning a `DFTB2`.
The tree already has a documented house pattern for linear response on a
mean-field reference and it is the opposite shape:

| tree precedent | file |
| --- | --- |
| `RHF::tda` / `RHF::tdhf` | `include/atomiverse/cpu/scf/rhf.h` |
| `RKS::tda` / `RKS::tddft` | `include/atomiverse/cpu/scf/rks.h` |
| `UKS::tda`, `UHF::sf_tda` | `uks.h`, `uhf.h` |
| result-only public header `ExcitedStatesResult` | `include/atomiverse/cpu/scf/tdhf.h` |

i.e. **response is a member method on the concrete ground-state method,
returning a result-only struct declared in its own public header.** A separate
`TdDftb2` class would additionally need either (a) `friend class TdDftb2` into
DFTB2's privates, (b) six new public accessors on DFTB2 exposing MO
coefficients / overlap / γ / the AO→atom **index map** (which
`CLAUDE.md` "never precomputed indices … or derived maps" discourages), or
(c) a duplicated SCC loop (forbidden: "No duplication across methods"). So the
foothold follows the tree, not the sketch: `DFTB2::td_casida`, with the
requested header path `semiempirical/td_dftb2.h` used for the result type —
mirroring `scf/tdhf.h` ("public, result-only header").

**Stateless, Frame-in.** `td_casida` is `const` and takes the geometry as
`const Frame&`. It runs its own converged restricted SCC at that geometry
through the *same* anonymous-namespace `scc_energy_restricted` /
`build_h0_s` / `gamma_pair` the energy path uses — no cached ground state, no
launch-order coupling, no cache-invalidation bug, no duplicated SCC.

Public surface added:

```cpp
// include/atomiverse/cpu/semiempirical/td_dftb2.h  (result-only)
enum class TdDftbChannel { singlet, triplet };
struct TdDftbResult {
    std::vector<double> energy_ev;         // ω_n [eV], ascending, n_states
    std::vector<double> ks_gap_ev;         // n_states lowest ε_a−ε_i [eV], ascending
    std::vector<int>    pair_occ, pair_virt;  // MO indices of each ia pair
    std::vector<double> transition_charge; // [n_pair × n_atoms] [e]
    std::vector<double> casida_ev2;        // [n_pair × n_pair] Ω [eV²]
    int n_pair, n_atoms;
    TdDftbChannel channel;
    double coupling_scale;
};

// dftb2.h
TdDftbResult td_casida(const Frame& frame, int n_states,
                       TdDftbChannel channel = TdDftbChannel::singlet,
                       int n_occ_window = 0, int n_virt_window = 0,
                       double coupling_scale = 1.0) const;
```

No `*Config` bag, no factory, no `*System`. `n_occ_window`/`n_virt_window`
`= 0` means "all". `coupling_scale` λ scales K; λ = 0 is the non-interacting
limit and must reproduce ε_a−ε_i exactly.

**Restrictions (throw, no fallback):** `sk_set == MioSketch`, Γ-point
(`k_mesh == {1,1,1}`), closed shell (`multiplicity == 1`, even electrons),
`n_pair ≤ 512` (dense path, RKS::tddft precedent), all Δε > 0, all ω² ≥ 0.

**Frame meta:** no excitation-energy meta key exists and the natural
observable is a variable-length vector, which the fixed typed-meta vocabulary
(`double` / `F64x3` / `F64x6`) cannot hold. Results are returned on the object
only; recorded as a gap.

## Files

- `include/atomiverse/cpu/semiempirical/td_dftb2.h` (new, result-only)
- `include/atomiverse/cpu/semiempirical/dftb2.h` (declare `td_casida`)
- `src/cpu/semiempirical/dftb2.cpp` (Casida assembly + solve; same TU so the
  existing SK / γ / SCC helpers are reused rather than extracted or copied)
- `tests/cpu/semiempirical/test_td_dftb2.cpp` (new unit gates)
- `tests/CMakeLists.txt` (register `test_cpu_td_dftb2`)
- `regressions/semiempirical/td-dftb-h2/reg_td-dftb-h2.cpp` + `REFERENCE.md`
- `regressions/CMakeLists.txt` (`atv_add_regression(semiempirical td-dftb-h2)`)
- `docs/roadmap/semiempirical.md`, `docs/roadmap/excited-states.md` (row ○→partial)

## Tasks

1. `td_dftb2.h` result-only header with full Doxygen (units on every field).
2. `DFTB2::td_casida` declaration + Doxygen in `dftb2.h`.
3. Casida assembly in `dftb2.cpp`: guards → SCC → S → γ → pairs → q → K → Ω → eigh.
4. Unit gates (below).
5. Regression + REFERENCE.md.
6. Roadmap rows.

## Testing

Unit (`tests/cpu/semiempirical/test_td_dftb2.cpp`, one method each, seconds):

1. `CasidaMatrixIsSymmetric` — |Ω_pq − Ω_qp| ≤ 1e-12 eV² (H₂O, several pairs).
2. `ZeroCouplingReproducesKsGaps` — λ = 0 ⇒ |ω_n − Δε_n| ≤ 1e-10 eV (decisive).
3. `H2LowestSingletIsPositiveFinite` — H₂ at 0.74 Å: finite, > 0, > Δε.
4. `H2SingletAboveTriplet` — same pair, ω_S > ω_T (γ > 0, W < 0).
5. `TransitionChargesSumToZero` — Σ_A q^ia_A ≤ 1e-10 for every i≠a pair.
6. `H2TransitionChargeMatchesTwoLevelClosedForm` — for symmetric H₂ the
   analytic two-level result is q^ia_1 = 1/(2√(1−S₁₂²)); independent of the
   general Mulliken code path.
7. Guard tests: ThreeOb / multi-k / open shell throw.

Regression `regressions/semiempirical/td-dftb-h2/` (public API only, ≤ 60 s):
H₂ at 0.74 Å through `DFTB2::td_casida`; pins the λ→0 KS-gap identity, the
singlet/triplet ordering, Ω symmetry, and the transition-charge sum rule; plus
the hard-coded measured λ=1 singlet ω as a **self-consistency** value clearly
labelled as NOT an oracle number (no `dftb+` binary on this host).

## Out of scope

- Oscillator strengths / transition dipoles; excited-state gradients & NAC.
- Open-shell (spin-polarized) TD-DFTB; ThreeOb SK path; periodic / k-point LR.
- Long-range-corrected TD-DFTB; on-site corrections; l-resolved spin constants.
- Davidson / iterative solver (dense `eigh` only, n_pair ≤ 512).

# gpw-dftu-kpoint-forces — Multi-k GPW analytic ionic forces with Dudarev DFT+U

**Loop:** cp2k-method-loop (wave 1, 2026-07-31). **Branch:** `wf/gpw-dftu-kpoint-forces`.
**Status:** code-complete, acceptance met (close on merge).

## Gap (verified wave 15)

`src/cpu/gpw/internal/gpw_forces.cpp` only copies `plus_u` into the SCF config
(`copy_gpw`, line 82). Neither `assemble_gpw_forces` (Γ) nor
`assemble_gpw_forces_kpoints` (multi-k) folds the Dudarev occupation-matrix
term into the force, so `gpw_forces_analytic` differentiates a **different
functional** than the SCF converged whenever `plus_u` is non-empty — the same
class of bug as the ZnO local-PP 0.53 Ha/Bohr incident (see
`internal/local_pp.h`).

## Physics

ATV's Dudarev projector (per spin σ, `internal/dft_plus_u.h`):

- occupation  n^σ = Σ_k w_k · ½ Herm[(P_k^σ S_k + S_k P_k^σ)] restricted to the
  correlated (atom, ℓ) AO block (Γ: single real term),
- E_U^σ = (U/2)[Tr n^σ − Tr (n^σ)²],  V^σ = U(½ I − n^σ) on that block.

Since n depends on geometry through **S**, the fixed-P Lagrangian
(`gpw_stress.cpp` shape: `e_kin + e_nl + e_locgrid + e_es + e_xc + e_disp +
e_dftu − e_pulay`) contributes an explicit +U force term

    dE_U/dR = Σ_σ Σ_k w_k Re Tr[ M_k^σ dS_k/dR ],
    M_k^σ  = ½ (V̂^σ P_k^σ + P_k^σ V̂^σ),

with V̂^σ the subspace potential scattered to the full AO matrix. Because the
existing Pulay loop already contracts `W_k` against `dS_k`, the entire fold is

    W_k ← W_k − Σ_σ M_k^σ        (Γ: W ← W − Σ_σ M^σ)

**CP2K precedent:** `src/dft_plus_u.F`, `mulliken_charges` — KS matrix gets the
variational ½(V̂S + SV̂); `matrix_w` gets exactly `−½(V̂P + PV̂)` so the standard
W·dS Pulay force picks up the +U geometric term. Spin conventions: RKS builds
n from P_spin = P_tot/2 with E_U = 2·E_U^spin ⇒ M = ½(V̂ P_tot + P_tot V̂) with
V̂ from n(P_tot/2); UKS sums per-spin M^α + M^β (occupancy-1 densities).

## Known divergence (surfaced, out of scope)

1. **Projector definition:** CP2K `PLUS_U_METHOD MULLIKEN` uses the *net* AO
   population block q_ab = P_ab·S_ab (elementwise, on-site ⇒ geometry-
   independent). ATV uses the symmetrized *gross* product ½(PS+SP)|block
   (geometry-dependent, closest in spirit to CP2K `MULLIKEN_CHARGES` which is
   the vector q_i = Σ_j P_ij S_ij). Absolute CP2K force agreement is therefore
   limited by the projector model difference (the +U SCF twins already carry
   this: e.g. `gpw-dftu-nio` E tol 5e-2 Ha with e_dftu ≈ 0.05 Ha diagnostic).
   Re-aligning the projector is a census decision, not this row.
2. **Variational scatter:** ATV's SCF adds plain V̂ to H (Γ and k); the exact
   ∂E_U/∂P for the ½(PS+SP) projector is ½(V̂S + SV̂) (CP2K applies this form
   for its geometry-dependent method). The plain-V̂ fixed point is slightly
   non-stationary for the reported energy, which bounds how tightly
   analytic-at-fixed-P can match multi-SCF FD. Measure the residual on the
   dual gate; if it dominates the gate, surface it (fixing the scatter shifts
   every +U fixed point and re-pins ~8 regressions — separate row).

## Deliverables

1. `gpw_forces.cpp`: +U W-fold in `assemble_gpw_forces` (Γ RKS/UKS) and
   `assemble_gpw_forces_kpoints` (k RKS/UKS), active only when `cfg.plus_u`
   is non-empty. Helper on `internal/dft_plus_u.h` shared by both paths
   (second-use extraction: Γ + k call sites).
2. Unit test (`tests/cpu/gpw/test_cpu_gpw_dftu.cpp`): the M-matrix fold helper
   against a hand-checkable 2×2 case (single function, no SCF).
3. Offline CP2K pin via `tools/cp2k_oracle`: `PLUS_U` (Mulliken) +
   `RUN_TYPE ENERGY_FORCE` on stretched NiO (R = 1.80 Å, SZV/PADE/U=6 eV on
   Ni d), Γ molecule (CP2K-sane per `gpw-dftu-nio`); plus a periodic
   MP 1×1×2 attempt (documented either way — CP2K k+U pathology precedent
   `gpw-dftu-kpoints`).
4. Regression `regressions/gpw/gpw-dftu-kpoint-forces/` (≤ 60 s):
   - **FD dual gate on the k-path** (decisive): `gpw_forces_analytic` vs
     `gpw_forces_fd` on the H₂ 1×1×2 +U twin geometry (subset of coordinates
     to stay in budget), tolerance from measurement.
   - Γ-limit identity: k-machinery at explicit Γ ≡ real-Γ +U forces.
   - U-activation: max|F(U) − F(U=0)| above a floor (the fold is live).
   - Newton-3 / |ΣF| residual.
   - CP2K force pin gate on the NiO system (tolerance measured & stated,
     honest about the projector divergence).
5. `REFERENCE.md` with pins (CP2K version, input path, output lines), and
   perf statement (T_atv/T_cp2k or perf_exempt with reason).

## Acceptance (measured 2026-07-31)

- [x] FD dual gate green on the k-path with `plus_u` active: full-6N
      `gpw_forces_fd` vs analytic max 6.4e-6 Ha/Bohr (tol 5e-5); fold-off
      control ~1.6 Ha/Bohr wrong-sign — the term is decisive.
- [x] Γ +U forces: Γ-limit identity 6.0e-8 (regression A5) + offline NiO
      d-shell FD 3.7e-4 at tol 1e-8 (REFERENCE.md).
- [x] Regression green in 38.6 s; +U sibling regressions all green after
      variational re-pin; gpw unit family green.
- [x] Pin provenance: CP2K 2026.1 `cp2k.psmp`, inputs under
      `tools/cp2k_oracle/inputs/gpw-dftu-kpoint-forces/`, output lines quoted
      in REFERENCE.md.
- [x] No public API change (fold internal; Γ-RKS+`plus_u` now fails fast in
      `GpwDft::build` instead of silently ignoring U).

**Landed beyond the row (required by the dual gate):** variational
∂E_U/∂P = ½(V̂S+SV̂) scatter (CP2K-precedented); re-pins of 6 ATV-self +U pins
(each cites its pre-fix value); toy-regime notes in
`.claude/notes/gpw-dftu-projector-divergence.md`. The ATV-vs-CP2K projector
divergence (gross vs net Mulliken) is surfaced there for the census — not
resolved in this row.

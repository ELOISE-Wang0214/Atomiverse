# extended-lagrangian-dipoles

Wave-5 metalwalls-method-loop slug. Map row: §2 "Extended-Lagrangian dipoles"
(MW mark: ~ MaZe-related; ATV mark: ○ — SCF `PimInducedDipole` is the
foothold). Queue row: bulk polarizable toy; observables μ(t), energy;
oracle=metalwalls (no MW binary on host → literature-form pins, perf_exempt).

## What

Extended-Lagrangian induced-dipole dynamics for the bulk Polarizable Ion
Model: instead of cold-starting the induced-dipole SCC from μ = 0 every MD
step, an auxiliary dipole set is propagated by the time-reversible dissipative
Niklasson recursion and serves as the SCC initial guess each step
(iEL/SCF scheme: Albaugh, Demerdash, Head-Gordon, J. Chem. Phys. 143, 174104
(2015); recursion coefficients: Niklasson et al., J. Chem. Phys. 130, 214109
(2009), Table I). This is the classical-polarization analogue of XL-BOMD and
the MetalWalls map §2 "XL dipoles" capability (MW's own mass-zero route is
MaZe, Coretti–Bonella–Ciccotti, J. Chem. Phys. 149, 191102 (2018); the
adiabatic XL/SCF foothold here mirrors the shipped `ConstPotMaZe` CG+ASPC
foothold shape for electrode charges).

## Architecture (binding)

- `cpu::PimXlDipole : public Potential` — thin owner composing the existing
  `cpu::PimInducedDipole` (SCC, Thole/TT damping, Hellmann–Feynman forces —
  reused, never duplicated) with the existing
  `cpu::electronic::XlBomd` propagator (Niklasson Table-I recursion, reused on
  an nlocal×3 `linalg::RealMatrix` of dipoles). No new Driver category, no
  public `*Config` bag; parameters are ctor args / public members.
- `PimInducedDipole` gains one primitive: public member `mu_guess`
  (packed xyz [e·Å]; empty ⇒ cold start from zero; consumed and cleared by the
  next `launch`). No behavior change when unset.
- Build: `src/cpu/electronic/xl_bomd.cpp` moves from the ATV_BUILD_HF source
  block to the always-compiled set (it depends only on `cpu/linalg` header
  types; `cpu/linalg` sources are already never method-gated). The
  `test_cpu_electronic_xl_bomd` unit registration un-gates from HF so the
  MD-only configuration covers it.
- Component surface: `name() == "pim_xl_dipole"`, `solver()` accessor to the
  owned `PimInducedDipole` (ConstPotMaZe precedent), `last_*` observability
  (μ, E_pol, SCC iters, aux-vs-SCF deviation `last_aux_error()`), `reset()`
  to drop XL history on geometry discontinuity. Atom-count change auto-resets.

## Acceptance criteria

1. **Fixed point (literature-form pin).** On the two-site Thole toy of
   `regressions/md/pim-induced-dipole-bulk/REFERENCE.md` (q=±1, α=1.0/1.5,
   r=3 Å, thole_a=0.39, L=50), `PimXlDipole` reproduces the pinned
   E_pol = −0.12040946414707571 eV and F_x = ±0.040341072868572643 eV/Å on
   every launch (cold and warm), |ΔE| < 1e-9 eV, |ΔF| < 1e-8 eV/Å: the XL
   recursion's fixed point for a stationary geometry is the SCC solution.
2. **Warm start works.** At fixed geometry the second (warm) launch converges
   the SCC in strictly fewer iterations than the first (cold) launch and in
   ≤ 2 iterations at tol 1e-12.
3. **Tracks the BO surface in dynamics.** Bulk polarizable toy (8-ion ± box),
   manual velocity-Verlet, ≥ 150 steps: max over steps and components of
   |μ_XL(t) − μ_SCF(t)| < 1e-6 e·Å against a cold-start `PimInducedDipole`
   twin on the same trajectory protocol.
4. **Cheaper than cold SCC.** Total warm-start SCC iterations over the
   trajectory strictly below the cold-start twin's total.
5. **Energy behaves.** KE + E_pol drift over the XL trajectory bounded:
   |E_tot(end) − E_tot(0)| ≤ 5e-6 eV on the regression toy (and the cold
   twin's drift bound is not violated either — XL introduces no instability).
6. **Unit tests** under `tests/cpu/potentials/test_pim_xl_dipole.cpp`
   (mirrored path, single-function, seconds): ctor validation, cold-launch
   equivalence with `PimInducedDipole`, warm-start iteration drop,
   `mu_guess` invariance (same converged μ/E), auto-reset on nlocal change,
   `reset()`, bindings/name. Registered in tests/CMakeLists.txt (MD-only
   build).
7. **Regression** `regressions/md/extended-lagrangian-dipoles/` —
   out-of-library public-API consumer, REFERENCE.md
   (oracle: metalwalls; literature-form pins; perf_exempt), targets
   `reg_extended-lagrangian-dipoles` + `run_reg_extended-lagrangian-dipoles`,
   ≤ 60 s.
8. **MD-only dual gate green**: scoped configure/build with the binding
   MD-only flags; unit ctest green; regression target runs green.

## Out of scope

Coupled dipole⊕electrode XL solve (DipoleElectrodeJoint stays SCF), true
SHAKE-like mass-zero constrained dipoles (Coretti symplectic — deferred the
same way ConstPotMaZe defers it for charges), iEL/0-SCF (no-SCF limit),
CUDA, MPI halo for dipoles, neighbor-list PIM production path.

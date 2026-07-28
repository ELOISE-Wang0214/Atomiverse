# Specs index

Live specs under `.claude/specs/`. Completed work is **deleted**; history is
the git log. Method census: `docs/roadmap/`.

## Live — driver registration + AIMD wrappers (2026-07-27)

Operator design (set_evaluator, no public `add`, CPMD as types). **No BC.**

| Spec | Status | One-line |
| ---- | ------ | -------- |
| [driver-api-01-evaluator-registration](driver-api-01-evaluator-registration/spec.md) | approved | `set_evaluator`; delete public `add` / `add_potential` |
| [driver-api-02-aimd-wrappers](driver-api-02-aimd-wrappers/spec.md) | approved | CarParrinello / XlBomdForce / SecondGenCp; strip method flags |

## Live — esmm campaign: ES × MM connection layer (2026-07-26, 16 specs)

Operator-approved three-wave campaign (plan-mode design review + Q&A 2026-07-26):
unify the force-pipeline composition contract, redesign QM/MM around an engine
slot + `qm_z` mask column, generalize CPMD into a shared `cpu/electronic/`
propagation layer. Execution order (operator, 2026-07-26): **esmm-cpmd wave
first** (`/mol:impl-all esmm-cpmd`), then contract + qmmm waves;
esmm-contract-05 rebases its gpw/scf `=`→`+=` sites over the cpmd wave.
The libcint/C++20 sweep gate was cleared by committing it (2026-07-26).

| # | Spec | Deps | One-line |
| - | ---- | ---- | -------- |
| W1-1 | esmm-contract-01-zeroing | — | DONE 2026-07-27 (closed; see git 66664744) — Zero* unification + k/2 oracle fix |
| W1-2 | esmm-contract-02-validation-walks | — | DONE 2026-07-27 (closed; see git 17f52238 + close commit) — node-aware validation walks |
| W1-3 | esmm-contract-03-runslot-nested | — | DONE 2026-07-27 (closed; see git 3480c6ab/6ea2ea23) — run_slot_ nested dispatch + R2 CCCBDB pin |
| W1-4 | esmm-contract-04-forces-requested | W1-2 | DONE 2026-07-27 (closed; see git fc204c63 + R3 commit) — energy-only phases in the Plan + ASPC-under-MC guard + first QM-in-MC regression |
| W1-5 | esmm-contract-05-accumulate-switch | W1-1 | DONE 2026-07-27 (closed) — one accumulation contract: 5 ES evaluators `+=`, kEnergyElectronic, ASPC-fallback 2x guard, first ES+MM co-bound regression |
| W2-1 | esmm-qmmm-01-embed-forces | — | analytic embedding gradient + MM-site forces (third law) |
| W2-2 | esmm-qmmm-02-qm-engine-protocol | W2-1 | QmEngine/EmbeddingEngine protocols; MFM×7 + AM1 + DFTB2 |
| W2-3 | esmm-qmmm-03-electrostatic-rewrite | W2-2, W1-4, W1-5 | qm_z mask column + engine slot + real MM coupling |
| W2-4 | esmm-qmmm-04-mech-pol-port | W2-3 | Mechanical/Polarizable onto QmSubsystem; quantitative re-pins |
| W2-5 | esmm-qmmm-05-md-closure | W2-4 | QM/MM-MD conservation regression; roadmap ○→✓; tutorials |
| W3-1 | esmm-cpmd-01-electronic-home | — | DONE 2026-07-26 (closed; see git c40acdb3) — cpu/electronic/ shared propagators + linalg orbital_tangent extraction |
| W3-2 | esmm-cpmd-02-timestep-protocol | W3-1 | DONE 2026-07-27 (closed; see git) — TimestepProvider/Consumer; cpmd_dt deleted; cpmd×stress gate |
| W3-3 | esmm-cpmd-03-observability | W3-2 | DONE 2026-07-27 (closed; see git) — kEnergyElectronKinetic/Thermostat + CpmdConserved → Thermo drift |
| W3-4 | esmm-cpmd-04-ao-cpmd | W3-1..3 | DONE 2026-07-27 (closed; see git 2de56c65) — Gaussian AO-CPMD on MeanFieldMethod |
| W3-5 | esmm-cpmd-05-scf-xlbomd | W3-4 | DONE 2026-07-27 (closed; see git) — XL-BOMD on the SCF methods (guess ladder; ASPC-vs-explicit-guess bug fixed) |
| W3-6 | esmm-cpmd-06-single-eval | W3-3 | DONE 2026-07-27 (closed; see git) — GpwDft::overlap(); CPMD single KS eval per step |

## Ticket — design review owed (operator TODO, 2026-07-27)

- **Minimize/Optimizer taxonomy is non-orthogonal** (operator): "CG 到底
  optimizer 还是 minimizer，还是某种方法？" Three near-synonyms span the
  layers — namespace `cpu::minimizer::`, Driver preset `Minimize` (the
  loop), component kind `Optimizer` (the step-proposal role), concretes
  `CG`/`LBFGS`/`FIRE` — and "method" is already taken by electronic
  structure. Live symptom: `add(CG{})` vs `set_optimizer(CG{})` confusion
  in tests during the concurrent no-compat API migration. Route:
  /mol:discuss for a vocabulary decision (rename namespace? fold kind
  naming into the driver's? document the tri-layer as-is?); overlaps the
  known "namespaces not yet aligned to directories" debt. No mid-wave
  redesign.

## Ticket — /mol:fix owed (named 2026-07-27, esmm-cpmd-04 close-out)

- **BOMD+d3bj/d4 forces are non-conservative**: `e_total_ += e_disp_`
  (mean_field_method.cpp DFT-D block) has zero dispersion force terms in
  `ScfGradient` — forces are not the gradient of the reported energy for any
  d3bj/d4 run. CPMD refuses d3bj/d4 to keep E_cons honest; BOMD does not.
  Fix = analytic D3(BJ)/D4 gradient legs in scf (GPW twin has them).
- Minor (same sweep): `reg_cpmd-water` single-threaded budget — now 54 s
  after esmm-cpmd-06's single-eval (was 62 s, under budget again);
  `test_cpu_scf_uhf_hessian` + `test_cpu_gpw_dftu` each showed one
  full-suite-context transient (3/3 green isolated) — flaky-test
  investigation. `reg_cpmd-uks` is nondeterministic under threading
  (pre-existing OpenMP reduction-order noise in the Γ-UKS grid/XC path,
  ~3e-5 relative; byte-compares must run OMP_NUM_THREADS=1) — /mol:fix if
  determinism is wanted.
- **Regressions run thread-unpinned; small-system cases pay catastrophic
  fork/join overhead** (named 2026-07-27, esmm-contract-04 R3): the
  `regressions` target runs binaries with no `OMP_NUM_THREADS`, and only ~3
  of ~250 registered cases pin their thread count (`set_threads(1)`;
  precedent `neb_muller_brown`, now also `mc_bias_es`). Measured on
  mc_bias_es (2 atoms, 2 AOs): 7 s CPU at 1 thread, 231 s at 8, unfinished
  after 11 min at the node default 48 — results bit-identical at every
  count, pure OpenMP fork/join waste. Every other small-system case is
  exposed. Fix options: per-case `set_threads` sweep, or the regressions
  CMake target exporting `OMP_NUM_THREADS=1` by default with per-case
  opt-out. Route /mol:fix (mechanical, wide).
- ~~**13 GPW-gated dynamics regressions do not compile at pristine HEAD**~~
  **FIXED 2026-07-28 (e30ea9a3)** — was 14, not 13 (`second-gen-cpmd-kpoints`
  was missing from the ticket). Call-site migration to the post-17f52238 API:
  order 0 is the `set_evaluator_ref` default; `cpmd_mu`/`second_gen_cp_k`/
  `_max_iter` are ctor args; `use_ot` stays a public `GPWEvaluator` member set
  before the move; `last_s2` via the wrapper's `method()`. `cmake --build`
  now exits 0 (703 targets, zero errors); all 14 PASS their drift envelopes.
  No pin or tolerance touched.
- **Two OT cases exceed the 60 s regression budget** (named 2026-07-28, same
  sweep): `second-gen-cpmd-ot` 86 s and `second-gen-cpmd-uks-ot` 95 s at the
  node's default thread count — both PASS; pre-existing case cost that was
  invisible while they did not compile. NOTE the fix is **not** `set_threads(1)`
  (opposite shape to the small-system thread trap below): single-threaded they
  take 110 s and >180 s. Real options: shrink the case (physics/spec change) or
  raise the budget for GPW/OT cases. Route /mol:discuss + /mol:fix.

## Ticket — /mol:fix owed (named 2026-07-28, repo-cleanup wave)

- ~~**DFTB3's SCC potential is not the derivative of its energy**~~
  **FIXED 2026-07-28** (`602c7d02` on the branch, merged as `68df18a8`;
  chain closed by `c0bfffef`..`148a17ed`). Euler identity now 4.4e-16, FD
  5.0e-10. The branch is merged and pruned. Original diagnosis kept below
  for the record.
  Measured on the merge candidate: `V_code` vs numerical `∂(E²+E³)/∂Δq` off by
  **0.457 eV** at Δq=(−0.6,0.3,0.3), while with `U_d=0` (pure second order) the
  residual is 4.5e-10 — the defect is isolated to the third-order block. Euler
  check at the converged Δq: `Σ Δq·V³(code) = 15.496` vs exact `3·E³ = 23.240`
  (a correct derivative satisfies it to 7e-15). The double-counting subtraction
  is also short a term: code subtracts `E²+E³`, a consistent V³ needs `E²+2E³`.
  Symptom: `|E_DFTB3 − E_DFTB2| = 46.96 eV` against a reported `E³ = 6.46 eV`.
  Its regression pin is ALSO stale (−84.5617 vs −89.0775 eV), proven to be
  100% the H-parameter upgrade — but re-pinning was refused: it would bless a
  non-variational solution as dev's reference DFTB3 number. Fix the physics
  first, then re-pin and merge the branch.
- ~~**DFTB2's charge-dependent gradient is wrong for polar systems**~~
  **FIXED 2026-07-28** — root cause was three-fold and none of it was the
  gradient formula alone: (a) the SCC never converged (undamped Picard
  limit cycle; energies were selected by `max_scc` parity), (b) the
  double-counting route used net instead of gross Mulliken population
  (a spurious term linear in dq, ~7.85 eV), (c) the analytic gamma gradient
  carried the wrong SIGN in all three assemblies. The earlier A/B that
  "ruled out the sign flip" was measuring a non-stationary solution, so it
  proved nothing. Now: Anderson/Pulay mixing (`mix_beta` public member,
  budget exhaustion throws), `E = Tr(rho H0) + E2 (+E3) (+E_spin) + E_rep`,
  and analytic-vs-FD forces at 2.5e-07 / 4.9e-09 eV/A (was 0.228 / 0.748).
  Gates tightened 1.0 -> 1e-5 eV/A and 5e-3 -> 1e-8 eV/A^3; both DFTB+
  invariant pins unmoved. Original text kept below.
  (pre-existing, now present in three functions). Every Δq=0 case is exact
  (O2 3.0e-07, H2 5e-03 eV/Å); every polar case is off (H2O mio 0.23/0.75,
  H2O 3ob **0.84/0.24** eV/Å — 8× relative on `F_H1_y`). A converged SCC-DFTB2
  gradient is stationary and must be FD-exact, so the "charge-response
  residual" framing is masking a bug; A/B rules out the cheap fixes (term
  dropped 1.020, sign flipped 2.268 eV/Å). **These gates pass only because
  their tolerance is 1.0 eV/Å** — fix, then tighten. Also: the spin gate
  cannot see it (O2 has Δq ≡ 0 by symmetry); an open-shell Δq ≠ 0 case is owed.
- **D3 dispersion gaps** (from the merged `dftb2-dispersion-d3`): D3 forces are
  ungated in-tree (both D3 gates run `with_forces = false`), and D3 is absent
  from `kStress` on both paths — `d3=true` + `with_stress=true` silently yields
  an electronic+repulsive-only stress. No periodic D3, no periodic stress.
- ~~**The mio-1-1 H–H reader exists twice**~~ **FIXED 2026-07-28**
  (`c0bfffef`): dftb2's private copy deleted, both dftb2 and dftb3 now read
  `internal/mio_hh_table.h`; `d_repulsive_dr` / `d_hopping_ss_dr` /
  `d_overlap_ss_dr` added. Proven a pure refactor — every DFTB2 number
  byte-identical across the rewire.
- **No open-shell Delta-q != 0 force gate** (named 2026-07-28): O2 is
  homonuclear, so `analytic_forces_unrestricted` carries the corrected gamma
  sign with nothing testing it — which is exactly how the sign bug survived
  in three functions. An open-shell polar case is owed.
- **3ob `with_stress` still falls back to multi-SCC strain-FD** although the
  analytic 3ob SK derivatives now exist.
- **MioSketch water charges are physically backwards** (Delta-q = +0.122 on O,
  -0.061 on each H) — a property of the invented exponential O/H sketch rows,
  not of the SCC. Real O/H mio-1-1 SK data is owed, like H-H already has.
- **Seven PW/bands/DOS regressions resolve `tests/data/` relative to CWD**, so
  they FAIL when the `regressions` target runs them from the build directory
  and PASS from the repo root (verified both ways on `reg_bands-dos-si`,
  `reg_dos-spin-si`). Harness defect, not physics. Fix = resolve the data path
  from a compiled-in root or set the target's WORKING_DIRECTORY.
- **Five genuinely failing regressions, attribution proven pre-existing**
  (rebuilt and re-run at `2b350442`, i.e. before the esmm contract wave, with
  byte-identical failure messages): `rigid_water`, `settle_tip3p` (no verdict
  line — SETTLE Frame-lookup budget), `npt_lj` ("<V> off the NPT ensemble
  prediction"), `remd_harmonic` ("accept ratio out of (0.01,0.99): 1.000"),
  `gaussian-electrode-charge` (no verdict line). Plus `python_pure_import`
  (missing `typer` in the env). Classical-MD/MC domain work, unrelated to the
  ES accumulation contract.
- **Four `.claude/specs/` dirs arrived with the merged branches**
  (`dftb2-dispersion-d3`, `dftb2-spin-analytic-forces`, `dftb2-analytic-stress`,
  `dftb2-analytic-forces-threeob`) — house rule says specs die at done. Decide
  their lifecycle (close + delete, or keep as the wave's record).

## Live — deferred PySCF engine gaps (2026-07-24)

Approved for **later** implementation (do **not** auto `/mol:impl` until
scheduled). Unblock `docs/developer/pyscf-method-queue.md` blocked rows.

| Spec | Status | Unblocks | One-line |
| ---- | ------ | -------- | -------- |
| [rsh-eri-cam-b3lyp](rsh-eri-cam-b3lyp/spec.md) | approved, deferred | `cam-b3lyp` | erf-attenuated ERI + CAM-B3LYP RKS energy foothold |
| [periodic-rhf-gamma](periodic-rhf-gamma/spec.md) | approved, deferred | `periodic-rhf-gamma` | Γ-point periodic RHF energy (lattice HF-X) foothold |

Acceptance: `rsh-eri-cam-b3lyp/acceptance.md`, `periodic-rhf-gamma/acceptance.md`.
Independent — no cross-deps. Roadmap: `docs/roadmap/dft.md` §6, `hartree-fock.md` §2.

## Closed — SCF review remediation chain (2026-07-24, all 8 done)

Fixed all 27 findings of the 2026-07-24 /mol:review of `cpu/scf` (verdict
was REQUEST CHANGES: 11 🔴 / 9 🟡 / 7 🟢). Both architecture scrubs restored
to green. Specs deleted on close; history in git (commits f9a6443a..HEAD).

| # | Spec | Deps | One-line |
| - | ---- | ---- | -------- |
| 01 | scf-review-01-linalg | — | DONE 2026-07-24 (closed; see git) — sym_power/eig_nonsymmetric/solve_pcg in cpu/linalg + OrbitalHessian::direct_rpa |
| 02 | scf-review-02-eri-hoists | 01 | DONE 2026-07-24 (closed; see git) — internal/ao2mo.h + internal/spin_orbital_eri.h + mp_denominator |
| 03 | scf-review-03-posthf-dedup | 01, 02 | DONE 2026-07-24 (closed; see git) — eom_ccsd_common (EomImds ctor) + ci_engine; build_imds scrub debt cleared |
| 04 | scf-review-04-meanfield-dedup | — | DONE 2026-07-24 (closed; see git) — FD-Hessian CRTP template + check_hessian_preconditions + GradientResult tag ctor |
| 05 | scf-review-05-api-hygiene | 01 | DONE 2026-07-24 (closed; see git) — grid ctor + throws + units::pi + .cpp moves; BOTH scrubs restored green |
| 06 | scf-review-06-scf-perf | 05 | DONE 2026-07-24 (closed; see git) — XCGridBasis cache + jk_rhf one-sweep + ERI micro-opts (~14% on fixture) |
| 07 | scf-review-07-posthf-perf | 02 | DONE 2026-07-24 (closed; see git) — block eri_dense (10x/54x) + RI-MP2 GEMM + copy elision |
| 08 | scf-review-08-ccsd-spatial-diis | 02, 03, 07 | DONE 2026-07-24 (closed; see git) — CcAmplitudeDiis + spatial CCSD (~87x, oracle parity 1e-9) |

## Live — deferred PySCF only (after harness wave 12)

| Spec | Branch / track | Status |
| ---- | -------------- | ------ |
| [periodic-rhf-gamma](periodic-rhf-gamma/) | pyscf / HF periodic | open (deferred) |
| [rsh-eri-cam-b3lyp](rsh-eri-cam-b3lyp/) | pyscf / RSH | open (deferred) |

Open CP2K specs on `dev` = **0**. Stock VV10 is **blocked**. Incomplete
`[~]` = **0** after close of `gpw-fd-phonon-si` (2026-07-26; reg registered
with 1D-curvature dual-gate + UKS density warm-start). Pending residual
**12** (GapCensus wave 13/14 list); **not** queue exhausted.

## Closed — CP2K harness wave 12 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 12 | xl-bomd-linearized-forces, second-gen-cpmd-kpoints, gpw-d3-atm, gpw-d4-atm, bands-spin-metallic-si, dos-spin-metallic-si | zero-SCF linearized XL; multi-k SGCP; D3/D4 ATM; spin+metallic bands/DOS |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 11 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 11 | ugw-uks-rpbe, gpw-ts-forces, gpw-xdm-forces, gpw-kpoint-uks-stress-gga, cpmd-electronic-thermostat, xl-bomd-uks-kpoints, bands-metallic-si, dos-spin-si, dfpt-finite-q-spin-smearing | UGW@UKS-RPBE; TS/XDM force pins; multi-k UKS stress GGA; CPMD e-T; multi-k UKS XL; metallic bands; spin DOS; finite-q DFPT spin+smear |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 10 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 10 | dfpt-spin-smearing, ot-uks-dftu-smearing, second-gen-cpmd-uks-ot, xl-bomd-kpoints, bands-spin-si, dos-metallic-si, periodic-wannier-uks, berry-spin-si, gpw-mbd-kpoint-forces | DFPT spin+smear; multi-k UKS OT+U+smear; SGCP UKS OT; multi-k XL-BOMD; spin bands; metallic DOS; spin MLWF; spin Berry; multi-k MBD forces |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 9 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 9 | ot-uks-smearing, ot-uks-dftu, ot-dftu-smearing, gpw-d4-forces, gpw-mbd-forces, gpw-uks-multigrid, dfpt-finite-q-smearing, second-gen-cpmd-ot | UKS OT+smear; multi-k UKS OT+U; OT+U+smear; D4/MBD forces; UKS multigrid; finite-q DFPT smear; SGCP OT |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 8 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 8 | periodic-urpa-kpoints, ot-smearing, ot-uks-kpoints, ot-dftu-kpoints, gpw-d4-stress, gpw-vdw-scf-potential, dfpt-metallic-si, dfpt-spin-incommensurate, gpw-optb88-vdw | multi-k URPA; OT+smear; multi-k UKS OT; multi-k OT+DFT+U; D4 analytic stress; nonlocal vdW SCF pot/forces; metallic DFPT; spin-incommensurate DFPT; optB88-vdW |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 7 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 7 | periodic-urpa-si, periodic-ugw-kpoints, gpw-uks-smearing, gpw-hirshfeld-uks, cpmd-uks, second-gen-cpmd-uks, ot-kpoints, ot-dftu, dfpt-spin-si | URPA@UKS; multi-k UGW; GPW UKS FD/MP smear; UKS Hirshfeld; CPMD UKS; SGCP UKS; OT+k; OT+DFT+U; collinear-spin DFPT |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 6 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 6 | gpw-metallic-smearing, gpw-xdm-ice, gpw-uks-rpbe, gpw-dftu-kpoint-stress, periodic-ump2-kpoints, periodic-ugw-si, gpw-hirshfeld-kpoints, dfpt-gga-si, xl-bomd-uks | GPW FD/MP smear; free-atom XDM; Γ-UKS RPBE; multi-k DFT+U stress; multi-k UMP2; UGW@UKS; multi-k Hirshfeld; DFPT-GGA; XL-BOMD UKS |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 5 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 5 | gpw-dftu-kpoints, periodic-gw-finite-q, gpw-uks-revpbe, cpmd-rattle-water, periodic-ump2-si, ot-uks-o2 | multi-k PLUS_U; finite-q G0W0 \(W\); Γ-UKS revPBE; CPMD RATTLE; UMP2@UKS; OT-UKS O₂ |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 4 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 4 | gpw-kpoint-uks-stress, periodic-gw-kpoints, gpw-vdw-df2-ice, gpw-rvv10-ice, gpw-mbd-periodic, pw-spin-gga, gpw-uks-blyp, periodic-rpa-finite-q, gpw-dftu-stress | multi-k UKS stress; multi-k G0W0; vdW-DF2 LMKLL; rVV10; periodic Ewald MBD; PW collinear spin-GGA; Γ-UKS BLYP; finite-q dRPA; DFT+U analytic stress |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 3 (2026-07-25)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 3 | gpw-rpbe-water, gpw-kpoint-uks-forces, gpw-uks-stress, gpw-vdw-df-ice, gpw-mbd-hirshfeld, periodic-rpa-kpoints, periodic-gw-si, dfpt-incommensurate-si | RPBE GGA; multi-k UKS forces; Γ-UKS stress; vdW-DF DRSLL; Hirshfeld MBD/TS; multi-k dRPA; Γ spectral G0W0; incommensurate DFPT |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 2 (2026-07-25)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 2 | gpw-kpoint-forces, gpw-revpbe-si, gpw-kpoints-uks, gpw-uks-forces, gpw-kpoint-stress, gpw-mbd-ice, periodic-mp2-kpoints, periodic-rpa-si, dfpt-finite-q-si, gpw-variable-cell-molopt | MP k-forces; revPBE GGA; k-UKS; Γ-UKS forces; k-stress; free-atom MBD; multi-k RI-MP2; Γ dRPA; finite-q DFPT X; MOLOPT CELL_OPT |

Do not re-create closed slugs without superseding intent.

MetalWalls wave 3 closed 2026-07-24 (specs deleted). **Wave-4 assign**
(2026-07-24) staffed 10 production-depth / residual map rows to
`## in_progress` (still unstarted as of **wave-5 gap-census** 2026-07-24:
MaZe, NeighborList salt/Gaussian, electrode forces, CV-SMD, aqueous EDL,
D-field, charge tools, …). Wave-5 census refilled **pending=10** (per-electrode
\(\eta\), const-charge control, MPI charge-halo, XL dipoles, QEq, PIM-NPT IL,
Voronoi TF, IL EDL box, nanoporous product, TI capacitance). No open
MD/electrochem specs on `dev` — workers write `.claude/specs/<slug>/` on
`wf/<slug>` before impl.

CP2K claude wave 1 closed 2026-07-25 (5 merged). **`xl-bomd-water` closed
2026-07-25** (resume of WIP — verified green, merged, spec deleted).


## Closed — CP2K harness wave 1 (2026-07-25)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 1 | gpw-blyp-water, gpw-uks-pbe, pw-uks-fe-smearing, pw-variable-cell-si, gpw-ts-ice, phonon-dispersion-si, periodic-wannier-si, second-gen-cpmd-water, periodic-ri-mp2-si | GPW BLYP; Γ-UKS PBE O₂; bulk Fe UKS+smear; PW CELL_OPT Si; free-atom TS; phonon ω(q) IFC Fourier; periodic MV MLWF; 2nd-gen CP water; Γ RI-GPW MP2 product |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K parallel wave (2026-07-25)

| Queue | Specs (deleted; see git) | Shipped |
| ----- | ------------------------ | ------- |
| cp2k-method-loop parallel | mpi-pw-scaling-4rank, gpw-analytic-stress, wannier-centers, mpi-gpw-scaling-4rank, cpmd-water | PW/GPW MPI k-RR gates; analytic GPW stress; Foster–Boys Wannier; CPMD H₂O |

## Closed — CP2K xl-bomd-water (2026-07-25)

| Queue | Spec (deleted; see git) | Shipped |
| ----- | ----------------------- | ------- |
| cp2k-method-loop | xl-bomd-water | `XlBomd` + `guess_P` + `GPWEvaluator::xl_bomd`; `reg_xl-bomd-water` NVE drift 1.2e-5 |

## Closed — CP2K capability claude wave 1 (2026-07-25)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop claude wave 1 | gpw-kpoints-si, gpw-variable-cell-ice, dfpt-phonon-si, berry-polarization, periodic-mp2-si | MP k-point GPW RKS; GPW variable-cell relax vs CP2K CELL_OPT; Γ-point DFPT phonons (PW); KSV Berry-phase polarization (PW); Γ-point canonical GPW MP2 foothold |

Do not re-create closed slugs without superseding intent.

## Closed — MetalWalls capability wave 3 (2026-07-24)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| metalwalls-method-loop wave 3 | constant-force-pull-ion, oscillating-efield-md, potential-profile-observer, surface-tension-coupling, virtual-site-multi-topology, dipole-electrode-joint | FixAddForce; FixOscillatingEfield; Ψ(z) PotentialProfile; SurfaceTensionNPT; VirtualSites Weighted3+OutOfPlane3; DipoleElectrodeJoint |

Do not re-create closed slugs without superseding intent.

## Closed — MetalWalls capability wave 2 (2026-07-24)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| metalwalls-method-loop wave 2 | const-pot-cg-large-electrode, electrode-charge-predictor, finite-field-single-electrode, pim-induced-dipole-bulk, piston-electrode-pz, thomas-fermi-electrode, fumi-tosi-nacl, steele-wall-argon, born-repulsion-foothold, capacitance-charge-fluctuation | ConstPotCG; ChargeASPC predictor; Dufils single-electrode finite-field; bulk PIM+Thole; PistonElectrode; TF metallicity; PairBorn/Fumi–Tosi; FixWallLJ1043 Steele; CapacitanceChargeFluctuation |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K capability wave 2 (2026-07-24)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop wave 2 | gpw-stress-water, gpw-bomd-water, pw-analytic-stress-si, fd-phonon-supercell, gpw-dftu-nio, bands-dos-si, pw-uks-kmesh-fe, gpw-pbesol-si, gpw-d4-water | GPW strain-FD stress; GPW BOMD NVE; PW Nielsen–Martin stress; FD supercell phonons; Dudarev DFT+U NiO; bands/DOS Si; k-mesh UKS O₂; PBEsol GPW; periodic D4 |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K capability wave 1 (2026-07-23)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop wave 1 | gpw-water-pbe-mt, gpw-forces-water, gpw-uks-o2-gamma, gpw-pbe-d3bj-ice, pw-lda-si, pw-pbe-gga-si, pw-smearing-al, pw-spin-polarized-o2, pw-bomd-h, ot-minimizer-water | MT water PBE pin; GPW forces; Γ-UKS O₂; GPW D3(BJ); PW LDA Si; true PBE GGA PW; Fermi/MP Al; PW spin O₂; PW BOMD H; OT water pin |

## Closed — MetalWalls capability wave 1 (2026-07-23)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| metalwalls-method-loop wave 1 | const-pot-bo-parallel-plate, freeze-electrode-atoms, const-pot-charge-observer, charging-work-energy, const-pot-md-nvt-short, gaussian-electrode-charge, slab-2d-ewald-foothold, external-efield-md, buckingham-pair-lammps, four-site-tip4p-foothold | ConstPotBO + FixFreeze + electrode charge observer + charging work + short const-pot NVT; Gaussian electrode Coulomb; slab Ewald (`slab_factor`); FixEfield; PairBuck; TIP4P `lj/cut/tip4p/cut` |

## Closed — CPU next campaign (2026-07-23, parallel worktrees)

| Branch (deleted) | Specs | Highlights |
| ------ | ----- | ---------- |
| `feat/cpu-next-md` | parrinello-rahman-npt, nph-ensemble | PR full-cell NPT, isotropic NPH |
| `feat/cpu-next-tb-arch` | tb-sk-product, concrete-pw-method, concrete-gpw-method, concrete-nddo-method | Multi-element DFTB2; concrete PW/GPW/AM1-RM1 |
| `feat/cpu-next-msamp` | qmmm-link-atoms, polarizable-qmmm, multi-cv-metad, multi-window-wham | Link atoms, pol QM/MM, ND metaD, WHAM (molrs) |
| `feat/cpu-next-response` | full-grid-alda, analytic-tda-forces, true-nac-vectors, ri-mp2 | Full MO-product ALDA, CIS TDA forces, NAC vectors, RI-MP2 |

## Closed earlier on `dev` (2026-07-23)

| depth-next | grid-alda-tda, qmmm-electrostatic, tda-nac-skeleton, triclinic-cell-relax | ALDA diagonal, ES QM/MM, NAC strength, triclinic cell |
| open scales | alda-tddft-tda … ri-mp2-skeleton (8) | GFN0, aniso cell, mech QM/MM, metaD, multi-DPD, TDA forces, MP2 |
| method ladder | rks-tda … hamiltonian-remd-lambda (8) | RKS-TDA, DFTB2, PW-cell, DihedralCV, FENE+WCA, DPD, Yukawa, H-REMD |
| breadth-next | soft-core-coulomb … brownian-dynamics (9) | soft Coul, AngleCV, TI, joint cell, string, IRC, Rg², SETTLE, BD |
| breadth-first | FIRE … dimer (6) | soft LJ, Widom, PIMD, isotropic cell, dimer |
| architecture | backend-mirror, cuda-resurrect, cuda-pme-accuracy | layout rewrite, CUDA live, PME Madelung |

## Closed — CUDA OpenMM MD engine (2026-07-23, from origin/dev)

| cuda-openmm-md 01–05 | contract, rebuild, tile-force, multigpu, delete | CS=32 clusters, real masks, tile LJ/Coul, multi-GPU classify; Verlet default removed on CUDA |

# Reference — FD vibrational spectrum from GPW analytic forces (Si₂)

**Status 2026-07-26: pins recorded, gate NOT registered.** Not because the
forces are wrong — see *History*, which retracts the earlier claim. What blocks
it is SCF conditioning plus wall cost for an open-shell dimer in a cluster box;
the numbers are in *Convergence and cost matrix*.

## Method

Finite-displacement (central difference, \(h = 10^{-3}\) Å) Cartesian force
constants, mass-weighted dynamical matrix, full \(3N\) spectrum in cm⁻¹, with
**GPW analytic forces** as the force engine:

```text
cpu::vib::Phonon(1,1,1)
  → vib::Phase::force += cpu::gpw::GPWEvaluator{…, compute_forces = true}
```

No `NeighborList` (GPW is not a `NeighborConsumer`). No public `*Config` /
`make_*`.

An **isolated diatomic** has \(3N = 6\) modes: one stretch plus five near-zero
(three translations, two rotations). That is the correct expectation for this
system — the first attempt compared it against a periodic-crystal pin expecting
three optical modes, which is a different physical system.

## System (the pinned twin)

| Item | Value |
|------|-------|
| System | Si₂, **UKS triplet** (\(^3\Sigma_g^-\), CP2K `LSD` + `MULTIPLICITY 3`; Atomiverse `n_alpha = 5`, `n_beta = 3`) |
| Bond | \(r = 2.246\) Å |
| Cell | 10 × 10 × 10 Å, **`PERIODIC NONE`** + `PSOLVER MT` |
| Placement | **centred** (5, 5, 5) ± 1.123 Å along \(x\) |
| Basis | `DZVP-MOLOPT-GTH` |
| PP | `GTH-PBE-q4` → 8 valence e⁻ |
| XC | PADE |
| Cutoff | 400 Ry (`REL_CUTOFF 40`, `NGRIDS 4`) |
| SCF | `EPS_SCF 1e-7`; CP2K `BROYDEN_MIXING` α = 0.3 |
| FD step | \(10^{-3}\) Å (CP2K `DX 0.001`) |

Centring is mandatory: `GpwDft::build` rejects off-centre `nonperiodic` systems
because Martyna–Tuckerman boundaries are not translation invariant.

## Offline CP2K pins (2026-07-26)

| Field | Value |
| ----- | ----- |
| Engine | **CP2K 2026.1** (`cp2k.psmp`, EasyBuild `2026.1-gmpflf-2026.03`) |
| Host | Arrhenius epyc login, offline only (`mpirun -np 4`, `OMP_NUM_THREADS=2`) |
| Inputs | `tools/cp2k_oracle/inputs/gpw-fd-phonon-si/si2_ef.inp`, `si2_vib.inp` |
| Walls | ENERGY_FORCE **14.3 s**; VIBRATIONAL_ANALYSIS **172.8 s** |

| Quantity | CP2K |
| -------- | ---- |
| Total energy | **−7.613235635300093 Ha** |
| \(F_x\) (atom 0) | **+1.99658355 × 10⁻² Ha/Bohr** (atom 1: −same) |
| \(\lvert\sum F\rvert\) | 9.19 × 10⁻¹⁴ Ha/Bohr |
| **Stretch frequency** | **463.489404 cm⁻¹** |
| Residual trans/rot | ≤ 4.46 cm⁻¹ |

A second twin at 9 Å / 200 Ry also converges in CP2K
(E = −7.609816439926913 Ha, \(F_x = 2.27528784\times10^{-2}\) Ha/Bohr) and is
kept as the cheaper option for when Atomiverse can converge there.

## Convergence and cost matrix (why the gate is not registered)

Measured this wave on the same geometry (\(r = 2.246\) Å, centred, UKS triplet):

| Setting | Atomiverse SCF | CP2K SCF | Atomiverse cost / SCF |
| ------- | -------------- | -------- | --------------------- |
| 8 Å / 150 Ry | **converges** (153 iters) | **fails** | 3.2 s |
| 9 Å / 150 Ry | — | **fails** | ≈ 4.5 s |
| 9 Å / 200 Ry | **fails** (400 iters) | **converges** | ≈ 7 s |
| 10 Å / 400 Ry | **fails** (300 iters) | **converges** | ≈ 38 s |
| 12 Å / 400 Ry (\(r = 2.351\)) | **converges** (45 iters) | — | ≈ 38 s |

Two independent blockers, neither about force correctness:

1. **No overlapping converged window.** Atomiverse converges where CP2K does
   not and vice versa. CP2K uses Broyden mixing; Atomiverse offers linear + DIIS
   + OT, and OT converges onto the *closed-shell* saddle rather than the triplet
   minimum. An open-shell cluster mixer is the actionable next step.
2. **Wall budget.** 12 displacements × 7–38 s = 84–456 s at every setting where
   a CP2K twin exists, versus the ≤ 60 s e2e tier. This is the known GPW
   performance gap, not a defect.

At 12 Å / 400 Ry (where Atomiverse *does* converge as a triplet) the force
quality is fully verified: analytic vs multi-SCF FD **2.5 × 10⁻⁶ Ha/Bohr**,
\(\lvert\sum F\rvert = 1.0\times10^{-7}\).

## History — the earlier block reason was wrong

The wave-12/13 tip recorded: *"GPW analytic forces on Si disagree with energy FD
by ~271× and with CP2K by ~40×, while total energy matches to ~0.014 Ha."* That
was **not reproducible**. Measured on `dev`:

| Case | analytic vs FD |
| ---- | -------------- |
| Si₂, 5.43 Å box, 150 Ry, grid local PP | **1.7e−6** Ha/Bohr |
| Si₂, same, analytic local PP | 3.6e−6 |
| Si₂, same, 400 Ry | 2.1e−6 |
| Si₂ UKS triplet, 12 Å, 400 Ry | **2.5e−6** |

Three **setup** defects produced the symptom:

1. **Forced closed shell.** Si₂'s ground state is a triplet; the case ran RKS.
   With OT the SCF reports `converged` on a saddle of the closed-shell
   functional, where the frozen-density force is not the energy gradient:
   analytic **0.1308** Ha/Bohr with a **7.5e−2** FD residual — exactly the 0.131
   the worker reported. The 0.014 Ha energy offset was the same saddle.
2. **Cluster boundary placement** — a real product defect, now fixed. An atom at
   (0,0,0) with `nonperiodic` boundaries clips the collocated density:
   \(\lvert\sum F\rvert = 1.98\) Ha/Bohr while the energy still looks plausible.
   `GpwDft::build` now throws; `gpw-d3-forces` measured the same trap at
   **31.9 Ha/Bohr** on a water dimer and asserts \(\lvert\sum F\rvert\).
3. **Twin mismatch.** An isolated-cluster spectrum was compared against a
   `FULLY_PERIODIC T` crystal pin expecting three optical modes.

The analytic force code was never at fault. Do not re-open it as debt.

## Acceptance (when the mixer lands)

| Check | Tol | Rationale |
| ----- | --- | --------- |
| Mode count | exactly 6 | \(3N\) for a dimer |
| Near-zero modes | 5 with \(\lvert\omega\rvert < 40\) cm⁻¹ | 3 translations + 2 rotations; CP2K's own residual is ≤ 4.5 cm⁻¹ |
| Stretch vs CP2K | **25 cm⁻¹** | same-method twin at identical grid/basis/PP; CP2K's own FD step and grid set the floor |
| \(F_x\) vs CP2K | 5e−3 Ha/Bohr | same justification as `reg_gpw-forces-water`'s 1e−3 tier, loosened for the coarser cluster box |
| Wall | ≤ 60 s | e2e tier |

## Notes

- Unit wiring gate (green, registered): `test_cpu_vib_phonon_gpw_wiring` —
  `Phonon` accepts `GPWEvaluator` without a `NeighborList`.
- Sibling FD product: `regressions/pw/fd-phonon-supercell` (SW potential).
- Sibling PW DFPT: `regressions/response/dfpt-phonon-si`.
- Water is the green GPW force foothold: `regressions/gpw/gpw-forces-water`.

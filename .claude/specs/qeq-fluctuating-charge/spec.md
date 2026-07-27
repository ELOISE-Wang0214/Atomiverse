# qeq-fluctuating-charge — QEq / fluctuating-charge equilibration (LAMMPS `fix qeq/point` twin)

- Loop: metalwalls-method-loop (wave 5) · Map: §2 "Fluctuating-charge / QEq" (○ → ●)
- Oracle: **LAMMPS** (`fix qeq/point`, QEQ package) — binding priority rule: feature exists in LAMMPS.
- Status: code-complete (2026-07-27 — all 6 acceptance criteria verified: unit
  gate `test_cpu_md_qeq` passed, regression pins OK at 1e-13, perf ratio 0.45,
  MD-only flags confirmed OFF in CMakeCache)

## What

A charge-equilibration Component `cpu::md::QEq` for **electrolyte / molecular**
fluctuating charges (electronegativity equalization) — explicitly *not* the
const-potential electrode solvers (`ConstPotBO`/`ConstPotCG`/`ConstPotMaZe`
stay untouched; they solve a different problem: Gaussian elastance at fixed
Ψ on electrode groups).

QEq (Rappé–Goddard 1991 style, LAMMPS `fix qeq/point` algorithm) minimizes

E(q) = Σ_i (χ_i q_i + ½ η_i q_i²) + Σ_{i<j, r≤rc} q_i q_j / r_ij

subject to Σ_i q_i = q_tot (default 0), where the Coulomb coupling kernel is
**exactly LAMMPS qeq/point's**: J_ij = 1/r_ij within a hard cutoff, with the
literal LAMMPS convention of **no qqrd2e prefactor** on the off-diagonal
(verified from `src/QEQ/fix_qeq_point.cpp`: `H.val = 0.5/r` on a full list,
`sparse_matvec` adds `eta_i x_i`; RHS `b_s = -chi_i`, `b_t = -1`; constraint
via `q = s − (Σs/Σt) t`). χ in eV, η in eV/e² per the LAMMPS param file.

Stationarity + constraint = the linear system
(η_i δ_ij + J_ij) q = −χ + λ·1, solved with the same two-RHS (s, t) scheme as
LAMMPS, then q = s − ((Σs − q_tot)/Σt)·t. Solve is dense Cholesky (toy/small-N
scope, ≤ few hundred atoms; matches ConstPotBO's dense path).

## Component shape (ARCH LAW)

- `cpu::md::QEq : public cpu::Component`, header
  `include/atomiverse/cpu/dynamics/qeq.h`, impl `src/cpu/dynamics/qeq.cpp`
  (sibling of `const_pot_bo`; registered in root `CMakeLists.txt` MD block).
- Constructor args + public members, **no Config bag, no factories**:
  `QEq(double cutoff)`; per-type params via
  `set_type(int type, double chi, double eta)` (lvalue & rvalue chaining like
  `ConstPotBO::add_electrode`); members `qtot` [e].
- Default binding `md::Phase::pre_force`, order −50 (before pair/kspace, same
  slot as ConstPotBO) — Component on `cpu::md::MD`, never a Driver category,
  never inside KSpace*.
- Writes equilibrated charges to `kAtom_charge` (Frame is sole I/O).
- Observables: `charges()`, `qeq_energy()` [eV, k_e≡1 kernel convention as
  LAMMPS], `net_charge()`, `n_solved()`.
- Distances via `frame.box.min_image` (primary-image kernel; the pinned toy is
  an isolated cluster in a large box, identical to LAMMPS `boundary f f f`).
- Errors throw: no `set_type` for a type present in the frame, cutoff ≤ 0,
  η ≤ 0.

## Oracle pin (offline, never at ctest time)

- Geometry: formaldehyde-like H₂CO 4-atom cluster (small organic toy),
  C=O 1.205 Å, C–H 1.111 Å, ∠HCH 116.13°, centered in a 40 Å f/f/f box.
- Params (RG91-flavored, identical in both codes): C χ=5.343 η=10.126,
  O χ=8.741 η=13.364, H χ=4.528 η=13.8904 (eV, eV/e²).
- LAMMPS: `units metal`, `atom_style charge`, `pair_style coul/cut 10.0`,
  `fix qeq/point 1 10.0 1e-10 500 param.qeq`, `run 0`; dump per-atom q,
  thermo `ecoul`.
- Inputs + pin notes: `tools/metalwalls_oracle/inputs/lammps/qeq-fluctuating-charge/`
  + `tools/metalwalls_oracle/pin_notes_qeq-fluctuating-charge.md`.

## Acceptance criteria

1. `cpu::md::QEq` compiles in the MD-only configure (HF/PW/GPW/CUDA/TORCH all
   OFF) and is exercised by unit tests under `tests/cpu/dynamics/test_qeq.cpp`
   (single-function, ms-scale): two-atom closed-form charge
   q = (χ₂−χ₁)/(η₁+η₂−2/r); Σq = qtot conservation (incl. nonzero qtot);
   beyond-cutoff decoupling; energy vs closed form; throw on missing type
   params; charges written to `kAtom_charge`; default binding at `pre_force`.
2. Regression `regressions/md/qeq-fluctuating-charge/` (out-of-library
   consumer, public headers only, ≤ 60 s, hard-coded pins):
   per-atom q_i match the LAMMPS `fix qeq/point` pins to ≤ 1e-6 e;
   Coulomb energy of the equilibrated system via `cpu::PairCoulCut` matches
   the LAMMPS `ecoul` pin to ≤ 1e-6 eV; `REFERENCE.md` states
   `oracle: lammps` + `oracle_version` + cited numbers/units.
3. Targets `reg_qeq-fluctuating-charge` + `run_reg_qeq-fluctuating-charge`
   exist and pass.
4. Dual gate green: scoped build (atomiverse + test_cpu_md_qeq +
   reg_qeq-fluctuating-charge), `ctest -R test_cpu_md_qeq`, and the regression
   run target.
5. Perf: fair wall clock — N repeated equilibrations of the same toy in ATV vs
   LAMMPS `run N` loop time with `fix qeq/point` every step; ratio recorded in
   REFERENCE.md; within 10× or perf_exempt with reason.
6. Existing electrode components untouched (no duplication of
   ConstPot*/PairGaussianCoul machinery).

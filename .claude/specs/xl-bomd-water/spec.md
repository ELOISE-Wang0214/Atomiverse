# xl-bomd-water — Extended-Lagrangian BOMD (Niklasson dissipative XL-BOMD)

- **Loop:** cp2k-method-loop (wave 3, resume `[~]`, write-first)
- **Domain:** dynamics
- **Dep:** `gpw-bomd-water` ✓ (GPW BOMD NVE foothold, merged wave 2)
- **Branch:** `wf/xl-bomd-water` off `dev`

## Problem

The MD roadmap lists XL-BOMD as open (`ASPC ≠ XL`): Atomiverse has the Kolafa
ASPC predictor (`cpu::scf::ASPCPredictor`, an *extrapolation*) and a plain
warm-start (previous MO coefficients) on the GPW BOMD path, but **no
extended-Lagrangian auxiliary electronic degree of freedom**. XL-BOMD
(Niklasson, PRL **100**, 123004 (2008); dissipative form Niklasson et al.,
JCP **130**, 214109 (2009)) propagates an auxiliary density matrix `P`
alongside the nuclei with a time-reversible Verlet-like recursion plus a
weak dissipative term, and uses it as the SCF initial guess each step:

```
P_{n+1} = 2 P_n − P_{n−1} + κ (D_n − P_n) + α Σ_{k=0..K} c_k P_{n−k}
```

where `D_n` is the SCF density at step n and `(κ, α, c_k)` come from Table I
of JCP 130, 214109 (2009) for `K ∈ [3, 9]`. Unlike ASPC extrapolation, the
auxiliary density is a *dynamical variable* — the propagation is
time-reversible up to the `O(α)` dissipation, which removes the systematic
energy drain/drift that plain guess-reuse (broken time reversibility) causes
in long BOMD trajectories.

## Scope (foothold)

- **In:** dissipative XL-BOMD auxiliary-density propagation on the **GPW RKS
  Γ-point diagonalization/DIIS** path (the `gpw-bomd-water` twin: H₂O,
  DZVP-MOLOPT-GTH, GTH-PBE, PBE, MT nonperiodic Poisson). SCF is still
  converged to `tol` each step (SCF-converged XL-BOMD; the guess quality is
  what XL improves).
- **Out:** UKS / OT / k-point XL; the zero-SCF linearized-functional XL-BOMD
  force (`∂U/∂R` at fixed `P`); thermostatted XL; CPMD (held slug
  `cpmd-water`). Refused loudly, never silently degraded.

## Design (architecture-conform)

1. **`cpu::gpw::XlBomd`** — new concrete class,
   `include/atomiverse/cpu/gpw/xl_bomd.h` + `src/cpu/gpw/xl_bomd.cpp`.
   Constructor `XlBomd(int k = 5)` (throws unless `3 ≤ k ≤ 9`); members are
   the Niklasson Table-I coefficients; methods
   `advance(const linalg::RealMatrix& d_scf)` (first call seeds history,
   later calls propagate), `ready()`, `aux()` (current auxiliary density =
   next SCF guess; throws before seeding), `reset()`, plus read-only
   `order() / kappa() / alpha() / c()` for tests. No `*Config` bag, no
   factory.
2. **`GpwDft::run(atoms, guess_C, guess_P)`** — third optional parameter:
   an RKS initial-density override (the auxiliary density). Non-null
   `guess_P` with UKS or OT configured throws `std::invalid_argument`
   (no-fallback); wrong shape throws. `guess_C` semantics unchanged.
3. **`GPWEvaluator`** members `bool xl_bomd = false; int xl_bomd_k = 5;`.
   `build()` refuses `xl_bomd` with UKS occupations or the OT path.
   `launch()`: guess SCF from `xl.aux()` when ready; after a converged SCF,
   `xl.advance(P_scf)`; on a non-converged SCF, `xl.reset()` (discontinuity
   policy). Force-path SCF (`gpw_forces_analytic`) warm-starts from the
   just-converged `C` in XL mode (same geometry, same fixed point).
4. **`gpw_forces_analytic(atoms, cfg, guess_C = nullptr)`** — internal
   header gains an optional warm-start pass-through (default keeps the cold
   behavior bit-for-bit for existing callers).

Frame stays the sole I/O currency; the auxiliary density is **private method
workspace** (never on Frame, never public state beyond the propagator
object). Public units unchanged (eV / Å / fs); the propagator itself is
unit-agnostic linear algebra on AO density matrices [Hartree-atomic
internals of GPW].

## Oracle / pin

Offline per `tools/cp2k_oracle/README.md`. **No CP2K binary or module on
this host** (re-checked 2026-07-25: `command -v cp2k{,.psmp,.ssmp,.sopt}`
empty, `module avail cp2k` empty), same status as the `gpw-bomd-water` pin.
The pin is therefore the **literature XL-BOMD/NVE conservation envelope**
(Niklasson 2008/2009; Marx & Hutter BOMD chapter; Frenkel & Smit §4.3.4)
plus the exact Table-I coefficient pin at unit level, with the CP2K twin
input recorded for a future re-pin. Force quality for this geometry/method
family remains offline-pinned vs CP2K 2026.1 in `gpw-forces-water`.

## Deliverables

- `include/atomiverse/cpu/gpw/xl_bomd.h`, `src/cpu/gpw/xl_bomd.cpp`
- `GpwDft::run` guess_P + `GPWEvaluator` XL wiring + forces warm-start
- Unit tests `tests/cpu/gpw/test_gpw_xl_bomd.cpp` (single-function, fast)
- Regression `regressions/dynamics/xl-bomd-water/reg_xl-bomd-water.cpp`
  + `REFERENCE.md` (registered `atv_add_regression(dynamics xl-bomd-water)`,
  gated `ATV_BUILD_GPW AND ATV_BUILD_MD`)
- Product page `docs/tutorial/methods/xl-bomd-water.md`

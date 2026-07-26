# berry-spin-si — acceptance

All criteria concrete: observable + tolerance vs pin. Gate = scoped
build/test: `atomiverse` + `test_cpu_pw_berry` + `reg_berry-spin-si`.

## Criteria

1. **Unit — existing berry math still green** (`test_cpu_pw_berry`):
   string geometry, `phase_from_link_dets`, `wrap_centered`, Hermitian Jacobi
   tests unchanged.
2. **Unit — spin path is accepted at the API boundary**:
   constructing `BerryPolarization` + documenting that
   `method.spin_polarized` no longer hard-rejects before SCF (the e2e gate
   is the regression; unit layer may assert wrap quantum helpers if any
   spin-specific free function is added). No full SCF in unit tests.
3. **Regression `regressions/pw/berry-spin-si`** (public API only,
   hard-coded literature pins, ≤ 60 s):
   a. `spin_polarized=true`, `multiplicity=1` on diamond Si; both SCF runs
      converge.
   b. Ionic: ideal \(p^\mathrm{ion}=(1,1,1)\) tol 1e-9; displaced
      \(\Delta p^\mathrm{ion}_3 = 0.08\) tol 1e-9.
   c. **Centrosymmetric pin** (Spaldin 2012): distance of every ideal
      total \(p_i\) to nearest multiple of 0.5 < **0.05**.
   d. **Born-charge pin** (ASR): \(|Z^*|=|\Delta p_3|/\delta < 0.5\) with
      \(\delta=0.02\).
4. **Architecture**: no public `*Config`/factory; closed-shell path
   behaviour unchanged; `run_collinear_scf` is private internal (not on
   public include path for users beyond what bands/Berry already share).
5. **Perf**: fair CP2K wall or `perf_exempt` with reason (host has no
   loadable CP2K/QE PW twin — same precedent as `berry-polarization`).

## Status (2026-07-26)

- [x] 1 unit suite green (`test_cpu_pw_berry`: strings, phase math, wrap,
      nspin=2 quantum-1 wrap, Hermitian Jacobi)
- [x] 2 API boundary: `spin_polarized` no longer hard-rejects; collinear
      insulator path runs to completion
- [x] 3 regression green (`reg_berry-spin-si`): ideal phases quantize to 0,
      p_ion=(1,1,1), Spaldin lattice pin, Z*=0.0000 ASR; wall ~32 s
- [x] 4 architecture: `run_collinear_scf` private internal; no public
      `*Config`/factory; closed-shell path retained
- [x] 5 perf: `perf_exempt=true` in REFERENCE.md (no CP2K/QE PW twin wall)

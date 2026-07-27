# Spec: dftb2-dftbplus-energy-pin

## Goal

Close the “finite energy only” hygiene gap on ATV `DFTB2` by pinning **total
energy** against offline **DFTB+** (family-locked oracle for DFTB) on H₂ with
**mio-1-1** H–H parameters.

## In scope

- Transcribe mio-1-1 `H-H.skf` (SK grid + spline repulsive + H onsite/U) into
  static headers under `basis_data/` (no runtime `.skf` parser).
- DFTB2 H–H path: DFTB+-compatible “new” SK interpolation + spline repulsive.
- Regression `regressions/semiempirical/dftb2-dftbplus-energy-pin/` with
  `REFERENCE.md` stating `oracle: dftb+` and `oracle_version`.
- Unit smoke: H₂ energy finite and near pin (loose unit tol); existing SK tests
  remain green.

## Out of scope

- Full multi-element mio-1-1 / 3ob tables (→ `dftb2-3ob-organic`).
- Analytic forces (→ `dftb2-analytic-forces`); forces stay re-SCC FD.
- DFTB3 / spin / periodic.
- Live DFTB+ inside ctest.

## Oracle pin (offline)

| Field | Value |
| ----- | ----- |
| oracle | dftb+ |
| oracle_version | 24.1 |
| SK set | mio-1-1 (`dftbparams/mio` skfiles) |
| System | H₂, r(H–H) = 0.74 Å, isolated, charge 0, closed shell |
| SCC | Yes, tol 1e-8 |
| Total energy | −0.6749509324 Ha = **−18.36635052 eV** (Ha→eV = 27.211386245988) |
| Wall (DFTB+) | ~0.15 s (cold start; fair ratio optional) |

Inputs archived under `tools/tb_oracle/inputs/dftbplus/h2_mio_scc/`.

## Acceptance

1. `reg_dftb2-dftbplus-energy-pin` asserts
   `|E_ATV − E_pin| ≤ 1e-3 eV` at the stated geometry (public Evaluate + DFTB2).
2. `REFERENCE.md` documents oracle, version, units, Hamiltonian, SK set.
3. SE-only cmake; scoped build targets only.
4. Existing `test_cpu_dftb2` / `test_cpu_dftb2_sk` pass.
5. No public `*Config` / `make_*` / runtime SK parsers.

## Perf

`T_atv / T_oracle ≤ 10` preferred; pin-fidelity slug may set `perf_exempt` if
oracle cold-start dominates (document wall times in REFERENCE).

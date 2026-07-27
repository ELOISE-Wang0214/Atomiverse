# Acceptance — `dftb0-h2`

## Dual gate

1. **Unit** (`test_cpu_dftb0`):
   - H2 energy finite; writes `frame.meta[kEnergy]`
   - Forces finite and non-zero along bond when `with_forces=true`
   - Unsupported Z throws
   - Force component matches central FD of E within 0.05 eV/Å (non-SCC FD is
     self-consistent with the same energy path)

2. **Regression** (`reg_dftb0-h2` / `run_reg_dftb0-h2`):
   - H2 R=0.74 Å total energy matches REFERENCE pin within tol
   - Energy differs at a second separation (pipeline not flat)
   - Optional force pin if REFERENCE states one

## REFERENCE hygiene

`regressions/semiempirical/dftb0-h2/REFERENCE.md` must include:

```text
oracle: dftb+
oracle_version: <string>
```

and document Hamiltonian (non-SCC), geometry, SK set, units.

## Perf

Prefer \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\). If offline DFTB+ binary
is unavailable on the host, set `perf_exempt` with reason (no fair wall).

## Build surface

SE-only cmake only. Scoped targets:

```text
atomiverse test_cpu_dftb0 reg_dftb0-h2
```

Never bare `cmake --build build`. Never enable HF/MD/CUDA/MINIMIZER.

## Done when

- Spec green, dual gate pass, commit on `wf/dftb0-h2`, `acceptance_met=true`

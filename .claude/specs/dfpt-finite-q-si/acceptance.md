# Acceptance — `dfpt-finite-q-si`

## Dual gate

| Gate | Criterion |
| ---- | --------- |
| Correctness | All regression checks in `reg_dfpt-finite-q-si` PASS; unit `test_cpu_pw_dfpt` PASS |
| Perf | `perf_exempt=true` (no CP2K DFPT wall on pin host) **or** ratio ≤ 10 |

## Regression checks (hard-coded pins in `REFERENCE.md`)

1. SCF + all response SCFs converged.
2. Six finite frequencies at \(X=(1/2,0,0)\); all finite / non-NaN.
3. Acoustic branch at X: TA-like pair and LA inside literature-biased windows.
4. Optical branch at X: TO/LO inside windows (minimal ecut bias vs Giannozzi LDA).
5. Hermiticity diagnostic \(\max|C-C^\dagger|/2\) ≤ pin.
6. Offline frozen-phonon energy second-difference pin for one Cartesian
   pattern (or diagonal Re \(C\)) within tol.
7. `q_frac={0,0,0}` still matches Γ C00 pin (smoke, or shared unit path).

## Unit

- `ewald_dynamical_matrix_q` at \(q=0\) matches `ewald_dynamical_matrix` (≤ 1e-10).
- Finite-q Ewald Hermitian to 1e-12; force-FD cross-check optional.
- Public guards: empty UPF / non-LSDA / smearing / spin still throw.
- `q_frac` default is Γ.

## Done when

- `acceptance_met=true`, branch `wf/dfpt-finite-q-si` committed, scoped build green.
- Spec remains until harness merges and deletes it.

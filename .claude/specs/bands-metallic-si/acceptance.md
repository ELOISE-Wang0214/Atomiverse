# Acceptance — bands-metallic-si

## Correctness

1. `Bands::compute` with `smearing = fermi_dirac` on diamond Si **does not
   throw**.
2. SCF converges (`br.scf_converged`).
3. Path has ≥ 2 k-points; each k has `n_bands` eigenvalues sorted ascending.
4. Path includes \(\Gamma\) and \(X\) (diamond path).
5. Integer-fill spectrum still **open** (`gap_ry > 0`) on this twin.
6. `e_fermi_ry` finite, \(|E_F| < 10^3\) Ry, and within
   \(\pm 4\,\mathrm{degauss}\) of \([E_\mathrm{vbm}, E_\mathrm{cbm}]\).
7. Local same-\(k\) gaps open for valence/conduction pair
   (\(n_\mathrm{occ}=4\) for Si).
8. Fixed-occupation `Bands` path remains unchanged (no reg of
   `bands-dos-si` / unit helpers broken).
9. Spin + smear still throws a clear message (not implemented this slug).

## Performance

- Prefer \(T_\mathrm{atv}/T_\mathrm{cp2k} \le 10\) when a fair twin wall
  exists.
- Else **`perf_exempt=true`** with reason (no CP2K module on host) —
  same policy as `dos-metallic-si` / `bands-dos-si`.

## Dual gate (worker)

```bash
cmake --build build -j$(nproc) --target atomiverse test_cpu_pw_bands_dos reg_bands-metallic-si
ctest --test-dir build -R '^test_cpu_pw_bands_dos$' --output-on-failure
cmake --build build -j$(nproc) --target run_reg_bands-metallic-si
```

## Done when

- Spec acceptance criteria pass on `wf/bands-metallic-si`.
- Regression + unit test green.
- Commit on branch; builds wiped after.

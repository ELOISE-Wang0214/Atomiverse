# Acceptance — `gpw-hirshfeld-uks`

## Dual gate

1. **Unit** (`test_cpu_gpw_hirshfeld`):
   - Γ UKS + `mbd` + `hirshfeld` does **not** throw; reports
     `volume_ratios.size()==n_atoms`, finite attractive `e_disp`, and
     \(v_A\in(0.2,2.0)\).
   - Closed-shell water UKS (\(N_\alpha=N_\beta\)) + hirshfeld-MBD E_disp
     matches RKS Hirshfeld-MBD within a loose numerical gate (same total \(\rho\)).
   - Free-atom multi-k RKS path unchanged.

2. **Regression** (`reg_gpw-hirshfeld-uks`):
   - O₂ triplet Γ-UKS GPW PBE, `mbd` + `hirshfeld`, MT isolated.
   - Hard gate on **E_disp** to **1e-8 Ha** vs offline pin in `REFERENCE.md`.
   - Volume ratios reported (2 O atoms) and within physical bounds.

## Perf

`perf_exempt = true` — primary deliverable is UKS volume-ratio + E_disp
parity (same policy as `gpw-mbd-hirshfeld` / `gpw-hirshfeld-kpoints`).

## Done when

- Spec acceptance met; unit + `run_reg_gpw-hirshfeld-uks` green.
- Commit on `wf/gpw-hirshfeld-uks`.

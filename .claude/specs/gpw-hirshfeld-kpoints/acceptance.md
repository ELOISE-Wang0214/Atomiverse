# Acceptance — `gpw-hirshfeld-kpoints`

## Dual gate

1. **Unit** (`test_cpu_gpw_hirshfeld`):
   - Multi-k RKS + `mbd` + `hirshfeld` does **not** throw; reports
     `volume_ratios.size()==n_atoms`, finite attractive `e_disp`, and
     \(v_A\in(0.2,2.0)\).
   - UKS + `hirshfeld` still throws (RKS-only gate retained).
   - Free-atom multi-k MBD (`hirshfeld=false`) path unchanged (no throw).

2. **Regression** (`reg_gpw-hirshfeld-kpoints`):
   - Water-dimer GPW PBE, **k-mesh 2×2×2**, `mbd` + `hirshfeld`.
   - Hard gate on **E_disp** to **1e-8 Ha** vs offline pin in `REFERENCE.md`.
   - Volume ratios reported and within physical bounds.

## Perf

`perf_exempt = true` — primary deliverable is multi-k volume-ratio + E_disp
parity (same policy as `gpw-mbd-hirshfeld` / `gpw-mbd-ice`).

## Done when

- Spec acceptance met; unit + `run_reg_gpw-hirshfeld-kpoints` green.
- Commit on `wf/gpw-hirshfeld-kpoints`.

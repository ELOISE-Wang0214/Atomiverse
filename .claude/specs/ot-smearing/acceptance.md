# Acceptance — ot-smearing

| Gate | Status | Evidence |
| ---- | ------ | -------- |
| A1 build accepts RKS OT+FD | **met** | `test_cpu_gpw_smearing`, reg A1 |
| A2 OT+FD converges | **met** | reg 12 OT iters |
| A3 \|F_OT−F_DIIS\| ≤ 1e-4 Ha | **met** | 1.15×10⁻⁷ Ha |
| A4 Ne = 18 | **met** | reg A4 |
| A5 frozen pin ≤ 1e-3 Ha | **met** | pin −92.1263304935 Ha |
| A6 multi-k UKS OT throws | **met** | reg A6 |
| Unit tests | **met** | `test_cpu_gpw_smearing`, `test_cpu_gpw_ot` |
| Perf ≤10× | **met** | r ≈ 0.20 (ATV 9–18 s / CP2K ~89.5 s) |

`acceptance_met = true`

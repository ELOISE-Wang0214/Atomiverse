# Acceptance — `ot-uks-o2`

| # | Criterion | Verify |
|---|-----------|--------|
| A1 | `use_ot=true` + UKS no longer throws; OT SCF converges on O₂ | regression |
| A2 | Closed-shell OT / DIIS UKS paths unchanged (unit suite subset) | `ctest -R test_cpu_gpw_ot\|test_cpu_gpw_uks` |
| A3 | \|E_ATV,OT − E_CP2K,OT\| ≤ **2.5e-3 Ha** (250 Ry grid twin class) | regression |
| A4 | \|E_ATV,OT − E_ATV,DIIS\| ≤ **1e-4 Ha** | regression |
| A5 | ⟨S²⟩ within **0.15** of CP2K pin (contamination class) | regression |
| A6 | No public `*Config` / `make_*` / `build_*` for this feature | review / scrub |
| A7 | Offline pin: CP2K version, input path, geometry in `REFERENCE.md` | REFERENCE.md |
| A8 | Perf: ratio ≤ 10 **or** `perf_exempt` with reason | REFERENCE.md |
| A9 | Frame energy eV consistent with Ha × conversion | regression |

**`acceptance_met=true` only when A1–A9 hold.**

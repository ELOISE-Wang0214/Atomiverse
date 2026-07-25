# Acceptance — `gpw-uks-forces`

## Correctness

| ID | Check | Gate |
|----|--------|------|
| A1 | Offline CP2K 2026.1 `ENERGY_FORCE` pin for O₂ triplet UKS-PBE / DZVP-MOLOPT-GTH / GTH-PBE-q6 / MT / 12 Å / R=1.2075 Å (same twin family as `gpw-uks-pbe`) | `REFERENCE.md` cites engine, input path, numbers |
| A2 | Regression `reg_gpw-uks-forces` max component \(\|F_\mathrm{ATV}-F_\mathrm{CP2K}\| \le 5\times 10^{-3}\) Ha/Bohr (first-green multigrid twin; unit analytic-vs-FD 5e-5) | dual gate |
| A3 | Energy soft check \(\|E-E_\mathrm{CP2K}\| \le 1\times 10^{-3}\) Ha | same as energy twin residual |
| A4 | Unit: analytic UKS forces vs FD energy gradient on a tiny open-shell system (max \(\|\Delta F\|\) small) | `test_cpu_gpw_forces` |
| A5 | Closed-shell RKS force path bit-stable / still green under existing unit + `gpw-forces-water` not required in worker gate | no RKS regression |
| A6 | Frame public forces eV/Å; pin compare converts with `units::kEvPerAngstromInHartreePerBohr` | code review |
| A7 | No public Config bags / factories | architecture scrub |

## Perf

| ID | Check |
|----|--------|
| P1 | Fair wall \(T_\mathrm{atv}/T_\mathrm{cp2k} \le 10\), or `perf_exempt` with reason |

## Done when

`acceptance_met=true` and `status=green` with commit on `wf/gpw-uks-forces`.

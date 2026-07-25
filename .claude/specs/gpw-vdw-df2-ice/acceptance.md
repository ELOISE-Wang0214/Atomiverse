# Acceptance — gpw-vdw-df2-ice

| # | Criterion | Evidence |
| - | --------- | -------- |
| 1 | `GpwDft::vdw_df2` evaluates post-SCF LMKLL \(E_c^{\mathrm{nl}}\) into `e_disp` | `reg_gpw-vdw-df2-ice` |
| 2 | Kernel matches libvdwxc DF2 on fixed synthetic density to ≤ 1e−6 Ha | `test_cpu_gpw_vdw_df` |
| 3 | Water-dimer regression hard pin \|E − E_ref\| ≤ 5e−5 Ha | REFERENCE.md + reg |
| 4 | Mutual exclusion vs d3/d3bj/d4/ts/mbd/vdw_df; throws on nonperiodic/multigrid | unit / code path |
| 5 | No public `*Config` / `make_*` / live oracle in ctest | review |
| 6 | perf_exempt (no fair CP2K wall on pin host) | REFERENCE.md |

`acceptance_met=true` when 1–5 green on `wf/gpw-vdw-df2-ice`.

# Spec — `gpw-vv10-ice`

**Slug:** `gpw-vv10-ice`  
**Method:** Stock VV10 nonlocal correlation (Vydrov & Van Voorhis 2010) ≠ rVV10  
**Domain:** GPW (periodic)  
**Branch:** `wf/gpw-vv10-ice`  
**Status:** **blocked** (reconfirmed 2026-07-26, CP2K capability loop wave 8)

## Intent

Ship post-SCF stock **VV10** nonlocal correlation on Γ-point GPW as a product
flag distinct from the already-productized **rVV10** path
(`GpwDft::rvv10` / `reg_gpw-rvv10-ice`). Target geometry: ice-relevant water
dimer (same cell as `gpw-rvv10-ice` / `gpw-vdw-df-ice`).

## Oracle contract (this loop)

This loop pins against **offline CP2K** (`tools/cp2k_oracle/`). Queue text
claimed:

```text
CP2K `VDW_POTENTIAL NON_LOCAL TYPE VV10`
```

## Gap finding (binding block)

**CP2K 2026.1 does not expose stock VV10 as a NON_LOCAL type.**

Verified on Arrhenius tree install (`work/cp2k/install`, libcp2k 2026.1):

| Probe | Result |
| ----- | ------ |
| Input `TYPE VV10` + `rVV10_kernel_table.dat` + `PARAMETERS 5.9 0.0093` | **ABORT** `invalid value for enumeration:VV10` |
| Source `input_cp2k_xc.F` NON_LOCAL `TYPE` enum | only `DRSLL`, `LMKLL`, `RVV10` |
| `input_constants.F` | `vdw_nl_DRSLL=1`, `vdw_nl_LMKLL=2`, `vdw_nl_RVV10=3` — no `vdw_nl_VV10` |
| `data/` | ships `rVV10_kernel_table.dat` + `vdW_kernel_table.dat` (DRSLL/LMKLL); **no** `VV10_kernel_table.dat` |
| Wave-6 re-probe (`cp2k.psmp` 2026.1, 2026-07-26) | same **ABORT** `invalid value for enumeration:VV10` |
| Wave-7 re-probe (`cp2k.psmp` 2026.1, 2026-07-26) | same **ABORT** before SCF; source enum still DRSLL\|LMKLL\|RVV10 only |
| Wave-8 re-probe (`cp2k.psmp` 2026.1, 2026-07-26) | same **ABORT** before SCF; source enum + data kernels unchanged |

Probe artefact: `tools/cp2k_oracle/inputs/gpw-vv10-ice/` (this tree).

## Relation to product rVV10

| | Stock VV10 (this slug) | rVV10 (product) |
| - | ---------------------- | ---------------- |
| Literature | Vydrov & Van Voorhis, *J. Chem. Phys.* **133**, 244103 (2010) | Sabatini, Gorni, de Gironcoli, *Phys. Rev. B* **87**, 041108(R) (2013) |
| CP2K surface | **none** | `TYPE RVV10`, kernel `rVV10_kernel_table.dat`, default b=6.3, C=0.0093 |
| ATV | not implemented | `GpwDft::rvv10` / `GPWEvaluator::rvv10`; `reg_gpw-rvv10-ice` |

Do **not** rebrand `rvv10` with b=5.9 as stock VV10 — that would be false
capability closure. The RPS kernel table and (in general) the local-gap
revision differ; rVV10 is the only VV-family nonlocal CP2K can pin today.

## What would unblock

Any one of:

1. **Upstream CP2K** grows `TYPE VV10` + a stock VV10 kernel table (or libxc
   NLC path usable from Quickstep NON_LOCAL) — then offline pin + ATV RPS/kernel
   embed as done for rVV10.
2. **Loop policy change** allowing a non-CP2K offline oracle for this row
   (e.g. molecular libxc/PySCF NLC on a fixed density, or a self-generated
   Vydrov-φ RPS table with dual-frozen ATV gate + literature formula unit
   tests). That is **outside** the current CP2K-oracle contract for this loop.
3. **Supersede** the queue row as not-a-CP2K-gap: roadmap `dft.md` §8 VV10
   stays ○ for “stock VV10”; product nonlocal VV-family on GPW is rVV10 ✓.
   Harness should move the queue row from `## in_progress` to `## blocked`
   (oracle-surface gap, same class as other non-twinable rows).

## Acceptance (when unblocked)

- Public flag `GpwDft::vv10` (and evaluator mirror), exclusive vs
  `d3`/`d3bj`/`d4`/`ts`/`mbd`/`vdw_df`/`vdw_df2`/`rvv10`.
- No public `*Config` / `make_*` / `build_*`.
- Post-SCF E_c^{nl} on final total density → `GpwResult::e_disp` [Ha], folded
  into `e_total`; Frame energy still eV on the public surface.
- `regressions/gpw/gpw-vv10-ice/` with offline pin + REFERENCE.md.
- Unit tests for q0 / kernel piece (mirroring `test_cpu_gpw_rvv10`).
- Dual gate: correctness pin; `perf_exempt` allowed (E_disp pin policy as
  rVV10/vdW-DF siblings) **or** \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\).

## Non-goals (this slug)

- Self-consistent nonlocal potential / forces (post-SCF energy only, same as
  product rVV10/vdW-DF).
- Molecular post-HF (PySCF loop).
- Claiming CP2K TYPE VV10 exists.
- Implementing `GpwDft::vv10` without an offline pin that satisfies the loop
  oracle contract.

## Decision

**`acceptance_met = false`.** Block reason: no CP2K oracle surface for stock
VV10; cannot fulfil offline-CP2K pin contract. Product nonlocal VV-family
capability on GPW remains **rVV10**. Wave-8 worker re-probed live
(`cp2k.psmp` ABORT `invalid value for enumeration:VV10`), re-read
`input_cp2k_xc.F` / `input_constants.F` / `data/` kernels, rebased onto
current `dev` (`be75a5ed`); no `GpwDft::vv10` impl landed (would be false
closure). Harness should move the queue row from `## in_progress` to
`## blocked` (stop re-staffing until unblock checklist clears).

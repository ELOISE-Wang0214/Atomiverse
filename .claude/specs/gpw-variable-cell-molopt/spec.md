# Spec — `gpw-variable-cell-molopt`

**Status:** acceptance met (dual gate green)  
**Loop:** cp2k-method-loop  
**Branch:** `wf/gpw-variable-cell-molopt`  
**Dep:** `gpw-variable-cell-ice` ✓

## Intent

MOLOPT re-twin of GPW CELL_OPT on the ice-like water chain. Product machinery
unchanged (`AnisotropicCellRelax` + `StressConverged` + analytic fixed-\(P\)
stress). Basis pivot: **DZVP-MOLOPT-SR-GTH** (full `DZVP-MOLOPT-GTH` periodic
collocation > 60 s/SCF — secondary CP2K pin only).

## Twin

| Knob | Value |
| ---- | ----- |
| Basis | DZVP-MOLOPT-SR-GTH |
| XC / PP | PBE / GTH-PBE-q6/q1 |
| a₀ | 3.24 Å; b = c = 8 Å fixed |
| a\* pin | 3.200237 Å |
| E\* pin | −17.184075549824836 Ha |

## Acceptance

| Gate | Result |
| ---- | ------ |
| \|Δa\| ≤ 0.02 Å | 0.0012 Å |
| \|ΔE\| ≤ 1e-3 Ha | 1.4e-4 Ha |
| wall ≤ 60 s | 42.5 s |
| T_atv/T_cp2k ≤ 10 | 7.5 |

# Spec — second-gen-cpmd-kpoints

**Slug:** `second-gen-cpmd-kpoints`  
**Branch:** `wf/second-gen-cpmd-kpoints`  
**Queue:** CP2K method loop (wave 12)  
**Domain:** dynamics / GPW AIMD  
**Parent:** `dev` @ assign wave 12

## Goal

Close multi-k **second-generation Car–Parrinello MD** (Kühne SGCP): Kolafa
**ASPC on complex Bloch densities \(P_k\)** + partial SCF corrector on the
GPW multi-k path.

Γ RKS/UKS SGCP (real \(P\) / dual \(P^\sigma\)) and OT correctors are
already product. Multi-k XL-BOMD already seeds `guess_P_k` via
`XlBomdCplx`. This slug wires the **same multi-k density-guess socket** to
ASPC extrapolation (not Niklasson Verlet).

## Non-goals

- Classical CPMD with fictitious mass on Bloch \(C_k\) (`cpmd-kpoints`)
- Zero-SCF linearized XL forces (`xl-bomd-linearized-forces`)
- Multi-k SGCP + OT corrector (solver rejects density guess + `use_ot`)
- Live CP2K in ctest; fair wall only if offline pin exists

## Design

### Public API (members only; no `*Config`)

```cpp
cpu::gpw::GPWEvaluator eval({H}, /*cutoff*/100, /*max_iter*/60, /*tol*/1e-5,
                            /*mix*/0.3, XcKind::PADE, /*nonperiodic*/false,
                            /*compute_forces*/true);
eval.k_mesh[0] = 1; eval.k_mesh[1] = 1; eval.k_mesh[2] = 2;
eval.k_time_reversal = false;
eval.second_gen_cp = true;
eval.second_gen_cp_k = 1;
eval.second_gen_cp_max_iter = 5;
```

### Internals

1. **`scf::ASPCPredictorCplx`** — same Kolafa \(\omega_k\) and corrector
   weights as `ASPCPredictor`, element-wise on `linalg::CplxMatrix`.
2. **`GPWEvaluator` multi-k SGCP branch:**
   - Cold: full `max_iter` multi-k SCF; record each \(P_k\) (and \(P_k^\beta\)
     if UKS) into one ASPC history per k.
   - Warm: predict \(\{P_k\}\), cap SCF at `second_gen_cp_max_iter`, pass
     `guess_P_k` (/`guess_P_k_beta`) into `GpwDft::run`, Kolafa-correct +
     record, forces via `gpw_forces_from_scf` with `converged=true` (forces
     at corrector state).
3. **Guards:** mutual exclusion with `cpmd` / `xl_bomd` unchanged; multi-k
   + `second_gen_cp` + `use_ot` throws; multi-k + real `guess_P` already
   throws inside the solver.

### Literature

- T. D. Kühne et al., *Phys. Rev. Lett.* **98**, 066401 (2007)
- J. Kolafa, *J. Comput. Chem.* **25**, 335 (2004)
- CP2K `EXTRAPOLATION ASPC` + `&KPOINTS`

## Acceptance (summary)

See `acceptance.md`. Dual gate: regression green + unit guards; perf_exempt
when no offline CP2K wall twin on host.

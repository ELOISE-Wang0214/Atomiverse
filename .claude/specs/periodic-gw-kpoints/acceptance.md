# Acceptance — periodic-gw-kpoints

## Dual gate

### Unit (`tests/cpu/gpw/test_cpu_gpw_g0w0.cpp`)

- Existing: dRPA amplitudes match energy path; G0W0 guards (before build /
  UKS / unconverged).
- New: multi-k kernel at \(N_k=1\) with real-valued complex MOs recovers
  the real Γ spectral path on the same orbitals (QP HOMO/LUMO within
  1e-8 Ha), **or** equivalent internal amplitude/Omega recovery if a
  full QP path is too heavy for unit tier.
- New: `GpwDft::g0w0` no longer documents / throws solely for multi-k
  (throw only when multi-k `C_k` missing or `mp2_ri=false`).

### Regression (`reg_periodic-gw-kpoints`)

System: diamond Si bct 4-atom, MP 2×1×1, PBE, DZVP-MOLOPT-GTH /
GTH-PBE-q4, 80 Ry, `n_virt=2` screening (same deck family as
`periodic-rpa-kpoints`).

Hard gates:

1. SCF converged; `n_kpts==2`; `C_k.size()==2`.
2. `g0w0` returns without throw; `gw.n_kpts==2`; `gw.ri==true`.
3. Newton QP converged (`gw.converged`).
4. Finite positive min-direct QP gap and KS gap; HOMO QP < LUMO QP on
   the reported spectrum.
5. First-green absolute band on `gap_ev` (setup integrity; freeze after
   first green; replace with fair CP2K when spectral twin exists).
6. Wall ≤ 60 s single process.

## Perf

`perf_exempt` until a fair wall twin exists (CP2K RI-G0W0 + minimax ≠
k-diagonal spectral + MT exchange). Target ratio ≤ 10 when fair twin
available. Never launch CP2K from ctest.

## Done when

- Spec acceptance met; regression green; unit green.
- `acceptance_met=true`; commit on `wf/periodic-gw-kpoints`.

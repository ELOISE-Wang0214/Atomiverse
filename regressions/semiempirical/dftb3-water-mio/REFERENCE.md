# H2O SCC-DFTB3 (static mio-style SK)

```text
oracle: dftb+
oracle_version: formula-Gaus2011 (DFTB+ module boundary; no host DFTB+ binary)
```

## Method

Self-consistent-charge **DFTB3** (third-order charge fluctuation):

- Gaus et al., *J. Chem. Theory Comput.* **7**, 931 (2011),
  DOI [10.1021/ct100684s](https://doi.org/10.1021/ct100684s)
- Second-order SCC backbone: Elstner et al., *Phys. Rev. B* **58**, 7260
  (1998), DOI [10.1103/PhysRevB.58.7260](https://doi.org/10.1103/PhysRevB.58.7260)

Family-locked numerical oracle package is **DFTB+** (never OpenMOPAC / xtb).

## Third-order convention (what the code computes)

With \(\gamma_{AB}\) the DFTB+ Elstner charge kernel and
\(\Gamma_{AB}=\partial\gamma_{AB}/\partial q_A\) (on-site \(\gamma_{AA}=U_A\),
\(\Gamma_{AA}=U^{\mathrm{d}}_A\); note \(\Gamma_{AB}\neq\Gamma_{BA}\)):

\[
E^{(2)}=\tfrac12\sum_{AB}\gamma_{AB}\,\Delta q_A\Delta q_B,\qquad
E^{(3)}=\tfrac13\sum_{AB}\Delta q_A^2\,\Delta q_B\,\Gamma_{AB}
\]

which is the Gaus/Cui/Elstner (2011) form. The SCC potential is its exact
derivative,

\[
V^{(3)}_C=\tfrac23\,\Delta q_C\sum_B\Gamma_{CB}\Delta q_B
         +\tfrac13\sum_B\Delta q_B^2\,\Gamma_{BC},
\]

so \(\sum_A\Delta q_A V_A=2E^{(2)}+3E^{(3)}\) and the band-structure
double-counting correction is \(E^{(2)}+2E^{(3)}\). Pinned by
`tests/cpu/semiempirical/test_dftb3_scc.cpp` (Euler identity, finite
differences, \(\Delta q\equiv 0\) reduction).

## Hamiltonian twin note

SK H⁰/S and repulsives for all ordered H/C/N/O pairs are transcribed mio-1-1
(`basis_data/dftb_mio_sk_tables.h`). Hubbard derivatives \(U^{\mathrm{d}}\)
remain sketches (not in the SKF). Absolute DFTB+ DFTB3 identity is still
deferred until \(U^{\mathrm{d}}\) are production parameters. This regression
pins the **ATV DFTB3 Hamiltonian twin** energy (real mio SK + sketch \(U^{\mathrm{d}}\))
at the stated geometry.

## Geometry

Water, slightly off-equilibrium (Å):

| Atom | x    | y    | z   |
|------|------|------|-----|
| O    | 0.00 | 0.00 | 0.0 |
| H    | 0.96 | 0.00 | 0.0 |
| H    |−0.24 | 0.93 | 0.0 |

Cubic box 20 Å (isolated).

## Parameters

- `scc_tol = 1e-10`, `max_scc = 50`, forces off (energy gate)
- Units: energy **eV**, positions **Å**

## Pin (ATV twin after mio-1-1 SK; \(U^{\mathrm{d}}\) still sketch)

| Quantity | Value | Unit | Tol |
| -------- | ----- | ---- | --- |
| \(E_{\mathrm{tot}}\) | **−111.190697917671** | eV | \(10^{-6}\) |

Re-pinned 2026-07-29 after full mio-1-1 H/C/N/O SK transcription. DFTB2 at
the same geometry is **−110.996433406106 eV** (vs DFTB+ **−110.9590 eV**).
\(E^{(3)}\) still uses sketch \(U^{\mathrm{d}}\).

**What invalidated the previous pin (−107.618793920072 eV).** Three defects,
all fixed 2026-07-28:

1. **DFTB3 was reading the SKETCH H–H row** while DFTB0/DFTB2 read the
   transcribed mio-1-1 grid, so the two methods were not even solving the same
   Hamiltonian. `dftb3.cpp` now reads `internal::MioHHTable`, the one shared
   reader (`dftb2.cpp`'s private copy of it was deleted in the same wave).
2. **The SCC did not converge.** The undamped Picard iteration was a period-2
   limit cycle (residual ≈3.6 e); `scc_tol` was never met, the loop always
   exhausted `max_scc`, and the reported energy was selected by the **parity**
   of the budget (`max_scc = 50` → −107.618793920072 eV;
   `max_scc = 49 / 51` → −337.034306211295 eV). The old pin was one branch of
   that cycle — an artefact of the iteration budget, not a physical energy.
   Charges are now Anderson-mixed (`internal::ChargeMixer`, public knob
   `mix_beta`, default 0.3) and an SCC that exhausts `max_scc` **throws**
   instead of returning its last iterate.
3. **The double-counting route was wrong at second order.**
   `E = E_band − E_dc` assumed \(\sum_A V_A q_A = \sum_A V_A \Delta q_A\),
   but the exact identity is
   \(E_{\mathrm{band}} = \mathrm{Tr}(\rho H^0) + \sum_A V_A q_A\) with the
   **gross** Mulliken population, so the route carried a spurious
   \(\sum_A V_A q^0_A\) that is *linear* in \(\Delta q\). Measured at the
   converged \(\Delta q\) of the DFTB2 twin on this geometry:

   | Quantity | Value [eV] |
   | -------- | ---------- |
   | \(\mathrm{Tr}(\rho H) - \mathrm{Tr}(\rho H^0)\) | 8.114966333050 |
   | \(\sum_A V_A q_A\) (gross) | 8.114966333050 |
   | \(\sum_A V_A \Delta q_A = 2E^{(2)}\) | 0.260077801279 |
   | **spurious** \(\sum_A V_A q^0_A\) | **7.854888531771** |

   Both methods now assemble
   \(E = \mathrm{Tr}(\rho H^0) + E^{(2)} + E^{(3)} + E_{\mathrm{rep}}\)
   (the 3ob path already used this route and its DFTB+ pin did not move).

**DFTB+ oracle still owed.** No DFTB+ binary exists on this host (`module
spider dftb` → nothing), so no oracle number could be generated; nothing here
claims DFTB+ identity. Replace this pin with a DFTB+ 24.1 + mio-1-1 value once
the offline oracle infrastructure lands under `tools/tb_oracle/inputs/dftbplus/`.

## Reduction check (a real oracle, not a snapshot)

For a homonuclear system the Mulliken fluctuations vanish by symmetry
(\(\Delta q\equiv 0\)), so \(E^{(2)}=E^{(3)}=0\) identically and DFTB3 must
reduce **exactly** to DFTB2/DFTB0. Verified on H2 at \(r=0.74\) Å: DFTB3 and
DFTB2 now agree at **−18.366350522958 eV**, which *is* the DFTB+ 24.1 +
mio-1-1 pin **−18.36635052 eV**
(`regressions/semiempirical/dftb2-dftbplus-energy-pin/`, residual 3.0e-09 eV).
Before the table rewire DFTB3 reported −18.137397986146 eV here — a 0.229 eV
gap that was pure sketch-vs-mio parameter mismatch. Gated by
`tests/cpu/semiempirical/test_dftb3.cpp::ReducesToDftb2ForHomonuclearH2`.

## Gates

- Finite, non-zero total energy [eV] on `frame.meta[kEnergy]`
- SCC residual < `scc_tol` (converged, not budget-exhausted)
- Hard-coded energy pin vs table above (tol \(10^{-6}\) eV)
- \(|E_{\mathrm{DFTB3}} - E_{\mathrm{DFTB2}}| > 10^{-6}\) eV (third-order active)
- \(0.5\,|E^{(3)}| \le |E_{\mathrm{DFTB3}} - E_{\mathrm{DFTB2}}| \le
  2\,|E^{(3)}|\) — the third-order correction must move the energy by its
  **own** magnitude

### Closed: the DFTB3−DFTB2 magnitude is now commensurate with \(E^{(3)}\)

At this geometry \(|E_{\mathrm{DFTB3}}-E_{\mathrm{DFTB2}}| =
2.179184\times10^{-3}\) eV against \(E^{(3)} = -2.240219\times10^{-3}\) eV
(ratio **0.973**). Before the fixes the same comparison read 69.07 eV against
\(E^{(3)} = 5.01\) eV — the two methods were sampled from two *different*
limit cycles and both carried the spurious \(\sum_A V_A q^0_A\) term.

## Perf

`perf_exempt: true` — host has no DFTB+ binary for fair wall ratio; formula
oracle only. Revisit when offline DFTB+ pin infrastructure lands under
`tools/tb_oracle/inputs/dftbplus/`.

# H2O SCC-DFTB2 (static mio-style SK)

```text
oracle: literature
oracle_version: Elstner et al. PRB 58, 7260 (1998) — framework gate (not DFTB+ E pin)
```

## Method

Self-consistent-charge DFTB2 (Elstner et al., *Phys. Rev. B* **58**, 7260
(1998), DOI [10.1103/PhysRevB.58.7260](https://doi.org/10.1103/PhysRevB.58.7260)).

## Data

In-tree static tables for H/C/N/O: H onsite/U + H–H from transcribed mio-1-1
(`dftb_mio_hh_sk.h`); other pairs exponential sketches (`dftb_mio_tables.h`).
No runtime `.skf` parser. Full multi-element mio E identity is a separate pin
(`dftb2-dftbplus-energy-pin` for H₂).

## Geometry

Water, slightly off-equilibrium (Å):

| Atom | x    | y    | z   |
|------|------|------|-----|
| O    | 0.00 | 0.00 | 0.0 |
| H    | 0.96 | 0.00 | 0.0 |
| H    |−0.24 | 0.93 | 0.0 |

Cubic box 20 Å (isolated).

## Gates

- Finite, non-zero total energy [eV] on `frame.meta[kEnergy]`
- Finite analytic forces [eV/Å] with |F|² > 0
- One Cartesian component matches central FD to **1×10⁻⁵ eV/Å** (measured
  4.5e-08 eV/Å). A converged SCC-DFTB2 gradient is stationary, so the analytic
  force is FD-exact; the old 1.5 eV/Å window was masking a real defect.

## The reported energy moved (2026-07-28)

`E = −42.11915111 eV` → **−94.27162200 eV** at this geometry. The old value was
a snapshot of a **non-convergent Picard limit cycle** — `scc_tol` was never
met, the loop always burned `max_scc`, and which branch of the period-2 cycle
got reported was decided by the parity of the budget — taken over a
double-counting route that was wrong at second order (`E = E_band − E²`
carries a spurious `Σ_A V_A q⁰_A`, **7.85 eV** here, because the band energy
contracts `V_A` with the *gross* Mulliken population).

The new value is a **converged** SCC solution: residual < `scc_tol` after ~6
Anderson-mixed iterations, identical for `max_scc` = 49 / 50 / 51, on the route
`E = Tr(ρH⁰) + E² + E_rep`. This case has **no absolute energy pin** and none is
added: it is not an oracle number. DFTB+ oracle identity for the MioSketch
multi-element path is still **owed** — no DFTB+ binary exists on this host
(`module spider dftb` → nothing). The one DFTB+-pinned MioSketch number,
H₂ at 0.74 Å (`dftb2-dftbplus-energy-pin`, −18.36635052 eV), is **unmoved**:
H₂ has Δq ≡ 0, so neither the mixing nor the double-counting change can touch
it.
